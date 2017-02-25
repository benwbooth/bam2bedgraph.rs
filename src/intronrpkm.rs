#![feature(ordering_chaining)]
#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::str;
use std::vec::Vec;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use std::ops::Range;

extern crate regex;
use regex::Regex;

extern crate rust_htslib;
use rust_htslib::bam::Read;
use rust_htslib::bam::IndexedReader;

extern crate bio;
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;
use bio::io::gff;
use bio::io::gff::GffType;

extern crate structopt;
#[macro_use]
extern crate structopt_derive;
use structopt::StructOpt;

extern crate bam2bedgraph;
use bam2bedgraph::errors::*;
use bam2bedgraph::cigar2exons;

#[derive(StructOpt, Debug)]
#[structopt(name = "intronrpkm", about = "Analyze RPKM values in intronic space")]
struct Options {
    #[structopt(help = "A genome annotation file in gtf/gff format", name="ANNOT_FILE")]
    annotfile: String,
    #[structopt(long="exon_type", help = "The exon type(s) to search for", name="TYPE")]
    exon_type: Vec<String>,
    #[structopt(long="transcript_type", help = "The transcript type(s) to search for", name="TYPE")]
    transcript_type: Vec<String>,
    #[structopt(long="gene_type", help = "The gene type(s) to search for", name="TYPE")]
    gene_type: Vec<String>,
    #[structopt(long="use_annotated_cassettes", help = "Should the annotated cassettes be used in reannotating?")]
    use_annotated_cassettes: bool,
    #[structopt(long="compare_constituitive", help = "Should constituitve exons be compared?")]
    compare_constituitive: bool,
    #[structopt(long="compare_cassette", help = "Should cassette exons be compared?")]
    compare_cassette: bool,
    #[structopt(help = "The set of .bam files to analyze", name="BAMFILE")]
    bamfiles: Vec<String>,
}

#[derive(Default)]
pub struct Record {
    seqname: String,
    source: String,
    feature_type: String,
    start: u64,
    end: u64,
    score: String,
    strand: String,
    frame: String,
    attributes: HashMap<String, String>,
}

impl PartialEq for Record {
    fn eq(&self, other: &Self) -> bool {
        self.seqname == other.seqname && self.source == other.source &&
        self.feature_type == other.feature_type && self.start == other.start &&
        self.end == other.end && self.score == other.score &&
        self.strand == other.strand && self.frame == other.frame
    }
}
impl Eq for Record {}

impl Hash for Record {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seqname.hash(state);
        self.source.hash(state);
        self.feature_type.hash(state);
        self.start.hash(state);
        self.end.hash(state);
        self.score.hash(state);
        self.strand.hash(state);
        self.frame.hash(state);
    }
}

impl Record {
    fn new() -> Record {
        Record { ..Default::default() }
    }
    fn from_record(record: &gff::Record) -> Record {
        let strand = match record.strand() {
            Some(bio::io::Strand::Forward) => "+",
            Some(bio::io::Strand::Reverse) => "-",
            _ => ".",
        };
        let score = match record.score() {
            Some(num) => num.to_string(),
            _ => String::from("."),
        };
        Record {
            seqname: String::from(record.seqname()),
            source: String::from(record.source()),
            feature_type: String::from(record.feature_type()),
            start: *record.start(),
            end: *record.end(),
            score: score,
            strand: String::from(strand),
            frame: String::from(record.frame()),
            attributes: record.attributes().clone(),
        }
    }
}

struct IndexedAnnotation {
    rows: Vec<Record>,
    id2row: HashMap<String, usize>,
    row2parents: HashMap<usize, Vec<usize>>,
    row2children: HashMap<usize, Vec<usize>>,
    tree: HashMap<String, IntervalTree<u64, usize>>,
}
impl IndexedAnnotation {
    fn from_file(annotfile: &str, gene_type: &str, transcript_type: &str) -> Result<IndexedAnnotation> {
        // figure out the file type of the annotation file
        let caps = Regex::new(r"\.([^.]*)$")?.captures(annotfile).r()?;
        let ext = caps.get(1).map_or("", |m| m.as_str()).to_lowercase();
        let filetype = match ext.as_ref() {
            "gff" | "gff3" => Ok(GffType::GFF3),
            "gtf" | "gtf2" => Ok(GffType::GTF2),
            ext => {
                Err(format!("Don't know how to parse file {} with extension {}!",
                            annotfile,
                            ext))
            }
        }?;
        // read all the records into a single vector
        let mut reader = gff::Reader::from_file(&annotfile, filetype)?;
        let mut id2row = HashMap::<String, usize>::new();
        let mut rows = Vec::<Record>::new();
        for r in reader.records().enumerate() {
            let row = r.0;
            let record = r.1?;
            if let Some(id) = record.attributes().get("ID") {
                id2row.insert(id.clone(), row);
            }
            rows.push(Record::from_record(&record));
        }
        // populate row2parents and row2children indices
        let mut row2parents = HashMap::<usize, BTreeSet<usize>>::new();
        let mut row2children = HashMap::<usize, BTreeSet<usize>>::new();
        // make fake rows for GTF gene_id and transcript_id attributes
        let mut fake_rows = Vec::<Record>::new();
        match filetype {
            GffType::GTF2 => {
                let mut firstrow = HashMap::<&Record, usize>::new();
                for (row, record) in rows.iter().enumerate() {
                    // in GTF, if two rows both have identical fields (except attributes), they
                    // are the same feature.
                    let row = *firstrow.entry(record).or_insert(row);
                    if let Some(transcript_id) = record.attributes.get("transcript_id") {
                        let mut transcript_record = Record::new();
                        transcript_record.feature_type = transcript_type.to_string();
                        transcript_record.attributes
                            .insert("ID".to_string(), transcript_id.clone());
                        if let Some(transcript_name) = record.attributes.get("transcript_name") {
                            transcript_record.attributes
                                .insert("Name".to_string(), transcript_name.clone());
                        }
                        let transcriptrow = rows.len() + fake_rows.len();
                        id2row.insert(transcript_id.clone(), transcriptrow);
                        row2parents.entry(row).or_insert_with(BTreeSet::new).insert(transcriptrow);
                        row2children.entry(transcriptrow)
                            .or_insert_with(BTreeSet::new)
                            .insert(row);
                        fake_rows.push(transcript_record);

                        if let Some(gene_id) = record.attributes.get("gene_id") {
                            let mut gene_record = Record::new();
                            gene_record.feature_type = gene_type.to_string();
                            gene_record.attributes.insert("ID".to_string(), gene_id.clone());
                            if let Some(gene_name) = record.attributes.get("gene_name") {
                                gene_record.attributes
                                    .insert("Name".to_string(), gene_name.clone());
                            }
                            let generow = rows.len() + fake_rows.len();
                            id2row.insert(gene_id.clone(), generow);
                            row2parents.entry(transcriptrow)
                                .or_insert_with(BTreeSet::new)
                                .insert(generow);
                            row2children.entry(generow)
                                .or_insert_with(BTreeSet::new)
                                .insert(transcriptrow);
                            fake_rows.push(gene_record);
                        }
                    }
                }
            }
            GffType::GFF3 => {
                for (row, record) in rows.iter().enumerate() {
                    if let Some(parent) = record.attributes.get("Parent") {
                        for p in parent.split(',') {
                            if let Some(parentrow) = id2row.get(p) {
                                row2parents.entry(row)
                                    .or_insert_with(BTreeSet::new)
                                    .insert(*parentrow);
                                row2children.entry(*parentrow)
                                    .or_insert_with(BTreeSet::new)
                                    .insert(row);
                            }
                        }
                    }
                }
            }
            _ => (),
        };
        rows.append(&mut fake_rows);

        // fill in missing information, build interval tree
        let mut updates = Vec::<(usize, Record)>::new();
        for (row, record) in rows.iter().enumerate() {
            // if the seqname is missing then the record is not complete
            if record.seqname.is_empty() {
                // recursively get all child nodes
                let mut children = HashSet::<usize>::new();
                let mut addchildren = BTreeSet::<usize>::new();
                if let Some(cs) = row2children.get_mut(&row) {
                    addchildren.append(cs);
                }
                while !addchildren.is_empty() {
                    let mut morechildren = BTreeSet::<usize>::new();
                    for ac in &mut addchildren.iter() {
                        if !children.contains(ac) {
                            children.insert(*ac);
                            if let Some(cs) = row2children.get_mut(ac) {
                                morechildren.append(cs);
                            }
                        }
                    }
                    addchildren.clear();
                    addchildren.append(&mut morechildren);
                }
                // aggregate the child nodes to get the missing information
                let mut seqname = HashMap::<&str, u64>::new();
                let mut source = HashMap::<&str, u64>::new();
                let mut start: Option<u64> = Option::None;
                let mut end: Option<u64> = Option::None;
                let mut strand = HashMap::<&str, u64>::new();
                for c in &children {
                    let child = &rows[*c];
                    if start.is_none() || child.start < start.r()? {
                        start = Some(child.start);
                    }
                    if end.is_none() || child.end > end.r()? {
                        end = Some(child.end);
                    }
                    *seqname.entry(&child.seqname).or_insert(0u64) += 1;
                    *source.entry(&child.source).or_insert(0u64) += 1;
                    *strand.entry(&child.strand).or_insert(0u64) += 1;
                }
                // enqueue the missing information
                let mut update = Record::new();
                if let Some((seqname, _)) = seqname.iter().max_by(|a, b| a.1.cmp(b.1)) {
                    update.seqname = String::from(*seqname);
                }
                if let Some((source, _)) = source.iter().max_by(|a, b| a.1.cmp(b.1)) {
                    update.source = String::from(*source);
                }
                if let Some(start) = start {
                    update.start = start;
                }
                if let Some(end) = end {
                    update.end = end;
                }
                if let Some((strand, _)) = strand.iter().max_by(|a, b| a.1.cmp(b.1)) {
                    update.strand = String::from(*strand);
                }
                updates.push((row, update));
            }
        }
        // update the missing information
        for u in &mut updates {
            let update = &mut u.1;
            let mut record = &mut rows[u.0];
            record.seqname = update.seqname.to_string();
            record.source = update.source.to_string();
            record.start = update.start;
            record.end = update.end;
            record.strand = update.strand.to_string();
        }
        // sort row2parents by coordinates
        let mut row2parents_update = HashMap::<usize, Vec<usize>>::new();
        for (row, parents) in row2parents {
            let mut up = row2parents_update.entry(row).or_insert_with(Vec::new);
            let mut vec: Vec<_> = parents.into_iter().collect();
            vec.sort_by(|a, b| {
                rows[*a]
                    .seqname
                    .cmp(&rows[*b].seqname)
                    .then_with(|| rows[*a].start.cmp(&rows[*b].start))
            });
            up.append(&mut vec);
        }
        // sort row2children by coordinates
        let mut row2children_update = HashMap::<usize, Vec<usize>>::new();
        for (row, children) in row2children {
            let mut up = row2children_update.entry(row).or_insert_with(Vec::new);
            let mut vec: Vec<_> = children.into_iter().collect();
            vec.sort_by(|a, b| {
                rows[*a]
                    .seqname
                    .cmp(&rows[*b].seqname)
                    .then_with(|| rows[*a].start.cmp(&rows[*b].start))
            });
            up.append(&mut vec);
        }
        // build the interval tree
        let mut tree = HashMap::<String, IntervalTree<u64, usize>>::new();
        for (row, record) in rows.iter().enumerate() {
            let mut t = tree.entry(record.seqname.clone()).or_insert_with(IntervalTree::new);
            t.insert(Interval::new((record.start - 1)..(record.end))?, row);
        }
        Ok(IndexedAnnotation {
            rows: rows,
            id2row: id2row,
            row2children: row2children_update,
            row2parents: row2parents_update,
            tree: tree,
        })
    }
}

#[derive(PartialEq, Eq, Hash)]
struct Cassette {
    range: Range<u64>,
    cassette_row: Option<usize>,
}

#[derive(PartialEq, Eq, Hash)]
struct ConstituitivePair {
    exon1_row: usize,
    exon2_row: usize,
    // start..end, optional cassette_row
    cassettes: Vec<Cassette>,
}

fn find_constituitive_exons(annot: &IndexedAnnotation,
                            options: &Options)
                            -> Result<Vec<ConstituitivePair>> {
    let gene_types: HashSet<_> = options.gene_type.clone().into_iter().collect();
    let transcript_types: HashSet<_> = options.transcript_type.clone().into_iter().collect();
    let exon_types: HashSet<_> = options.exon_type.clone().into_iter().collect();
    
    // set default feature types
    let mut exonpairs = Vec::<ConstituitivePair>::new();
    for (gene_row, gene) in annot.rows.iter().enumerate() {
        if gene_types.contains(&gene.feature_type) {
            if let Some(transcript_rows) = annot.row2children.get(&gene_row) {
                // (exon1_row, exon2_row) -> (cassette_rows, transcript_rows)
                let mut constituitive_pairs = HashMap::<(usize,usize),(HashSet<usize>,HashSet<usize>)>::new();
                // look in each transcript
                'transcript: for transcript_row in transcript_rows {
                    let transcript = &annot.rows[*transcript_row];
                    // make sure transcript coordinates are consistent with the gene coordinates
                    if transcript_types.contains(&transcript.feature_type) &&
                       transcript.seqname == gene.seqname &&
                       transcript.strand == gene.strand {
                        if let Some(feature_rows) = annot.row2children.get(transcript_row) {
                            // get the set of exon features we care about
                            let mut exon_rows = Vec::<usize>::new();
                            for feature_row in feature_rows {
                                let feature = &annot.rows[*feature_row];
                                if exon_types.contains(&feature.feature_type) {
                                    // skip any transcripts with weird exons
                                    if feature.seqname != transcript.seqname || feature.strand != transcript.strand {
                                        exon_rows.clear();
                                        continue 'transcript;
                                    }
                                    exon_rows.push(*feature_row);
                                }
                            }
                            // get the constituitive exons
                            let mut constituitive = HashSet::<usize>::new();
                            'exon: for exon_row in &exon_rows {
                                if let Some(parent_rows) = annot.row2parents.get(exon_row) {
                                    let parents: HashSet<_> = parent_rows.into_iter().collect();
                                    // are there any transcripts in this gene that this exon
                                    // does NOT have as a parent?
                                    for transcript_row in transcript_rows {
                                        if !parents.contains(transcript_row) {
                                            continue 'exon;
                                        }
                                    }
                                    constituitive.insert(*exon_row);
                                }
                            }
                            // get pairs of constituitive exons
                            for (rank, exon_row) in exon_rows.iter().enumerate() {
                                if rank < exon_rows.len()-1 && constituitive.contains(exon_row) {
                                    // two adjacent constitutive exons
                                    let adjacent_row = exon_rows[rank+1];
                                    if constituitive.contains(&adjacent_row) {
                                        let cp = constituitive_pairs.entry(
                                            (*exon_row, adjacent_row)).
                                            or_insert_with(|| (HashSet::new(), HashSet::new()));
                                        cp.1.insert(*transcript_row);
                                    }
                                    // two constituitive exons with a single cassette inbetween
                                    else if rank < exon_rows.len()-2 {
                                        let next_row = exon_rows[rank+2];
                                        if constituitive.contains(&next_row) {
                                            let cp = constituitive_pairs.entry(
                                                (*exon_row, next_row)).
                                                or_insert_with(|| (HashSet::new(), HashSet::new()));
                                            cp.0.insert(adjacent_row);
                                            cp.1.insert(*transcript_row);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // filter out constituitive pairs that do not cover every transcript
                let mut valid_constituitive_pairs: Vec<_> = constituitive_pairs.into_iter().
                    filter(|p| (p.1).1.len() == transcript_rows.len()).
                    map(|p| ConstituitivePair {
                        exon1_row: (p.0).0,
                        exon2_row: (p.0).1,
                        cassettes: (p.1).0.into_iter().
                            map(|c| Cassette {
                                range: annot.rows[c].start-1..annot.rows[c].end, 
                                cassette_row: Some(c),
                            }).collect(),
                    }).
                    collect();
                exonpairs.append(&mut valid_constituitive_pairs);
            }
        }
    }
    Ok(exonpairs)
}

fn reannotate_regions(
    annot: &IndexedAnnotation,
    pairs: &[ConstituitivePair], 
    bamfiles: &[String], 
    use_annotated_cassettes: bool) 
    -> Result<Vec<ConstituitivePair>>
{
    let reannotated = Vec::<ConstituitivePair>::new();
    for pair in pairs {
        for bamfile in bamfiles {
            let mut bam = IndexedReader::from_path(bamfile)?;
            let chr = &annot.rows[pair.exon1_row].seqname;
            if let Some(tid) = bam.header.tid(&chr.clone().into_bytes()) {
                let start = &annot.rows[pair.exon1_row].start;
                let end = &annot.rows[pair.exon2_row].end;
                bam.seek(tid, *start as u32, *end as u32)?;
                for (i, record) in bam.records().enumerate() {
                }
            }
        }
    }
    
    Ok(reannotated)
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    // set defaults for feature types
    options.gene_type =
        (if options.gene_type.is_empty() { vec!["gene".to_string()] } 
         else { options.gene_type.clone() }).into_iter().collect();
    options.transcript_type = 
        (if options.transcript_type.is_empty() { vec!["mRNA".to_string()] } 
         else { options.transcript_type.clone()
        }).into_iter().collect();
    options.exon_type =
        (if options.exon_type.is_empty() { vec!["CDS".to_string()] }
         else { options.exon_type.clone() }).into_iter().collect();

    let annot = IndexedAnnotation::from_file(&options.annotfile, 
        options.gene_type.get(0).r()?, 
        options.transcript_type.get(0).r()?)?;
    let exonpairs = find_constituitive_exons(&annot, &options)?;
    let reannotated = reannotate_regions(
        &annot,
        &exonpairs, 
        &options.bamfiles, 
        options.use_annotated_cassettes);

    Ok(())
}

fn main() {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "1");

    if let Err(ref e) = run() {
        println!("error: {}", e);
        for e in e.iter().skip(1) {
            println!("caused by: {}", e);
        }
        if let Some(backtrace) = e.backtrace() {
            println!("backtrace: {:?}", backtrace);
        }
        ::std::process::exit(1);
    }
}
