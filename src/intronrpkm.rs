#![feature(plugin)]
#![plugin(indoc)]
#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::str;
use std::vec::Vec;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use std::ops::Range;
use std::process::Command;
use std::cmp::Ordering::Less;
use std::fs::File;
use std::io::{BufReader, BufWriter, BufRead, Write};
use std::io::{stdout, stderr, sink};
use std::path::Path;

#[macro_use] 
extern crate lazy_static;

#[macro_use]
extern crate url;
use url::percent_encoding::{percent_decode, utf8_percent_encode, SIMPLE_ENCODE_SET};

extern crate regex;
use regex::Regex;
use regex::Captures;

extern crate rust_htslib;
use rust_htslib::bam::Read;
use rust_htslib::bam::IndexedReader;

extern crate bio;
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;

extern crate structopt;
#[macro_use]
extern crate structopt_derive;
use structopt::StructOpt;

extern crate bam2bedgraph;
use bam2bedgraph::errors::*;
use bam2bedgraph::cigar2exons;

extern crate linked_hash_map;
use linked_hash_map::LinkedHashMap;

extern crate serde;
use serde::ser::SerializeStruct;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

extern crate itertools;
use itertools::Itertools;

extern crate sprs;
use sprs::{CsVec, CsVecOwned};

extern crate string_cache;
use string_cache::DefaultAtom as Atom;

#[derive(StructOpt, Debug)]
#[structopt(name = "intronrpkm", about = "Analyze RPKM values in intronic space")]
struct Options {
    // input files
    #[structopt(long="gff", help = "A genome annotation file in gff3 format", name="ANNOT_GFF_FILE")]
    annotfile_gff: Option<String>,
    #[structopt(long="gtf", help = "A genome annotation file in gtf format", name="ANNOT_GTF_FILE")]
    annotfile_gtf: Option<String>,
    #[structopt(long="bam1", short="1", help = "The set of stranded .bam files to analyze where read1 indicates strand", name="BAMFILE1")]
    bam1: Vec<String>,
    #[structopt(long="bam2", short="2", help = "The set of stranded .bam files to analyze where read2 indicates strand", name="BAMFILE2")]
    bam2: Vec<String>,
    #[structopt(long="bam", short="u", help = "The set of unstranded .bam files to analyze", name="BAMFILE")]
    bam: Vec<String>,
    #[structopt(long="chrmap", help = "Optional tab-delimited chr name mapping file", name="CHRMAP_FILE")]
    chrmap_file: Option<String>,
    #[structopt(long="vizchrmap", help = "Optional tab-delimited chr name mapping file for bigwig/bigbed exports", name="VIZCHRMAP_FILE")]
    vizchrmap_file: Option<String>,
    #[structopt(long="sizes", help = "Optional chr sizes file", name="SIZES_FILE")]
    sizes_file: Option<String>,
    // output file
    #[structopt(long="out", short="o", help = "Output file", name="OUT_FILE", default_value="-")]
    outfile: String,
    // feature types filter
    #[structopt(long="exon_type", help = "The exon type(s) to search for", name="EXON_TYPE")]
    exon_type: Vec<String>,
    #[structopt(long="transcript_type", help = "The transcript type(s) to search for", name="TRANSCRIPT_TYPE")]
    transcript_type: Vec<String>,
    #[structopt(long="gene_type", help = "The gene type(s) to search for", name="GENE_TYPE")]
    gene_type: Vec<String>,
    #[structopt(long="cds_type", help = "The CDS type(s) to search for", name="CDS_TYPE")]
    cds_type: Vec<String>,
    
    // flags
    #[structopt(long="max_slippage", help = "How many bases (or % of exon length) to search for start/stop histogram peaks", name="MAX_SLIPPAGE", default_value="10")]
    max_slippage: String,
    #[structopt(long="fill_incomplete_exons", help = "Should we try to fill in exons with only one splice junction?")]
    fill_incomplete_exons: bool,
    
    // debug output files
    #[structopt(long="debug", help = "Output all debug files?")]
    debug: bool,
    #[structopt(long="debug_annot_gff", help = "Write the input annotation as a gff debug file", name="DEBUG_ANNOT_GFF_FILE")]
    debug_annot_gff: Option<String>,
    #[structopt(long="debug_annot_gtf", help = "Write the input annotation as a gtf debug file", name="DEBUG_ANNOT_GTF_FILE")]
    debug_annot_gtf: Option<String>,
    #[structopt(long="debug_annot_bigbed", help = "Write the input annotation as a bigBed file", name="DEBUG_ANNOT_BIGBED_FILE")]
    debug_annot_bigbed: Option<String>,
    #[structopt(long="debug_annot_json", help = "Dump the indexed annotation as a debug test", name="DEBUG_ANNOT_JSON")]
    debug_annot_json: Option<String>,
    #[structopt(long="debug_exon_bigbed", help = "Output the constituitive pairs to a bigbed file", name="DEBUG_EXON_BIGBED")]
    debug_exon_bigbed: Option<String>,
    #[structopt(long="debug_bigwig", help = "Output the reannotation splice start/end sites as bigwig files", name="DEBUG_BIGWIG_FILE_PREFIX")]
    debug_bigwig: Option<String>,
    #[structopt(long="debug_reannot_bigbed", help = "Output the reannotated constituitive pairs as a bigbed file", name="DEBUG_REANNOT_BIGBED_FILE")]
    debug_reannot_bigbed: Option<String>,
    #[structopt(long="debug_rpkm_region_bigbed", help = "Output the rpkm region features as a bigbed file", name="DEBUG_RPKM_REGION_BIGBED_FILE")]
    debug_rpkm_region_bigbed: Option<String>,
    #[structopt(long="debug_mapped_reads_bigbed", help = "Output the rpkm mapped reads as a bigbed file", name="DEBUG_MAPPED_READS_BIGBED_FILE")]
    debug_mapped_reads_bigbed: Option<String>,
    #[structopt(long="debug_rpkmstats_json", help = "Dump the RpkmStats to a JSON file", name="DEBUG_RPKMSTATS_JSON")]
    debug_rpkmstats_json: Option<String>,
    #[structopt(long="debug_trackdb", help = "Write a UCSC trackDb.txt file with all the bigwigs/bigbeds", name="DEBUG_TRACKDB_FILE")]
    debug_trackdb: Option<String>,
}

#[derive(Default, Serialize, Deserialize)]
struct Record {
    row: usize,
    seqname: Atom,
    source: Atom,
    feature_type: Atom,
    start: u64,
    end: u64,
    score: Atom,
    strand: Atom,
    frame: Atom,
    attributes: LinkedHashMap<Atom, Atom>,
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
    fn from_row(row: usize, line: &str, filetype: &str, chrmap: &HashMap<Atom,Atom>) -> Result<Record> {
        lazy_static! {
            static ref COMMENT: Regex = Regex::new(r"^#").unwrap();
            static ref GTF_ATTR: Regex = Regex::new(r#"^(?P<key>\S+)\s+(?:"(?P<qval>[^"]*)"|(?P<val>\S+));\s*"#).unwrap();
        }
        if COMMENT.is_match(line) {
            return Err("Comment".into());
        }
        let fields: Vec<_> = line.split('\t').collect();
        let seqname = Atom::from(*fields.get(0).unwrap_or(&""));
        let seqname = chrmap.get(&seqname).unwrap_or(&seqname);
        Ok(Record {
            row: row,
            seqname: seqname.clone(),
            source: Atom::from(*fields.get(1).unwrap_or(&"")),
            feature_type: Atom::from(*fields.get(2).unwrap_or(&"")),
            start: fields.get(3).unwrap_or(&"0").parse::<u64>()?,
            end: fields.get(4).unwrap_or(&"0").parse::<u64>()?,
            score: Atom::from(*fields.get(5).unwrap_or(&"")),
            strand: Atom::from(*fields.get(6).unwrap_or(&"")),
            frame: Atom::from(*fields.get(7).unwrap_or(&"")),
            attributes: if filetype == "gff" {
                fields.get(8).unwrap_or(&"").split(';').map(|a| {
                    let kv: Vec<_> = a.splitn(2, '=').collect();
                    let key = Atom::from(percent_decode(kv[0].as_bytes()).decode_utf8().unwrap());
                    let value = Atom::from(percent_decode(kv.get(1).unwrap_or(&"").as_bytes()).decode_utf8().unwrap());
                    (key, value)
                }).collect()
            } else if filetype == "gtf" {
                GTF_ATTR.captures_iter(fields.get(8).unwrap_or(&"")).map(|caps| {
                    (Atom::from(caps["key"].to_string()), 
                     Atom::from(caps.name("qval").unwrap_or_else(|| caps.name("val").unwrap()).as_str().to_string()))
                }).collect()
            } else {
                return Err(format!("Don't know how to read filetype {}", filetype).into());
            },
        })
    }
    
    fn to_gtf(&self) -> Result<String> {
        Ok(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            self.seqname, 
            self.source, 
            self.feature_type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attributes.iter().map(|(k,v)| format!("{} \"{}\";", k, v)).join(" ")))
        
    }
    
    fn to_gff(&self) -> Result<String> {
        define_encode_set! {
            pub GFF_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'\t', '\r', '\n', ';', '%', '='}
        }
        Ok(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            self.seqname, 
            self.source, 
            self.feature_type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attributes.iter().map(|(k,v)| 
                format!("{}={}", 
                    utf8_percent_encode(k, GFF_ENCODE_SET), 
                    utf8_percent_encode(v, GFF_ENCODE_SET))).join(";")))
    }
}

struct IndexedAnnotation {
    rows: Vec<Record>,
    id2row: HashMap<Atom, usize>,
    row2parents: HashMap<usize, Vec<usize>>,
    row2children: HashMap<usize, Vec<usize>>,
    tree: HashMap<Atom, IntervalTree<u64, usize>>,
    chrmap: HashMap<Atom, Atom>,
    vizchrmap: HashMap<Atom, Atom>,
    refs: LinkedHashMap<Atom, u64>,
}
impl serde::Serialize for IndexedAnnotation
{
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
        where S: serde::Serializer
    {
        let mut ia = serializer.serialize_struct("IndexedAnnotation", 5)?;
        ia.serialize_field("rows", &self.rows)?;
        ia.serialize_field("id2row", &self.id2row)?;
        ia.serialize_field("row2parents", &self.row2parents)?;
        ia.serialize_field("row2children", &self.row2children)?;
        let mut tree = HashMap::<Atom, Vec<((u64,u64),usize)>>::new();
        for (chr,interval_tree) in &self.tree {
            let mut nodes = Vec::<((u64, u64),usize)>::new();
            for node in interval_tree.find(0..std::u64::MAX) {
                nodes.push(((node.interval().start,node.interval().end), *node.data()));
            }
            tree.insert(chr.clone(), nodes);
        }
        ia.serialize_field("tree", &tree)?;
        ia.serialize_field("chrmap", &self.chrmap)?;
        ia.serialize_field("vizchrmap", &self.chrmap)?;
        ia.serialize_field("refs", &self.refs)?;
        ia.end()
    }
}
impl IndexedAnnotation {
    fn from_gtf(
        annotfile: &str, 
        gene_type: &str, 
        transcript_type: &str,
        chrmap_file: &Option<String>,
        vizchrmap_file: &Option<String>) 
        -> Result<IndexedAnnotation> 
    {
        IndexedAnnotation::from_file(annotfile, "gtf", gene_type, transcript_type, chrmap_file, vizchrmap_file)
    }
    fn from_gff(
        annotfile: &str, 
        gene_type: &str, 
        transcript_type: &str, 
        chrmap_file: &Option<String>,
        vizchrmap_file: &Option<String>) 
        -> Result<IndexedAnnotation> 
    {
        IndexedAnnotation::from_file(annotfile, "gff", gene_type, transcript_type, chrmap_file, vizchrmap_file)
    }
    fn from_file(
        annotfile: &str, 
        filetype: &str, 
        gene_type: &str, 
        transcript_type: &str, 
        chrmap_file: &Option<String>,
        vizchrmap_file: &Option<String>) 
        -> Result<IndexedAnnotation> 
    {
        // read in the optional chr map
        let mut chrmap = HashMap::<Atom,Atom>::new();
        if let Some(charmap_file) = chrmap_file.clone() {
            let f = File::open(&charmap_file)?;
            let file = BufReader::new(&f);
            for line in file.lines() {
                let line = line?;
                let cols: Vec<&str> = line.split('\t').collect();
                if let Some(key) = cols.get(0) {
                    if let Some(value) = cols.get(1) {
                        chrmap.insert(Atom::from(*key), Atom::from(*value));
                    }
                }
            }
        }
        // read in the optional visualization chr map
        let mut vizchrmap = HashMap::<Atom,Atom>::new();
        if let Some(vizcharmap_file) = vizchrmap_file.clone() {
            let f = File::open(&vizcharmap_file)?;
            let file = BufReader::new(&f);
            for line in file.lines() {
                let line = line?;
                let cols: Vec<&str> = line.split('\t').collect();
                if let Some(key) = cols.get(0) {
                    if let Some(value) = cols.get(1) {
                        vizchrmap.insert(Atom::from(*key), Atom::from(*value));
                    }
                }
            }
        }
        
        let mut id2row = HashMap::<Atom, usize>::new();
        let mut rows = Vec::<Record>::new();
        let f = File::open(&annotfile)?;
        let file = BufReader::new(&f);
        let id_atom = Atom::from("ID");
        let name_atom = Atom::from("Name");
        let mut refs = HashMap::<Atom,u64>::new();
        for line in file.lines() {
            let row = rows.len();
            if let Ok(record) = Record::from_row(row, &line?, filetype, &chrmap) {
                if let Some(id) = record.attributes.get(&id_atom) {
                    id2row.insert(id.clone(), row);
                }
                // get the max ref lengths
                let mut reflength = refs.entry(record.seqname.clone()).or_insert(record.end);
                if *reflength < record.end { *reflength = record.end }
                rows.push(record);
            }
        }
        // sort the refs
        let mut sorted_refs = LinkedHashMap::<Atom,u64>::new();
        let chrs = refs.keys().collect::<Vec<_>>();
        for chr in chrs {
            let length = refs[chr];
            sorted_refs.insert(chr.clone(), length);
        }
        
        // populate row2parents and row2children indices
        let mut row2parents = HashMap::<usize, BTreeSet<usize>>::new();
        let mut row2children = HashMap::<usize, BTreeSet<usize>>::new();
        // make fake rows for GTF gene_id and transcript_id attributes
        let mut fake_rows = Vec::<Record>::new();
        let transcript_id_atom = Atom::from("transcript_id");
        let transcript_name_atom = Atom::from("transcript_name");
        let gene_id_atom = Atom::from("gene_id");
        let gene_name_atom = Atom::from("gene_name");
        
        match filetype {
            "gtf" => {
                let mut firstrow = HashMap::<&Record, usize>::new();
                for (row, record) in rows.iter().enumerate() {
                    // in GTF, if two rows both have identical fields (except attributes), they
                    // are the same feature.
                    let row = *firstrow.entry(record).or_insert(row);
                    if let Some(transcript_id) = record.attributes.get(&transcript_id_atom) {
                        let mut transcript_record = Record::new();
                        transcript_record.feature_type = Atom::from(transcript_type);
                        transcript_record.attributes
                            .insert(id_atom.clone(), transcript_id.clone());
                        if let Some(transcript_name) = record.attributes.get(&transcript_name_atom) {
                            transcript_record.attributes
                                .insert(name_atom.clone(), transcript_name.clone());
                        }
                        let transcriptrow = rows.len() + fake_rows.len();
                        id2row.insert(transcript_id.clone(), transcriptrow);
                        row2parents.entry(row).or_insert_with(BTreeSet::new).insert(transcriptrow);
                        row2children.entry(transcriptrow)
                            .or_insert_with(BTreeSet::new)
                            .insert(row);
                        fake_rows.push(transcript_record);

                        if let Some(gene_id) = record.attributes.get(&gene_id_atom) {
                            let mut gene_record = Record::new();
                            gene_record.feature_type = Atom::from(gene_type);
                            gene_record.attributes.insert(id_atom.clone(), gene_id.clone());
                            if let Some(gene_name) = record.attributes.get(&gene_name_atom) {
                                gene_record.attributes
                                    .insert(name_atom.clone(), gene_name.clone());
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
            "gff" => {
                let parent_atom = Atom::from("Parent");
                for (row, record) in rows.iter().enumerate() {
                    if let Some(parent) = record.attributes.get(&parent_atom) {
                        for p in parent.split(',') {
                            if let Some(parentrow) = id2row.get(&Atom::from(p)) {
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
                let mut seqname = HashMap::<Atom, u64>::new();
                let mut source = HashMap::<Atom, u64>::new();
                let mut start: Option<u64> = Option::None;
                let mut end: Option<u64> = Option::None;
                let mut strand = HashMap::<Atom, u64>::new();
                for c in &children {
                    let child = &rows[*c];
                    if start.is_none() || child.start < start.r()? {
                        start = Some(child.start);
                    }
                    if end.is_none() || child.end > end.r()? {
                        end = Some(child.end);
                    }
                    *seqname.entry(child.seqname.clone()).or_insert(0u64) += 1;
                    *source.entry(child.source.clone()).or_insert(0u64) += 1;
                    *strand.entry(child.strand.clone()).or_insert(0u64) += 1;
                }
                // enqueue the missing information
                let mut update = Record::new();
                if let Some((seqname, _)) = seqname.iter().max_by_key(|a| a.1) {
                    update.seqname = seqname.clone();
                }
                if let Some((source, _)) = source.iter().max_by_key(|a| a.1) {
                    update.source = source.clone();
                }
                if let Some(start) = start {
                    update.start = start;
                }
                if let Some(end) = end {
                    update.end = end;
                }
                if let Some((strand, _)) = strand.iter().max_by_key(|a| a.1) {
                    update.strand = strand.clone();
                }
                updates.push((row, update));
            }
        }
        // update the missing information
        for u in &mut updates {
            let update = &mut u.1;
            let mut record = &mut rows[u.0];
            record.seqname = update.seqname.clone();
            record.source = update.source.clone();
            record.start = update.start;
            record.end = update.end;
            record.strand = update.strand.clone();
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
        let mut tree = HashMap::<Atom, IntervalTree<u64, usize>>::new();
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
            chrmap: chrmap,
            vizchrmap: vizchrmap,
            refs: sorted_refs,
        })
    }
    
    fn to_gtf(&self, filename: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if filename == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(filename)?) });
            
        let id_atom = Atom::from("ID");
        let gene_id_atom = Atom::from("gene_id");
        let transcript_id_atom = Atom::from("transcript_id");
        for (row, record) in self.rows.iter().enumerate() {
            if let Some(transcript_rows) = self.row2parents.get(&row) {
                for transcript_row in transcript_rows {
                    let transcript = &self.rows[*transcript_row];
                    if let Some(transcript_id) = transcript.attributes.get(&id_atom) {
                        if let Some(gene_rows) = self.row2parents.get(transcript_row) {
                            for gene_row in gene_rows {
                                let gene = &self.rows[*gene_row];
                                if let Some(gene_id) = gene.attributes.get(&id_atom) {
                                    let mut rec = Record {
                                        row: row,
                                        seqname: record.seqname.clone(),
                                        source: record.source.clone(),
                                        feature_type: record.feature_type.clone(),
                                        start: record.start,
                                        end: record.end,
                                        score: record.score.clone(),
                                        strand: record.strand.clone(),
                                        frame: record.frame.clone(),
                                        attributes: LinkedHashMap::new(),
                                    };
                                    rec.attributes.insert(gene_id_atom.clone(), gene_id.clone());
                                    rec.attributes.insert(transcript_id_atom.clone(), transcript_id.clone());
                                    for (k,v) in &record.attributes {
                                        if k != "gene_id" && k != "transcript_id" {
                                            rec.attributes.insert(k.clone(), v.clone());
                                        }
                                    }
                                    output.write_fmt(format_args!("{}\n", rec.to_gtf()?))?;
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
    
    fn to_gff(&self, filename: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if filename == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(filename)?) });
            
        let id_atom = Atom::from("ID");
        for (row, record) in self.rows.iter().enumerate() {
            let mut rec = Record {
                row: row,
                seqname: record.seqname.clone(),
                source: record.source.clone(),
                feature_type: record.feature_type.clone(),
                start: record.start,
                end: record.end,
                score: record.score.clone(),
                strand: record.strand.clone(),
                frame: record.frame.clone(),
                attributes: LinkedHashMap::new(),
            };
            let mut parents = Vec::<Atom>::new();
            if let Some(parent_rows) = self.row2parents.get(&row) {
                for parent_row in parent_rows {
                    let parent = &self.rows[*parent_row];
                    if let Some(parent_id) = parent.attributes.get(&id_atom) {
                        parents.push(parent_id.clone());
                    }
                }
            }
            for (k,v) in &record.attributes {
                if k == "Parent" {
                    rec.attributes.insert(k.clone(), Atom::from(parents.iter().map(|s| s.as_ref()).join(",")));
                }
                else {
                    rec.attributes.insert(k.clone(), v.clone());
                }
            }
            output.write_fmt(format_args!("{}\n", rec.to_gff()?))?;
        }
        Ok(())
    }
    
    fn to_json(&self, json_file: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if json_file == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(json_file)?) });
        
        serde_json::to_writer_pretty(&mut output, &self)?;
        Ok(())
    }
    
    fn to_bed(&self, 
        bed_file: &str, 
        exon_types: &[String],
        cds_types: &[String],
        transcript_types: &[String])
        -> Result<()>
    {
        let exon_types = exon_types.iter().map(|t| Atom::from(t.as_ref())).collect::<HashSet<Atom>>();
        let mut cds_types = cds_types.iter().map(|t| Atom::from(t.as_ref())).collect::<HashSet<Atom>>();
        if cds_types.is_empty() { cds_types.insert(Atom::from("CDS")); }
        let transcript_types = transcript_types.iter().map(|t| Atom::from(t.as_ref())).collect::<HashSet<Atom>>();
        let mut seen_transcript = HashSet::<usize>::new();
        let mut transcript_names = HashSet::<Atom>::new();
        
        let id_atom = Atom::from("ID");
        let name_atom = Atom::from("Name");
        let unsorted_bed_file = format!("{}.unsorted.bed", bed_file);
        {   let mut bw = BufWriter::new(File::create(&unsorted_bed_file)?);
            let mut chrs = self.tree.keys().collect::<Vec<_>>();
            chrs.sort_by_key(|a| a.as_bytes());
            for chr in chrs {
                for node in self.tree[chr].find(0..std::u64::MAX) {
                    let transcript_row = node.data();
                    if !seen_transcript.contains(transcript_row) {
                        let transcript = &self.rows[*transcript_row];
                        if transcript_types.is_empty() || transcript_types.contains(&transcript.feature_type) {
                            seen_transcript.insert(*transcript_row);
                            
                            // get the exon and CDS features
                            let mut exons = Vec::<usize>::new();
                            let mut cds_features = Vec::<usize>::new();
                            let mut exon_starts = HashSet::<u64>::new();
                            let mut exon_ends = HashSet::<u64>::new();
                            if let Some(child_rows) = self.row2children.get(transcript_row) {
                                for child_row in child_rows {
                                    let child = &self.rows[*child_row];
                                    if child.seqname == transcript.seqname {
                                        if cds_types.contains(&child.feature_type) {
                                            cds_features.push(*child_row);
                                        }
                                        else if exon_types.is_empty() || exon_types.contains(&child.feature_type) {
                                            exon_starts.insert(self.rows[*child_row].start-1);
                                            exon_ends.insert(self.rows[*child_row].end);
                                            exons.push(*child_row);
                                        }
                                    }
                                }
                            }
                            if exons.is_empty() { exons.push(*transcript_row); }
                            exons.sort_by_key(|a| self.rows[*a].start);
                            cds_features.sort_by_key(|a| self.rows[*a].start);
                            
                            // get the cds ranges by unsplicing the CDS features
                            let mut cdss = Vec::<Range<u64>>::new();
                            for cds_feature_row in cds_features {
                                let cds_feature = &self.rows[cds_feature_row];
                                if !cdss.is_empty() && 
                                    exon_ends.contains(&cdss[cdss.len()-1].end) && 
                                    exon_starts.contains(&(cds_feature.start-1)) 
                                {
                                    // extend current CDS if this is a splice
                                    if let Some(last) = cdss.last_mut() {
                                        (*last).end = cds_feature.end;
                                    }
                                }
                                else {
                                    // otherwise add new CDS
                                    cdss.push(cds_feature.start-1..cds_feature.end);
                                }
                            }
                            cdss.sort_by_key(|a| a.start);
                            
                            // write the bed record
                            let score = 0;
                            let item_rgb = &[0,0,0].iter().map(|v| v.to_string()).join(",");
                            let start = self.rows[exons[0]].start-1;
                            let end = self.rows[exons[exons.len()-1]].end;
                            // make a separate bed record for each cds
                            if cdss.is_empty() { cdss.push(start..start) }
                            for cds in cdss {
                                // choose a unique transcript name
                                let transcript_name = transcript.attributes.get(&id_atom).or_else(||
                                    transcript.attributes.get(&name_atom));
                                let mut transcript_name = match transcript_name {
                                    Some(t) => t.clone(),
                                    None => Atom::from(format!("{}:{}..{}:{}", 
                                        transcript.seqname, transcript.start-1, transcript.end, transcript.strand)),
                                };
                                while transcript_names.contains(&transcript_name) {
                                    lazy_static! {
                                        static ref TRANSCRIPT_RENAME: Regex = Regex::new(r"(?:\.([0-9]+))?$").unwrap();
                                    }
                                    transcript_name = Atom::from(TRANSCRIPT_RENAME.replace(&transcript_name.as_ref(), |caps: &Captures| {
                                        format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                                }
                                transcript_names.insert(transcript_name.clone());
                                
                                let chr = self.vizchrmap.get(&transcript.seqname).unwrap_or(&transcript.seqname);
                                // write the bed line
                                let line = &[
                                    chr.as_ref(),
                                    &start.to_string(),
                                    &end.to_string(),
                                    &transcript_name,
                                    &score.to_string(),
                                    &transcript.strand,
                                    &cds.start.to_string(),
                                    &cds.end.to_string(),
                                    item_rgb,
                                    &exons.len().to_string(),
                                    &exons.iter().map(|e| (self.rows[*e].end-self.rows[*e].start+1).to_string()).join(","),
                                    &exons.iter().map(|e| ((self.rows[*e].start-1)-start).to_string()).join(","),
                                ].join("\t");
                                writeln!(bw, "{}", line)?;
                            }
                        }
                    }
                }
            }
        }
        
        // sort the bed file
        let mut command = Command::new("sort");
        command.args(&["-k1,1","-k2,2n","-o", &bed_file, &unsorted_bed_file]);
        command.env("LC_COLLATE", "C");
        let status = command.spawn()?.wait()?;
        if !status.success() {
            return Err(format!(
                "Exit code {:?} returned from command {:?}", status.code(), command).into());
        }
        // remove the unsorted bed file
        std::fs::remove_file(&unsorted_bed_file)?;
        
        Ok(())
    }
    
    fn to_bigbed(&self, 
        file: &str, 
        exon_types: &[String],
        cds_types: &[String],
        transcript_types: &[String],
        trackdb: &mut BufWriter<Box<Write>>)
        -> Result<()> 
    {
        // write the bed file
        let bed_file = format!("{}.bed", file);
        self.to_bed(&bed_file, exon_types, cds_types, transcript_types)?;
        
        // write the genome file
        let genome_filename = format!("{}.genome", file);
        {   let mut genome_fh = File::create(&genome_filename)?;
            for (chr, length) in &self.refs {
                let chr = self.vizchrmap.get(chr).unwrap_or(chr);
                writeln!(genome_fh, "{}\t{}", chr, length)?;
            }
        }
        
        // call bedToBigBed
        let mut command = Command::new("bedToBigBed");
        command.args(&["-type=bed12", "-tab", "extraIndex=name", &bed_file, &genome_filename, file]);
        let status = command.spawn()?.wait()?;
        if !status.success() {
            return Err(format!(
                "Exit code {:?} returned from command {:?}", status.code(), command).into());
        }
        // remove the bed file
        std::fs::remove_file(&bed_file)?;
        // remove the genome file
        std::fs::remove_file(&genome_filename)?;
        
        define_encode_set! {
            pub PATH_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'+', '?', '&'}
        }
        let url = utf8_percent_encode(file, PATH_ENCODE_SET);
        let track_name = Path::new(file).file_stem().r()?.to_str().r()?;
        // write to the trackDb.txt file
        trackdb.write_fmt(format_args!(indoc!(r##"
        
            track {}
            type bigBed 12
            bigDataUrl {}
            shortLabel {}
            longLabel {}
            itemRgb On
            visibility hide
        "##), track_name, url, track_name, track_name))?;
        
        Ok(())
    }
}

#[derive(Clone)]
struct Cassette {
    range: Range<u64>,
    cassette_row: Option<usize>,
}
impl serde::Serialize for Cassette
{
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
        where S: serde::Serializer
    {
        let mut ia = serializer.serialize_struct("Cassette", 2)?;
        let range = vec![self.range.start, self.range.end];
        ia.serialize_field("range", &range)?;
        ia.serialize_field("cassette_row", &self.cassette_row)?;
        ia.end()
    }
}
struct ConstituitivePair {
    exon1_row: usize,
    exon2_row: usize,
    // start..end, optional cassette_row
    cassettes: Vec<Cassette>,
    mapped_reads: IntervalTree<u64,Atom>,
}
impl serde::Serialize for ConstituitivePair
{
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
        where S: serde::Serializer
    {
        let mut ia = serializer.serialize_struct("ConstituitivePair", 4)?;
        ia.serialize_field("exon1_row", &self.exon1_row)?;
        ia.serialize_field("exon2_row", &self.exon2_row)?;
        ia.serialize_field("cassettes", &self.cassettes)?;
        let mut mapped_reads = Vec::<((u64, u64),Atom)>::new();
        for node in self.mapped_reads.find(0..std::u64::MAX) {
            mapped_reads.push(((node.interval().start,node.interval().end), node.data().clone()));
        }
        ia.serialize_field("mapped_reads", &mapped_reads)?;
        ia.end()
    }
}

fn find_constituitive_exons(annot: &IndexedAnnotation,
                            options: &Options)
                            -> Result<Vec<ConstituitivePair>> {
    let gene_types: HashSet<_> = options.gene_type.iter().map(|t| Atom::from(t.as_ref())).collect();
    let transcript_types: HashSet<_> = options.transcript_type.iter().map(|t| Atom::from(t.as_ref())).collect();
    let exon_types: HashSet<_> = options.exon_type.iter().map(|t| Atom::from(t.as_ref())).collect();
    
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
                                    let parents: HashSet<_> = parent_rows.iter().collect();
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
                                range: (annot.rows[c].start-1)..annot.rows[c].end, 
                                cassette_row: Some(c),
                            }).collect(),
                        mapped_reads: IntervalTree::new(),
                    }).
                    collect();
                exonpairs.append(&mut valid_constituitive_pairs);
            }
        }
    }
    // sort the pairs, by chr, start
    // this is required by bedToBigBed
    exonpairs.sort_by(|a,b| 
        annot.rows[a.exon1_row].seqname.as_bytes().cmp(annot.rows[b.exon1_row].seqname.as_bytes()).
        then_with(|| annot.rows[a.exon1_row].start.cmp(&annot.rows[b.exon1_row].start)));
    Ok(exonpairs)
}

fn bed2bigbed(
    bed_file: &str, 
    bigbed_file: &str, 
    refs: &LinkedHashMap<Atom,u64>, 
    sort_bed: bool,
    remove_bed: bool,
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<()> 
{
    // write the genome file
    let genome_filename = format!("{}.genome", bigbed_file);
    {   let mut genome_fh = File::create(&genome_filename)?;
        for (chr, length) in refs {
            writeln!(genome_fh, "{}\t{}", chr, length)?;
        }
    }
    
    // sort the bed file
    let sorted_bed_file = format!("{}.sorted.bed", bed_file);
    if sort_bed {
        let mut command = Command::new("sort");
        command.args(&["-k1,1","-k2,2n","-o", &sorted_bed_file, bed_file]);
        command.env("LC_COLLATE", "C");
        let status = command.spawn()?.wait()?;
        if !status.success() {
            return Err(format!(
                "Exit code {:?} returned from command {:?}", status.code(), command).into());
        }
    }
    
    // call bedToBigBed
    let mut command = Command::new("bedToBigBed");
    command.args(&["-type=bed12", "-tab", "extraIndex=name", 
        if sort_bed {&sorted_bed_file} else {bed_file}, 
        &genome_filename, bigbed_file]);
    let status = command.spawn()?.wait()?;
    if !status.success() {
        return Err(format!(
            "Exit code {:?} returned from command {:?}", status.code(), command).into());
    }
    // remove the bed file
    if remove_bed { std::fs::remove_file(&bed_file)?; }
    if remove_bed && sort_bed { std::fs::remove_file(&sorted_bed_file)?; }
    // remove the genome file
    std::fs::remove_file(&genome_filename)?;
    
    // write to the trackDb file
    define_encode_set! {
        pub PATH_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'+', '?', '&'}
    }
    let url = utf8_percent_encode(bigbed_file, PATH_ENCODE_SET);
    let track_name = Path::new(bigbed_file).file_stem().r()?.to_str().r()?;
    // write to the trackDb.txt file
    trackdb.write_fmt(format_args!(indoc!(r##"
    
        track {}
        type bigBed 12
        bigDataUrl {}
        shortLabel {}
        longLabel {}
        itemRgb On
        visibility hide
    "##), track_name, url, track_name, track_name))?;
        
    Ok(())
}

fn read_sizes_file(sizes_file: &str, chrmap: &HashMap<Atom,Atom>) -> Result<LinkedHashMap<Atom,u64>> {
    let mut refs = HashMap::<Atom,u64>::new();
    let f = File::open(&sizes_file)?;
    let file = BufReader::new(&f);
    for line in file.lines() {
        let line = line?;
        let cols: Vec<&str> = line.split('\t').collect();
        if let Some(chr) = cols.get(0) {
            let chr = Atom::from(*chr);
            let chr = chrmap.get(&chr).unwrap_or(&chr);
            if let Some(size) = cols.get(1) {
                let size = size.parse::<u64>()?;
                refs.insert(chr.clone(), size);
            }
        }
    }
    let mut chrs = refs.keys().collect::<Vec<_>>();
    chrs.sort_by_key(|a| a.as_bytes());
    let mut sorted_refs = LinkedHashMap::<Atom,u64>::new();
    for chr in chrs {
        let size = refs[chr];
        sorted_refs.insert(chr.clone(), size);
    }
    Ok(sorted_refs)
}


fn get_bam_refs(bamfile: &str, chrmap: &HashMap<Atom,Atom>) -> Result<LinkedHashMap<Atom,u64>> {
    let mut refs = LinkedHashMap::<Atom,u64>::new();
    let bam = IndexedReader::from_path(bamfile)?;
    let header = bam.header();
    let target_names = header.target_names();
    for target_name in target_names {
        let tid = header.tid(target_name).r()?;
        let target_len = header.target_len(tid).r()? as u64;
        let target_name = Atom::from(std::str::from_utf8(target_name)?);
        let chr = chrmap.get(&target_name).unwrap_or(&target_name);
        refs.insert(chr.clone(), target_len);
    }
    Ok(refs)
}

fn reannotate_regions(
    annot: &IndexedAnnotation,
    pairs: &[ConstituitivePair], 
    bamfiles: &[String], 
    bamstrand: &[Option<bool>], 
    options: &Options,
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<Vec<ConstituitivePair>>
{
    let mut reannotated = Vec::<ConstituitivePair>::new();
    // bigwig histograms
    let mut plus_bw_histo = HashMap::<Atom,CsVecOwned<u64>>::new();
    let mut minus_bw_histo = HashMap::<Atom,CsVecOwned<u64>>::new();
    let mut start_plus_bw_histo = HashMap::<Atom,CsVecOwned<u64>>::new();
    let mut start_minus_bw_histo = HashMap::<Atom,CsVecOwned<u64>>::new();
    let mut end_plus_bw_histo = HashMap::<Atom,CsVecOwned<u64>>::new();
    let mut end_minus_bw_histo = HashMap::<Atom,CsVecOwned<u64>>::new();
    // allocate the bigwig histograms
    for (chr, length) in &annot.refs {
        if options.debug_bigwig.is_some() {
            plus_bw_histo.insert(chr.clone(), CsVec::empty(*length as usize));
            minus_bw_histo.insert(chr.clone(), CsVec::empty(*length as usize));
            start_plus_bw_histo.insert(chr.clone(), CsVec::empty(*length as usize));
            start_minus_bw_histo.insert(chr.clone(), CsVec::empty(*length as usize));
            end_plus_bw_histo.insert(chr.clone(), CsVec::empty(*length as usize));
            end_minus_bw_histo.insert(chr.clone(), CsVec::empty(*length as usize));
        }
    }
    // store the refseq names and sizes
    // iterate through each constituitive pair
    for pair in pairs {
        let exon1 = &annot.rows[pair.exon1_row];
        let exon2 = &annot.rows[pair.exon2_row];
        let start = annot.rows[pair.exon1_row].end;
        let end = annot.rows[pair.exon2_row].start-1;
        let reflength = annot.refs[&exon1.seqname];
        let mut histo = CsVecOwned::<u64>::empty(reflength as usize);
        let mut start_histo = CsVecOwned::<u64>::empty(reflength as usize);
        let mut end_histo = CsVecOwned::<u64>::empty(reflength as usize);
        let mut exon_regions = CsVecOwned::<u64>::empty(reflength as usize);
        let mut incomplete_starts = CsVecOwned::<u64>::empty(reflength as usize);
        let mut incomplete_ends = CsVecOwned::<u64>::empty(reflength as usize);
        let mut mapped_reads = IntervalTree::<u64, Atom>::new();
        let mut cassettes = Vec::<Cassette>::new();
        for (b, bamfile) in bamfiles.iter().enumerate() {
            let read1strand = bamstrand[b];
            let mut bam = IndexedReader::from_path(bamfile)?;
            let chr = &annot.rows[pair.exon1_row].seqname;
            let strand = &annot.rows[pair.exon1_row].strand;
            let strand_is_plus = strand == "+";
            let mut bw_histogram = 
                if strand_is_plus { plus_bw_histo.get_mut(chr).r()? } 
                else { minus_bw_histo.get_mut(chr).r()? };
            let mut start_bw_histogram = 
                if strand_is_plus { start_plus_bw_histo.get_mut(chr).r()? } 
                else { start_minus_bw_histo.get_mut(chr).r()? };
            let mut end_bw_histogram = 
                if strand_is_plus { end_plus_bw_histo.get_mut(chr).r()? }
                else { end_minus_bw_histo.get_mut(chr).r()? };
            
            // build the tid map for this bam file
            let mut tidmap = HashMap::<Atom,u32>::new();
            {   let header = bam.header();
                for target_name in header.target_names() {
                    let tid = header.tid(target_name).r()?;
                    let target_name = Atom::from(std::str::from_utf8(target_name)?);
                    let chr = annot.chrmap.get(&target_name).unwrap_or(&target_name);
                    tidmap.insert(chr.clone(), tid);
                }
            }
            if let Some(tid) = tidmap.get(chr) {
                // group reads by pairs
                let mut read_pairs = HashMap::<Atom,Vec<Vec<Range<u64>>>>::new();
                bam.seek(*tid, start as u32, end as u32)?;
                for read in bam.records() {
                    let read = read?;
                    let read_name = Atom::from(str::from_utf8(read.qname())?);
                    // make sure the read's strand matches
                    if let Some(is_read1strand) = read1strand {
                        if !(((is_read1strand == read.is_first_in_template()) == !read.is_reverse()) == strand_is_plus) {
                            continue;
                        }
                    }
                    let exons = cigar2exons(&read.cigar(), read.pos() as u64)?;
                    read_pairs.entry(read_name).or_insert_with(Vec::new).push(exons);
                }
                for (read_name,read_pair) in read_pairs {
                    // at least one of the read pairs must match a constituitive splice junction
                    let mut matches_splice = false;
                    for exons in &read_pair {
                        matches_splice = exons.iter().enumerate().any(|(i,exon)| 
                            (i < exons.len()-1 && exon.end == exon1.end) || 
                            (i > 0 && exon.start == exon2.start-1));
                        if matches_splice { break }
                    }
                    if !matches_splice { continue }
                    
                    for exons in &read_pair {
                        // write the internal reads histogram and mapped_reads
                        for exon in exons {
                            mapped_reads.insert(Interval::new(exon.clone())?, read_name.clone());
                            for pos in exon.start..exon.end {
                                histo[pos as usize] += 1;
                                bw_histogram[pos as usize] += 1;
                            }
                        }
                        
                        // write the start/stop histograms
                        let matches_splice = exons.iter().enumerate().any(|(i,exon)| 
                            (i < exons.len()-1 && exon.end == exon1.end) || 
                            (i > 0 && exon.start == exon2.start-1));
                        if matches_splice {
                            for (i, exon) in exons.iter().enumerate() {
                                // write start/stop histograms
                                if i > 0 {
                                    start_histo[exon.start as usize] += 1;
                                }
                                if i < exons.len()-1 {
                                    end_histo[exon.end as usize] += 1;
                                }
                                if options.debug_bigwig.is_some() {
                                    if i > 0 {
                                        start_bw_histogram[exon.start as usize] += 1;
                                    }
                                    if i < exons.len()-1 {
                                        end_bw_histogram[exon.end as usize] += 1;
                                    }
                                }
                                
                                // write the exon_regions histogram
                                if i > 0 && i < exons.len()-1 && start < exon.start && exon.end < end {
                                    for pos in exon.start..exon.end {
                                        exon_regions[pos as usize] += 1;
                                    }
                                }
                                // write incomplete starts/ends histograms
                                if i == 0 && start < exon.end && exon.end < end {
                                    incomplete_ends[exon.end as usize] += 1;
                                }
                                if i == exons.len()-1 && start < exon.start && exon.start < end {
                                    incomplete_starts[exon.start as usize] += 1;
                                }
                            }
                        }
                    }
                }
            }
            else {
                writeln!(stderr(), "Could not map tid for chr {} in bam file {}, skipping!", chr, bamfile)?;
            }
        }
        // fill in the incompletes
        if options.fill_incomplete_exons {
            for (i, _) in incomplete_starts.iter() {
                let mut last_end = None;
                for j in i..end as usize {
                    if end_histo.get(j).is_some() {
                        last_end = Some(j);
                    }
                    if histo.get(j).is_none() { break }
                }
                if let Some(last_end) = last_end {
                    for pos in i..last_end {
                        exon_regions[pos as usize] += 1;
                    }
                }
            }
            for (i, _) in incomplete_ends.iter() {
                let mut last_start = None;
                for j in (start as usize..i-1).rev() {
                    if start_histo.get(j).is_some() {
                        last_start = Some(j);
                    }
                    if j == 0 || histo.get(j-1).is_none() { break }
                }
                if let Some(last_start) = last_start {
                    for pos in last_start..i {
                        exon_regions[pos as usize] += 1;
                    }
                }
            }
        }
        // iterate through exon regions
        lazy_static! {
            static ref MAX_SLIPPAGE: Regex = Regex::new(r"^([0-9]+)(%?)$").unwrap();
        }
        let caps = MAX_SLIPPAGE.captures(&options.max_slippage).r()?;
        let max_slippage = caps.get(1).r()?.as_str().parse::<u64>()?;
        let mut exon_start = 0u64;
        let mut prev_i = std::usize::MAX;
        for (i, _) in exon_regions.iter() {
            if i-1 != prev_i {
                if prev_i != std::usize::MAX {
                    let exon_end = prev_i as u64;
                    let max_slippage = if caps.get(2).is_some() 
                        { ((max_slippage as f64 / 100f64)*((exon_end-exon_start) as f64)) as u64 }
                        else { max_slippage };
                    // find histogram peaks of the start and end
                    let mut max = 0u64;
                    let mut reannot_exon_start = exon_start;
                    for j in exon_start..std::cmp::min(exon_end, exon_start+max_slippage) {
                        if let Some(val) = start_histo.get(j as usize) {
                            if max < *val {  
                                max = *val;
                                reannot_exon_start = j;
                            }
                        }
                    }
                    let mut max = 0u64;
                    let mut reannot_exon_end = exon_end;
                    for j in std::cmp::max(reannot_exon_start, exon_end-max_slippage)..exon_end {
                        if let Some(val) = end_histo.get(j as usize) {
                            if max < *val {
                                max = *val;
                                reannot_exon_end = j;
                            }
                        }
                    }
                    let mut max = 0u64;
                    let mut reannot_exon_end2 = exon_end;
                    for j in std::cmp::max(exon_start, exon_end-max_slippage)..exon_end {
                        if let Some(val) = end_histo.get(j as usize) {
                            if max < *val {
                                max = *val;
                                reannot_exon_end2 = j;
                            }
                        }
                    }
                    let mut max = 0u64;
                    let mut reannot_exon_start2 = exon_start;
                    for j in exon_start..std::cmp::min(reannot_exon_end2, exon_start+max_slippage) {
                        if let Some(val) = start_histo.get(j as usize) {
                            if max < *val {  
                                max = *val;
                                reannot_exon_start2 = j;
                            }
                        }
                    }
                    let (exon_start, exon_end) = 
                        if reannot_exon_end-reannot_exon_start < reannot_exon_end2-reannot_exon_start2
                        {(reannot_exon_start2, reannot_exon_end2)}
                        else {(reannot_exon_start, reannot_exon_end)};
                    // store the reannotated cassette exon
                    cassettes.push(Cassette {
                        range: exon_start..exon_end,
                        cassette_row: None,
                    });
                }
                exon_start = i as u64;
            }
            prev_i = i;
        }
        reannotated.push(ConstituitivePair {
            exon1_row: pair.exon1_row,
            exon2_row: pair.exon2_row,
            cassettes: cassettes,
            mapped_reads: mapped_reads,
        });
    }
    if options.debug_bigwig.is_some() { 
        let start_prefix = format!("{}.start", &options.debug_bigwig.clone().r()?);
        let end_prefix = format!("{}.end", &options.debug_bigwig.clone().r()?);
        let vizrefs = annot.refs.iter().
            map(|(k,v)| (annot.vizchrmap.get(k).unwrap_or(k).clone(), *v)).
            collect::<LinkedHashMap<Atom,u64>>();
        write_bigwig(&options.debug_bigwig.clone().r()?, &plus_bw_histo, &vizrefs, &annot.vizchrmap, "+", trackdb, "debug_bigwig_+", true)?;
        write_bigwig(&start_prefix, &start_plus_bw_histo, &vizrefs, &annot.vizchrmap, "+", trackdb, "debug_bigwig_+", false)?;
        write_bigwig(&end_prefix, &end_plus_bw_histo, &vizrefs, &annot.vizchrmap, "+", trackdb, "debug_bigwig_+", false)?;
        write_bigwig(&options.debug_bigwig.clone().r()?, &minus_bw_histo, &annot.refs, &annot.vizchrmap, "-", trackdb, "debug_bigwig_-", true)?;
        write_bigwig(&start_prefix, &start_minus_bw_histo, &vizrefs, &annot.vizchrmap, "-", trackdb, "debug_bigwig_-", false)?;
        write_bigwig(&end_prefix, &end_minus_bw_histo, &vizrefs, &annot.vizchrmap, "-", trackdb, "debug_bigwig_-", false)?;
    }
    Ok(reannotated)
}

fn write_bigwig(
    file: &str, 
    histogram: &HashMap<Atom,CsVecOwned<u64>>, 
    refs: &LinkedHashMap<Atom,u64>,
    vizchrmap: &HashMap<Atom,Atom>,
    strand: &str,
    trackdb: &mut BufWriter<Box<Write>>,
    parent: &str,
    write_parent: bool) 
    -> Result<()> 
{
    // write the bedgraph file
    let bedgraph_file = format!("{}.{}.bedgraph", file, strand);
    {   let mut bw = BufWriter::new(File::create(&bedgraph_file)?);
        let mut start: usize = 0;
        let mut end: usize = 0;
        // sort the chrs by ascii values. This is required by bedGraphToBigWig
        let mut chrs = histogram.keys().map(|c| (c, vizchrmap.get(c).unwrap_or(c))).collect::<Vec<_>>();
        chrs.sort_by_key(|a| a.1.as_bytes());
        for (chr, vizchr) in chrs {
            let histo = &histogram[chr];
            let ref_length = refs[chr.into()] as usize;
            while start < ref_length {
                while histo[end] == histo[start] && end < ref_length {
                    end += 1
                }
                if histo[start] > 0 {
                    if strand == "-" {
                        writeln!(bw, "{}\t{}\t{}\t{}\n",
                                 vizchr, start, end, -(histo[start] as i64))?;
                    }
                    else {
                        writeln!(bw, "{}\t{}\t{}\t{}\n",
                                 vizchr, start, end, histo[start])?;
                    }
                }
                start = end;
            }
        }
    }
    
    // write the genome file
    let genome_filename = format!("{}.{}.bw.genome", file, strand);
    {   let mut genome_fh = File::create(&genome_filename)?;
        for (chr, length) in refs {
            writeln!(genome_fh, "{}\t{}", chr, length)?;
        }
    }
    
    // call bedGraphToBigWig
    let bigwig_file = format!("{}.{}.bw", file, strand);
    let mut command = Command::new("bedGraphToBigWig");
    command.args(&[&bedgraph_file, &genome_filename, &bigwig_file]);
    let status = command.spawn()?.wait()?;
    if !status.success() {
        return Err(format!(
            "Exit code {:?} returned from command {:?}", status.code(), command).into());
    }
    // remove the bedgraph file
    std::fs::remove_file(&bedgraph_file)?;
    // remove the genome file
    std::fs::remove_file(&genome_filename)?;
    
    // write to the trackDb file
    define_encode_set! {
        pub PATH_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'+', '?', '&'}
    }
    let url = utf8_percent_encode(&bigwig_file, PATH_ENCODE_SET);
    let track_name = Path::new(&bigwig_file).file_stem().r()?.to_str().r()?;
    // write to the trackDb.txt file
    if write_parent {
        let viewlimits = if strand == "-" { "-50:0" } else { "0:50" };
        trackdb.write_fmt(format_args!(indoc!(r##"

            track {}
            shortLabel {}
            longLabel {}
            visibility hide
            type bigWig
            container multiWig
            aggregate none
            maxHeightPixels: 100:50:8
            viewLimits {}
            priority 1
            "##), parent, parent, parent, viewlimits))?;
    }
    trackdb.write_fmt(format_args!(r##"
    track {}
    shortLabel {}
    longLabel {}
    type bigWig
    parent {}
    bigDataUrl {}
    "##, track_name, track_name, track_name, parent, url))?;
    
    Ok(())
}

fn get_bam_total_reads(bamfiles: &[String]) -> Result<u64> {
    let mut total_reads = 0u64;
    for bamfile in bamfiles {
        let mut cmd = Command::new("samtools");
        cmd.args(&["idxstats",bamfile]);
        let output = cmd.output()?;
        if !output.status.success() {
            return Err(format!("Command {:?} returned exit status {}", 
                cmd, output.status.code().unwrap_or(-1)).into());
        }
        let stdout = String::from_utf8(output.stdout)?;
        for line in stdout.lines() {
            let cols: Vec<&str> = line.split('\t').collect();
            if let Some(reads_str) = cols.get(2) {
                if let Ok(reads) = reads_str.parse::<u64>() {
                    total_reads += reads;
                }
            }
        }
    }
    Ok(total_reads)
}

#[derive(Serialize)]
struct RpkmStats {
    gene_row: usize,
    intron_rpkm: f64,
    max_cassette_rpkm: f64,
    total_constituitive_rpkm: f64,
    total_cassette_rpkm: f64,
}

fn compute_rpkm( 
    annot: &IndexedAnnotation, 
    pairs: &[ConstituitivePair],
    total_reads: u64,
    debug_rpkm_region_bigbed: &Option<String>,
    debug_mapped_reads_bigbed: &Option<String>,
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<Vec<RpkmStats>> 
{
    let mut rpkmstats = Vec::<RpkmStats>::new();
    {
        let mut region_bb: BufWriter<Box<Write>> = BufWriter::new(
            match debug_rpkm_region_bigbed.as_ref().map(String::as_ref) {
                Some("-") => Box::new(stdout()),
                Some(f) => Box::new(File::create(format!("{}.bed",f))?),
                None => Box::new(sink()),
            });
        let mut reads_bb: BufWriter<Box<Write>> = BufWriter::new(
            match debug_mapped_reads_bigbed.as_ref().map(String::as_ref) {
                Some("-") => Box::new(stdout()),
                Some(f) => Box::new(File::create(format!("{}.bed",f))?),
                None => Box::new(sink()),
            });
        let mut seen_read = HashSet::<String>::new();

        for pair in pairs {
            let mut transcript_parents = HashSet::<usize>::new();
            if let Some(parents) = annot.row2parents.get(&pair.exon1_row) {
                for parent in parents {
                    transcript_parents.insert(*parent);
                }
            }
            if let Some(parents) = annot.row2parents.get(&pair.exon2_row) {
                for parent in parents {
                    transcript_parents.insert(*parent);
                }
            }
            let mut gene_parents = BTreeSet::<usize>::new();
            for transcript_parent in &transcript_parents {
                if let Some(parents) = annot.row2parents.get(transcript_parent) {
                    for parent in parents {
                        gene_parents.insert(*parent);
                    }
                }
            }
            let chr = &annot.rows[pair.exon1_row].seqname;
            let chr = annot.vizchrmap.get(&chr).unwrap_or(&chr);
            // rpkm = (10^9 * num_mapped_reads_to_target)/(total_mapped_reads * mappable_target_size_bp)
            for gene_row in &gene_parents {
                let mapped_reads = &pair.mapped_reads;
                let constituitive_bases = annot.rows[pair.exon1_row].end-annot.rows[pair.exon1_row].start+1;
                let mut constituitive_reads = HashSet::<Atom>::new();
                for read in mapped_reads.find(annot.rows[pair.exon1_row].start-1..annot.rows[pair.exon1_row].end) {
                    constituitive_reads.insert(read.data().clone());
                }
                for read in mapped_reads.find(annot.rows[pair.exon2_row].start-1..annot.rows[pair.exon2_row].end) {
                    constituitive_reads.insert(read.data().clone());
                }
                let total_constituitive_rpkm = if total_reads == 0 || constituitive_bases == 0 { 0f64 } 
                else {
                    (1e10f64 * constituitive_reads.len() as f64) 
                        / (total_reads as f64 * constituitive_bases as f64)
                };
                
                let mut cassette_reads = HashMap::<Atom,Vec<Range<u64>>>::new();
                let mut cassette_features = Vec::<(Range<u64>,f64)>::new();
                let mut intron_reads = HashMap::<Atom,Vec<Range<u64>>>::new();
                let mut intron_features = Vec::<Range<u64>>::new();
                if pair.cassettes.is_empty() {
                    let intron_range = annot.rows[pair.exon2_row].start-1..annot.rows[pair.exon1_row].end;
                    for read in mapped_reads.find(&intron_range) {
                        let mut intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                        intron_read.push(read.interval().start..read.interval().end);
                    }
                    intron_features.push(intron_range);
                }
                else {
                    let intron_range = pair.cassettes[0].range.start..annot.rows[pair.exon1_row].end;
                    for read in mapped_reads.find(&intron_range) {
                        let mut intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                        intron_read.push(read.interval().start..read.interval().end);
                    }
                    intron_features.push(intron_range);
                    
                    for (i, cassette) in pair.cassettes.iter().enumerate() {
                        let reads: Vec<_> = mapped_reads.find(&cassette.range).collect();
                        for read in &reads {
                            let mut cassette_read = cassette_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                            cassette_read.push(read.interval().start..read.interval().end);
                        }
                        let exon_rpkm = (1e10f64 * reads.len() as f64) / 
                                        (total_reads as f64 * (cassette.range.end-cassette.range.start) as f64);
                        cassette_features.push((cassette.range.clone(), exon_rpkm));
                        if i < pair.cassettes.len()-1 {
                            let intron_range = cassette.range.end..pair.cassettes[i+1].range.start;
                            for read in mapped_reads.find(&intron_range) {
                                let mut intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                                intron_read.push(read.interval().start..read.interval().end);
                            }
                            intron_features.push(intron_range);
                        }
                    }
                    
                    let intron_range = pair.cassettes[pair.cassettes.len()-1].range.start..annot.rows[pair.exon2_row].start;
                    for read in mapped_reads.find(&intron_range) {
                            let mut intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                            intron_read.push(read.interval().start..read.interval().end);
                    }
                    intron_features.push(intron_range);
                }
                let cassette_bases = cassette_features.iter().map(|f| f.0.end-f.0.start).sum::<u64>();
                let intron_bases = intron_features.iter().map(|f| f.end-f.start).sum::<u64>();
                let max_exon_rpkm = cassette_features.iter().map(|f| f.1).fold(std::f64::NAN, f64::max);
                let total_cassette_rpkm = (1e10f64 * cassette_reads.len() as f64) / (total_reads as f64 * cassette_bases as f64);
                let intron_rpkm = (1e10f64 * intron_reads.len() as f64) / (total_reads as f64 * intron_bases as f64);
                rpkmstats.push(RpkmStats {
                    gene_row: *gene_row,
                    intron_rpkm: intron_rpkm,
                    max_cassette_rpkm: max_exon_rpkm,
                    total_constituitive_rpkm: total_constituitive_rpkm,
                    total_cassette_rpkm: total_cassette_rpkm,
                });
                
                // write cassette features to the region bed file
                let pair_name = get_pair_name(pair, annot);
                cassette_features.sort_by_key(|f| f.0.start);
                let cassette_name = format!("{}.cassettes", pair_name);
                let cassette_color = [255,0,0];
                if !cassette_features.is_empty() {
                    let max_exon = 
                        if let Some(me) = cassette_features.iter().max_by(|a,b| a.1.partial_cmp(&b.1).unwrap_or(Less)) 
                        { me.clone() } 
                        else { (cassette_features[0].0.start..cassette_features[0].0.start, 0f64) };
                    
                    region_bb.write_fmt(format_args!("{}\n", [
                        chr.as_ref(),
                        &cassette_features[0].0.start.to_string(),
                        &cassette_features[cassette_features.len()-1].0.end.to_string(),
                        &cassette_name,
                        &0.to_string(),
                        &annot.rows[pair.exon1_row].strand,
                        &max_exon.0.start.to_string(),
                        &max_exon.0.end.to_string(),
                        &cassette_color.iter().map(|v| v.to_string()).join(","),
                        &cassette_features.len().to_string(),
                        &cassette_features.iter().map(|f| (f.0.end-f.0.start).to_string()).join(","),
                        &cassette_features.iter().map(|f| (f.0.start-cassette_features[0].0.start).to_string()).join(","),
                    ].join("\t")))?;
                }
                // write intron features to the region bed file
                intron_features.sort_by_key(|f| f.start);
                let intron_name = format!("{}.introns", pair_name);
                let intron_color = [0,0,255];
                region_bb.write_fmt(format_args!("{}\n", [
                    chr.as_ref(),
                    &intron_features[0].start.to_string(),
                    &intron_features[intron_features.len()-1].end.to_string(),
                    &intron_name,
                    &0.to_string(),
                    &annot.rows[pair.exon1_row].strand,
                    &intron_features[0].start.to_string(),
                    &intron_features[0].start.to_string(),
                    &intron_color.iter().map(|v| v.to_string()).join(","),
                    &intron_features.len().to_string(),
                    &intron_features.iter().map(|f| (f.end-f.start).to_string()).join(","),
                    &intron_features.iter().map(|f| (f.start-intron_features[0].start).to_string()).join(","),
                ].join("\t")))?;
                
                // write to the mapped reads bed file
                for (read_name, ranges) in &mut cassette_reads {
                    let read_name = format!("{}.c", read_name);
                    if !seen_read.contains(&read_name) {
                        seen_read.insert(read_name.clone());
                        ranges.sort_by_key(|e| e.start);
                        reads_bb.write_fmt(format_args!("{}\n", [
                            chr.as_ref(),
                            &ranges[0].start.to_string(),
                            &ranges[ranges.len()-1].end.to_string(),
                            &read_name,
                            &0.to_string(),
                            &annot.rows[pair.exon1_row].strand,
                            &ranges[0].start.to_string(),
                            &ranges[0].start.to_string(),
                            &cassette_color.iter().map(|v| v.to_string()).join(","),
                            &ranges.len().to_string(),
                            &ranges.iter().map(|r| (r.end-r.start).to_string()).join(","),
                            &ranges.iter().map(|r| (r.start-ranges[0].start).to_string()).join(","),
                        ].join(",")))?;
                            
                    }
                }
                for (read_name, ranges) in &mut intron_reads {
                    let read_name = format!("{}.i", read_name);
                    if !seen_read.contains(&read_name) {
                        seen_read.insert(read_name.clone());
                        ranges.sort_by_key(|e| e.start);
                        reads_bb.write_fmt(format_args!("{}\n", [
                            chr.as_ref(),
                            &ranges[0].start.to_string(),
                            &ranges[ranges.len()-1].end.to_string(),
                            &read_name,
                            &0.to_string(),
                            &annot.rows[pair.exon1_row].strand,
                            &ranges[0].start.to_string(),
                            &ranges[0].start.to_string(),
                            &intron_color.iter().map(|v| v.to_string()).join(","),
                            &ranges.len().to_string(),
                            &ranges.iter().map(|r| (r.end-r.start).to_string()).join(","),
                            &ranges.iter().map(|r| (r.start-ranges[0].start).to_string()).join(","),
                        ].join(",")))?;
                    }
                }
            }
        }
    }
    // convert bed files to bigbed
    let vizrefs = annot.refs.iter().
        map(|(k,v)| (annot.vizchrmap.get(k).unwrap_or(k).clone(), *v)).
        collect::<LinkedHashMap<Atom,u64>>();
    match debug_rpkm_region_bigbed.as_ref().map(String::as_ref) {
        None | Some("-") => (),
        Some(f) => bed2bigbed(&format!("{}.bed", f), f, &vizrefs, true, true, trackdb)?
    };
    match debug_mapped_reads_bigbed.as_ref().map(String::as_ref) {
        None | Some("-") => (),
        Some(f) => bed2bigbed(&format!("{}.bed", f), f, &vizrefs, true, true, trackdb)?
    };
    Ok(rpkmstats)
}

fn write_rpkm_stats(outfile: &str, annot: &IndexedAnnotation, rpkmstats: &mut[RpkmStats]) -> Result<()> {
    let mut output: BufWriter<Box<Write>> = BufWriter::new(
        if outfile == "-" { Box::new(stdout()) } 
        else { Box::new(File::create(outfile)?) });
        
    // sort by intron_rpkm / max_exon_rpkm, descending
    rpkmstats.sort_by(|a, b|
        (b.intron_rpkm / b.max_cassette_rpkm).partial_cmp(&(a.intron_rpkm / a.max_cassette_rpkm)).unwrap_or(Less));
        
    // write the header
    output.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\n", 
        "gene_name",
        "intron_rpkm",
        "max_cassette_rpkm",
        "total_constituitive_rpkm",
        "total_cassette_rpkm",
        "intron_rpkm/max_cassette_rpkm",
    ))?;
    for rpkm in rpkmstats {
        let gene = &annot.rows[rpkm.gene_row];
        let gene_name = &gene.seqname;
        let ratio = rpkm.intron_rpkm / rpkm.max_cassette_rpkm;
        
        output.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\n", 
            gene_name, 
            rpkm.intron_rpkm, 
            rpkm.max_cassette_rpkm, 
            rpkm.total_constituitive_rpkm, 
            rpkm.total_cassette_rpkm, 
            ratio))?;
    }
    Ok(())
}

fn get_pair_name(pair: &ConstituitivePair, annot: &IndexedAnnotation) -> String {
    let mut gene_id = None;
    let id_atom = Atom::from("ID");
    let name_atom = Atom::from("Name");
    let gene_id_atom = Atom::from("gene_id");
    let transcript_id_atom = Atom::from("transcript_id");
    'FIND_GENE_ID:
    for transcript_row in &annot.row2parents[&pair.exon1_row] {
        for gene_row in &annot.row2parents[transcript_row] {
            if annot.rows[*gene_row].feature_type == *"gene" {
                if let Some(gene_name_attr) = annot.rows[*gene_row].attributes.get(&name_atom) {
                    gene_id = Some(gene_name_attr);
                    break 'FIND_GENE_ID;
                } else if let Some(gene_id_attr) = annot.rows[*gene_row].attributes.get(&id_atom) {
                    gene_id = Some(gene_id_attr);
                    break 'FIND_GENE_ID;
                }
            }
        }
        if let Some(transcript_name_attr) = annot.rows[*transcript_row].attributes.get(&name_atom) {
            gene_id = Some(transcript_name_attr);
            break 'FIND_GENE_ID;
        } else if let Some(transcript_id_attr) = annot.rows[*transcript_row].attributes.get(&id_atom) {
            gene_id = Some(transcript_id_attr);
            break 'FIND_GENE_ID;
        }
    }
    gene_id = gene_id.or_else(||
        annot.rows[pair.exon1_row].attributes.get(&gene_id_atom).or_else(||
        annot.rows[pair.exon1_row].attributes.get(&transcript_id_atom).or_else(||
        annot.rows[pair.exon1_row].attributes.get(&name_atom).or_else(|| 
        annot.rows[pair.exon1_row].attributes.get(&id_atom)))));
    // format the feature name    
    let name = format!("{}{}:{}..{}:{}", 
        match gene_id {
            Some(g) => format!("{}:", g),
            None => "".to_string(),
        },
        annot.rows[pair.exon1_row].seqname,
        annot.rows[pair.exon1_row].start-1,
        annot.rows[pair.exon2_row].end,
        annot.rows[pair.exon2_row].strand);
    name
}

fn write_exon_bigbed(
    pairs: &[ConstituitivePair], 
    annot: &IndexedAnnotation,
    file: &str, 
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<()> 
{
    // write the bed file
    let bed_file = format!("{}.bed", file);
    
    {   let mut bw = BufWriter::new(File::create(&bed_file)?);
        for pair in pairs {
            // find the gene ID
            let name = get_pair_name(pair, annot);
            // write the bed record
            let score = 0;
            let item_rgb = &[0,0,0].iter().map(|v| v.to_string()).join(",");
            let mut block_sizes = Vec::<u64>::new();
            let mut block_starts = Vec::<u64>::new();
            block_sizes.push(annot.rows[pair.exon1_row].end-annot.rows[pair.exon1_row].start+1);
            for cassette in &pair.cassettes {
                block_sizes.push(cassette.range.end-cassette.range.start);
                block_starts.push(cassette.range.start-annot.rows[pair.exon1_row].start+1);
            }
            block_sizes.push(annot.rows[pair.exon2_row].end-annot.rows[pair.exon2_row].start+1);
            let chr = &annot.rows[pair.exon1_row].seqname;
            let chr = annot.vizchrmap.get(&chr).unwrap_or(&chr);
            let line = &[
                chr.as_ref(),
                &(annot.rows[pair.exon1_row].start-1).to_string(),
                &annot.rows[pair.exon2_row].end.to_string(),
                &name,
                &score.to_string(),
                &annot.rows[pair.exon2_row].strand,
                &if pair.cassettes.is_empty() { (annot.rows[pair.exon1_row].start-1).to_string() } 
                else { pair.cassettes[0].range.start.to_string() },
                &if pair.cassettes.is_empty() { (annot.rows[pair.exon1_row].start-1).to_string() } 
                else { pair.cassettes[pair.cassettes.len()-1].range.end.to_string() },
                item_rgb,
                &(pair.cassettes.len()+2).to_string(),
                &block_sizes.iter().map(|v| v.to_string()).join(","),
                &block_starts.iter().map(|v| v.to_string()).join(","),
            ].join("\t");
            writeln!(bw, "{}", line)?;
        }
    }
    let vizrefs = annot.refs.iter().
        map(|(k,v)| (annot.vizchrmap.get(k).unwrap_or(k).clone(), *v)).
        collect::<LinkedHashMap<Atom,u64>>();
    bed2bigbed(&bed_file,file,&vizrefs,true,true,trackdb)?;
    Ok(())
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    // set defaults for feature types
    options.gene_type =
        (if options.gene_type.is_empty() { vec!["gene".to_string()] } 
         else { options.gene_type.clone() }).into_iter().collect();
    //options.transcript_type = 
    //    (if options.transcript_type.is_empty() { vec!["mRNA".to_string()] } 
    //     else { options.transcript_type.clone() }).into_iter().collect();
    options.exon_type =
        (if options.exon_type.is_empty() { vec!["exon".to_string()] }
         else { options.exon_type.clone() }).into_iter().collect();
    options.cds_type =
        (if options.cds_type.is_empty() { vec!["CDS".to_string()] }
         else { options.cds_type.clone() }).into_iter().collect();
    //  validate arguments
    lazy_static! {
        static ref VALIDATE_MAX_SLIPPAGE: Regex = Regex::new(r"^[0-9]+%?$").unwrap();
    }
    if !VALIDATE_MAX_SLIPPAGE.is_match(&options.max_slippage) {
        return Err(format!("Invalid max_slippage value {}: must be number of bases or % of exon length", options.max_slippage).into());
    }
    // set debug options if --debug flag is set
    if options.debug {
        if options.debug_annot_gff.is_none() { options.debug_annot_gff = Some(String::from("debug_annot.gff")) }
        if options.debug_annot_gtf.is_none() { options.debug_annot_gtf = Some(String::from("debug_annot.gtf")) }
        if options.debug_annot_bigbed.is_none() { options.debug_annot_bigbed = Some(String::from("debug_annot.bb")) }
        if options.debug_annot_json.is_none() { options.debug_annot_json = Some(String::from("debug_annot.json")) }
        if options.debug_exon_bigbed.is_none() { options.debug_exon_bigbed = Some(String::from("debug_exon.bb")) }
        if options.debug_bigwig.is_none() { options.debug_bigwig = Some(String::from("debug_intronrpkm")) }
        if options.debug_reannot_bigbed.is_none() { options.debug_reannot_bigbed = Some(String::from("debug_reannot.bb")) }
        if options.debug_rpkm_region_bigbed.is_none() { options.debug_rpkm_region_bigbed = Some(String::from("debug_rpkm.bb")) }
        if options.debug_mapped_reads_bigbed.is_none() { options.debug_mapped_reads_bigbed = Some(String::from("debug_mapped_reads.bb")) }
        if options.debug_rpkmstats_json.is_none() { options.debug_rpkmstats_json = Some(String::from("debug_rpkmstats.json")) }
        if options.debug_trackdb.is_none() { options.debug_trackdb = Some(String::from("trackDb.txt")) }
    }
    
    // set up the trackdb writer
    let mut trackdb: BufWriter<Box<Write>> = BufWriter::new(
        match options.debug_trackdb.as_ref().map(String::as_ref) {
            Some("-") => Box::new(stdout()),
            Some(f) => Box::new(File::create(f)?),
            None => Box::new(sink()),
        });

    // organize bamfiles by read strand type
    let mut bamfiles = Vec::<String>::new();
    bamfiles.append(&mut options.bam1.clone());
    bamfiles.append(&mut options.bam2.clone());
    bamfiles.append(&mut options.bam.clone());
    let mut bamstrand = Vec::<Option<bool>>::new();
    bamstrand.append(&mut options.bam1.iter().map(|_| Some(true)).collect());
    bamstrand.append(&mut options.bam2.iter().map(|_| Some(false)).collect());
    bamstrand.append(&mut options.bam.iter().map(|_| None).collect());
    // get the chromosome names and sizes from the first bam file
    if bamfiles.is_empty() {
        Options::clap().print_help()?;
        return Err("No bam files were passed in!".into());
    }
    let mrna_transcript_type = String::from("mRNA");
    let mut annot = if let Some(annotfile_gff) = options.annotfile_gff.clone() {
        writeln!(stderr(), "Reading annotation file {:?}", &annotfile_gff)?;
        IndexedAnnotation::from_gff(
            &annotfile_gff, 
            options.gene_type.get(0).r()?, 
            options.transcript_type.get(0).unwrap_or(&mrna_transcript_type), 
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else if let Some(annotfile_gtf) = options.annotfile_gtf.clone() {
        writeln!(stderr(), "Reading annotation file {:?}", &annotfile_gtf)?;
        IndexedAnnotation::from_gtf(&annotfile_gtf, 
            options.gene_type.get(0).r()?, 
            options.transcript_type.get(0).unwrap_or(&mrna_transcript_type), 
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else {
        Options::clap().print_help()?;
        return Err("No annotation file was given!".into());
    };
    writeln!(stderr(), "Getting refseq lengths from bam file {:?}", &bamfiles[0])?;
    let refs = match options.sizes_file.clone() {
        Some(sizes_file) => read_sizes_file(&sizes_file, &annot.chrmap)?,
        None => get_bam_refs(&bamfiles[0], &annot.chrmap)?,
    };
    annot.refs = refs;
    
    if let Some(debug_annot_gff) = options.debug_annot_gff.clone() {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_gff)?;
        annot.to_gff(&debug_annot_gff)?; 
    }
    if let Some(debug_annot_gtf) = options.debug_annot_gtf.clone() {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_gtf)?;
        annot.to_gtf(&debug_annot_gtf)?; 
    }
    if let Some(debug_annot_bigbed) = options.debug_annot_bigbed.clone() {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_bigbed)?;
        annot.to_bigbed(
            &debug_annot_bigbed, 
            &options.exon_type, 
            &options.cds_type, 
            &options.transcript_type,
            &mut trackdb)?; 
    }
    if let Some(debug_annot_json) = options.debug_annot_json.clone() {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_json)?;
        annot.to_json(&debug_annot_json)?; 
    }
    
    // get the total bam reads
    writeln!(stderr(), "Running samtools idxstats to get total bam read counts")?;
    let total_reads = get_bam_total_reads(&bamfiles)?;
    writeln!(stderr(), "Found {} total reads", total_reads)?;
    
    // find the constituitive exons
    writeln!(stderr(), "Searching for constituitive exons in the annotation")?;
    let exonpairs = find_constituitive_exons(&annot, &options)?;
    if let Some(debug_exon_bigbed) = options.debug_exon_bigbed.clone() {
        writeln!(stderr(), "Writing constituitive pairs to a bigbed file")?;
        write_exon_bigbed(&exonpairs, &annot, &debug_exon_bigbed, &mut trackdb)?;
    }
        
    writeln!(stderr(), "Reannotating cassette regions")?;
    let reannotated_pairs = reannotate_regions(
        &annot,
        &exonpairs, 
        &bamfiles, 
        &bamstrand,
        &options,
        &mut trackdb)?;
    if let Some(debug_reannot_bigbed) = options.debug_reannot_bigbed.clone() {
        writeln!(stderr(), "Writing reannotation to bigbed")?;
        write_exon_bigbed(&reannotated_pairs, &annot, &debug_reannot_bigbed, &mut trackdb)?;
    }
    
    writeln!(stderr(), "Computing RPKM stats")?;
    let mut rpkmstats = compute_rpkm(
        &annot, 
        &reannotated_pairs, 
        total_reads, 
        &options.debug_rpkm_region_bigbed,
        &options.debug_mapped_reads_bigbed,
        &mut trackdb)?;
    if let Some(debug_rpkmstats_json) = options.debug_rpkmstats_json.clone() {
        writeln!(stderr(), "Writing RPKM stats to {:?}", &debug_rpkmstats_json)?;
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if debug_rpkmstats_json == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(debug_rpkmstats_json)?) });
            serde_json::to_writer_pretty(&mut output, &rpkmstats)?;
    }
    writeln!(stderr(), "Writing RPKM stats to {:?}", &options.outfile)?;
    write_rpkm_stats(&options.outfile, &annot, &mut rpkmstats)?;
    Ok(())
}

fn main() {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "full");

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
