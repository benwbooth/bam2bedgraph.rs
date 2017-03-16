#![feature(ordering_chaining)]
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

// TODO: write trackdb, region bigbed, and log to stderr

#[derive(StructOpt, Debug)]
#[structopt(name = "intronrpkm", about = "Analyze RPKM values in intronic space")]
struct Options {
    // input files
    #[structopt(long="annot", short="a", help = "A genome annotation file in gtf/gff format", name="ANNOT_FILE")]
    annotfile: String,
    #[structopt(long="bam1", short="1", help = "The set of stranded .bam files to analyze where read1 indicates strand", name="BAMFILE1")]
    bam1: Vec<String>,
    #[structopt(long="bam2", short="2", help = "The set of stranded .bam files to analyze where read2 indicates strand", name="BAMFILE2")]
    bam2: Vec<String>,
    #[structopt(long="bam", short="u", help = "The set of unstranded .bam files to analyze", name="BAMFILE")]
    bam: Vec<String>,
    // output file
    #[structopt(long="out", short="o", help = "Output file", name="ANNOT_FILE", default_value="-")]
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
    #[structopt(long="use_annotated_cassettes", help = "Should the annotated cassettes be used in reannotating?")]
    use_annotated_cassettes: bool,
    #[structopt(long="splice_must_match", help = "Do we require splice junctions to match the constituitive splice starts?")]
    splice_must_match: bool,
    //#[structopt(long="compare_constituitive", help = "Should constituitve exons be compared?")]
    //compare_constituitive: bool,
    //#[structopt(long="compare_cassette", help = "Should cassette exons be compared?")]
    //compare_cassette: bool,
    
    // debug output files
    #[structopt(long="debug_annot", help = "Write the input annotation as a gtf/gff debug file", name="DEBUG_ANNOT_FILE")]
    debug_annot: Option<String>,
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
    #[structopt(long="debug_region_bigbed", help = "Output the reannotated intron/exon region features as a bigbed file", name="DEBUG_REGION_BIGBED_FILE")]
    debug_region_bigbed: Option<String>,
    #[structopt(long="debug_rpkmstats_json", help = "Dump the RpkmStats to a JSON file", name="DEBUG_RPKMSTATS_JSON")]
    debug_rpkmstats_json: Option<String>,
    #[structopt(long="debug_trackdb", help = "Write a UCSC trackDb.txt file with all the bigwigs/bigbeds", name="DEBUG_TRACKDB_FILE")]
    debug_trackdb: Option<String>,
}

#[derive(Default, Serialize, Deserialize)]
pub struct Record {
    row: usize,
    seqname: String,
    source: String,
    feature_type: String,
    start: u64,
    end: u64,
    score: String,
    strand: String,
    frame: String,
    attributes: LinkedHashMap<String, String>,
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
    fn from_row(row: usize, line: &str, filetype: &str) -> Result<Record> {
        lazy_static! {
            static ref COMMENT: Regex = Regex::new(r"^#").unwrap();
            static ref GTF_ATTR: Regex = Regex::new(r#"^(?P<key>\S+)\s+(?:"(?P<qval>[^"]*)"|(?P<val>\S+));\s*"#).unwrap();
        }
        if COMMENT.is_match(line) {
            return Err("Comment".into());
        }
        let fields: Vec<_> = line.split('\t').collect();
        Ok(Record {
            row: row,
            seqname: String::from(*fields.get(0).unwrap_or(&"")),
            source: String::from(*fields.get(1).unwrap_or(&"")),
            feature_type: String::from(*fields.get(2).unwrap_or(&"")),
            start: fields.get(3).unwrap_or(&"0").parse::<u64>()?,
            end: fields.get(4).unwrap_or(&"0").parse::<u64>()?,
            score: String::from(*fields.get(5).unwrap_or(&"")),
            strand: String::from(*fields.get(6).unwrap_or(&"")),
            frame: String::from(*fields.get(7).unwrap_or(&"")),
            attributes: if filetype == "gff" {
                fields.get(8).unwrap_or(&"").split(';').map(|a| {
                    let kv: Vec<_> = a.splitn(2, '=').collect();
                    let key = String::from(percent_decode(kv[0].as_bytes()).decode_utf8().unwrap());
                    let value = String::from(percent_decode(kv.get(1).unwrap_or(&"").as_bytes()).decode_utf8().unwrap());
                    (key, value)
                }).collect()
            } else if filetype == "gtf" {
                GTF_ATTR.captures_iter(fields.get(8).unwrap_or(&"")).map(|caps| {
                    (caps["key"].to_string(), 
                     caps.name("qval").unwrap_or_else(|| caps.name("val").unwrap()).as_str().to_string())
                }).collect()
            } else {
                return Err(format!("Don't know how to read filetype {}", filetype).into());
            },
        })
    }
    fn to_line(&self, filetype: &str) -> Result<String> {
        define_encode_set! {
            pub GFF_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'\t', '\r', '\n', ';', '%', '='}
        }
        match filetype {
            "gff" => Ok(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
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
                        utf8_percent_encode(v, GFF_ENCODE_SET))).
                    collect::<Vec<_>>().join(";"))),
            "gtf" => Ok(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                self.seqname, 
                self.source, 
                self.feature_type,
                self.start,
                self.end,
                self.score,
                self.strand,
                self.frame,
                self.attributes.iter().map(|(k,v)| 
                    format!("{} \"{}\";", k, v)).collect::<Vec<_>>().join(" "))),
            f => Err(format!("Unknown file format: {}", f).into()),
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
        let mut tree = HashMap::<String, Vec<((u64,u64),usize)>>::new();
        for (chr,interval_tree) in &self.tree {
            let mut nodes = Vec::<((u64, u64),usize)>::new();
            for node in interval_tree.find(0..std::u64::MAX) {
                nodes.push(((node.interval().start,node.interval().end), *node.data()));
            }
            tree.insert(chr.clone(), nodes);
        }
        ia.serialize_field("tree", &tree)?;
        ia.end()
    }
}
impl IndexedAnnotation {
    fn from_file(annotfile: &str, gene_type: &str, transcript_type: &str) -> Result<IndexedAnnotation> {
        // figure out the file type of the annotation file
        let caps = Regex::new(r"\.([^.]*)$")?.captures(annotfile).
            ok_or_else(|| format!("Could not determine extension of file {}", annotfile))?;
        let ext = caps.get(1).map_or("", |m| m.as_str()).to_lowercase();
        let filetype = match ext.as_ref() {
            "gff" | "gff3" => Ok("gff"),
            "gtf" | "gtf2" => Ok("gtf"),
            ext => {
                Err(format!("Don't know how to parse file {} with extension {}!",
                            annotfile,
                            ext))
            }
        }?;
        let mut id2row = HashMap::<String, usize>::new();
        let mut rows = Vec::<Record>::new();
        let f = File::open(&annotfile)?;
        let file = BufReader::new(&f);
        for (row, line) in file.lines().enumerate() {
            if let Ok(record) = Record::from_row(row, &line?, filetype) {
                if let Some(id) = record.attributes.get("ID") {
                    id2row.insert(id.clone(), row);
                }
                rows.push(record);
            }
        }
        // populate row2parents and row2children indices
        let mut row2parents = HashMap::<usize, BTreeSet<usize>>::new();
        let mut row2children = HashMap::<usize, BTreeSet<usize>>::new();
        // make fake rows for GTF gene_id and transcript_id attributes
        let mut fake_rows = Vec::<Record>::new();
        match filetype {
            "gtf" => {
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
            "gff" => {
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
    
    fn to_file(&self, filename: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if filename == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(filename)?) });
            
        let caps = Regex::new(r"\.([^.]*)$")?.captures(filename).
            ok_or_else(|| format!("Could not determine extension of file {}", filename))?;
        let ext = caps.get(1).map_or("", |m| m.as_str()).to_lowercase();
        let filetype = match ext.as_ref() {
            "gtf" | "gtf2" => "gtf",
            _ => "gff",
        };
        for (row, record) in self.rows.iter().enumerate() {
            match filetype {
                "gtf" => {
                    if let Some(transcript_rows) = self.row2parents.get(&row) {
                        for transcript_row in transcript_rows {
                            let transcript = &self.rows[*transcript_row];
                            if let Some(transcript_id) = transcript.attributes.get("ID") {
                                if let Some(gene_rows) = self.row2parents.get(transcript_row) {
                                    for gene_row in gene_rows {
                                        let gene = &self.rows[*gene_row];
                                        if let Some(gene_id) = gene.attributes.get("ID") {
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
                                            rec.attributes.insert("gene_id".to_string(), gene_id.clone());
                                            rec.attributes.insert("transcript_id".to_string(), transcript_id.clone());
                                            for (k,v) in &record.attributes {
                                                if k != "gene_id" && k != "transcript_id" {
                                                    rec.attributes.insert(k.clone(), v.clone());
                                                }
                                            }
                                            output.write_fmt(format_args!("{}\n", rec.to_line(&filetype)?))?;
                                        }
                                    }
                                }
                            }
                        }
                    }
                },
                "gff" => {
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
                    let mut parents = Vec::<String>::new();
                    if let Some(parent_rows) = self.row2parents.get(&row) {
                        for parent_row in parent_rows {
                            let parent = &self.rows[*parent_row];
                            if let Some(parent_id) = parent.attributes.get("ID") {
                                parents.push(parent_id.clone());
                            }
                        }
                    }
                    for (k,v) in &record.attributes {
                        if k == "Parent" {
                            rec.attributes.insert(k.clone(), parents.join(","));
                        }
                        else {
                            rec.attributes.insert(k.clone(), v.clone());
                        }
                    }
                    output.write_fmt(format_args!("{}\n", rec.to_line(&filetype)?))?;
                },
                _ => panic!(),
            }
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
    
    fn to_bigbed(&self, 
        file: &str, 
        refs: &LinkedHashMap<String,(u32,u32)>,
        exon_types: &[String],
        cds_types: &[String],
        transcript_types: &[String])
        -> Result<()> 
    {
        let exon_types = exon_types.iter().cloned().collect::<HashSet<String>>();
        let mut cds_types = cds_types.iter().cloned().collect::<HashSet<String>>();
        if cds_types.is_empty() { cds_types.insert("CDS".to_string()); }
        let transcript_types = transcript_types.iter().cloned().collect::<HashSet<String>>();
        let mut seen_transcript = HashSet::<usize>::new();
        let mut transcript_names = HashSet::<String>::new();
        
        let bed_file = format!("{}.bed", file);
        {   let mut bw = BufWriter::new(File::create(&bed_file)?);
            let mut chrs = self.tree.keys().collect::<Vec<_>>();
            chrs.sort_by(|a,b| a.as_bytes().cmp(b.as_bytes()));
            for chr in chrs {
                for node in self.tree[chr].find(0..std::u64::MAX) {
                    let transcript_row = node.data();
                    if !seen_transcript.contains(transcript_row) {
                        let transcript = &self.rows[*transcript_row];
                        if transcript_types.is_empty() || transcript_types.contains(&transcript.feature_type) {
                            seen_transcript.insert(*transcript_row);
                            
                            // choose a unique transcript name
                            let transcript_name = transcript.attributes.get("ID").or_else(||
                                transcript.attributes.get("Name"));
                            let mut transcript_name = match transcript_name {
                                Some(t) => t.to_string(),
                                None => format!("{}:{}..{}:{}", 
                                    transcript.seqname, transcript.start-1, transcript.end, transcript.strand),
                            };
                            while transcript_names.contains(&transcript_name) {
                                lazy_static! {
                                    static ref TRANSCRIPT_RENAME: Regex = Regex::new(r"(?:\.([0-9]+))?$").unwrap();
                                }
                                transcript_name = TRANSCRIPT_RENAME.replace(&transcript_name, |caps: &Captures| {
                                    format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}).to_string();
                            }
                            transcript_names.insert(transcript_name.clone());
                            
                            // get the exons and cds features
                            let mut exons = Vec::<usize>::new();
                            let mut cds = Vec::<usize>::new();
                            for child_row in &self.row2children[transcript_row] {
                                let child = &self.rows[*child_row];
                                if child.seqname == transcript.seqname {
                                    if exon_types.is_empty() || exon_types.contains(&child.feature_type) {
                                        exons.push(*child_row);
                                    }
                                    if cds_types.contains(&child.feature_type) {
                                        cds.push(*child_row);
                                    }
                                }
                            }
                            if exons.is_empty() { exons.push(*transcript_row); }
                            exons.sort_by(|a,b| self.rows[*a].start.cmp(&self.rows[*b].start));
                            cds.sort_by(|a,b| self.rows[*a].start.cmp(&self.rows[*b].start));
                            
                            // write the bed record
                            let score = 0;
                            let item_rgb = &[0,0,0].iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
                            let start = self.rows[exons[0]].start-1;
                            let end = self.rows[exons[exons.len()-1]].end;
                            let line = &[
                                transcript.seqname.clone(),
                                start.to_string(),
                                end.to_string(),
                                transcript_name,
                                score.to_string(),
                                transcript.strand.clone(),
                                (if cds.is_empty() { start } else { self.rows[cds[0]].start-1 }).to_string(),
                                (if cds.is_empty() { start } else { self.rows[cds[cds.len()-1]].end }).to_string(),
                                item_rgb.clone(),
                                exons.len().to_string(),
                                exons.iter().map(|e| (self.rows[*e].end-self.rows[*e].start+1).to_string()).collect::<Vec<_>>().join(","),
                                exons.iter().map(|e| ((self.rows[*e].start-1)-start).to_string()).collect::<Vec<_>>().join(","),
                            ].join("\t");
                            writeln!(bw, "{}", line)?;
                        }
                    }
                }
            }
        }
        
        // write the genome file
        let genome_filename = format!("{}.bb.genome", file);
        {   let mut genome_fh = File::create(&genome_filename)?;
            for (chr, &(_, length)) in refs {
                writeln!(genome_fh, "{}\t{}", chr, length)?;
            }
        }
        
        // call bedToBigBed
        let bigbed_file = format!("{}.bb", file);
        let mut command = Command::new("bedToBigBed");
        command.args(&["-type=bed12", "-tab", "extraIndex=name", &bed_file, &genome_filename, &bigbed_file]);
        let status = command.spawn()?.wait()?;
        if !status.success() {
            return Err(format!(
                "Exit code {:?} returned from command {:?}", status.code(), command).into());
        }
        // remove the bed file
        std::fs::remove_file(&bed_file)?;
        // remove the genome file
        std::fs::remove_file(&genome_filename)?;
        
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
    mapped_reads: IntervalTree<u64, String>,
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
        let mut mapped_reads = Vec::<((u64, u64),String)>::new();
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
    let gene_types: HashSet<_> = options.gene_type.iter().collect();
    let transcript_types: HashSet<_> = options.transcript_type.iter().collect();
    let exon_types: HashSet<_> = options.exon_type.iter().collect();
    
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

fn get_bam_refs(bamfile: &str) -> Result<LinkedHashMap<String,(u32,u32)>> {
    let mut refs = LinkedHashMap::<String,(u32,u32)>::new();
    let bam = IndexedReader::from_path(bamfile)?;
    let header = bam.header();
    let target_names = header.target_names();
    for target_name in target_names {
        let tid = header.tid(target_name).r()?;
        let target_len = header.target_len(tid).r()?;
        let target_name = std::str::from_utf8(target_name)?;
        refs.insert(target_name.to_string(), (tid, target_len));
    }
    Ok(refs)
}

fn reannotate_regions(
    annot: &IndexedAnnotation,
    pairs: &[ConstituitivePair], 
    bamfiles: &[String], 
    bamstrand: &[Option<bool>], 
    refs: &LinkedHashMap<String,(u32,u32)>,
    options: &Options) 
    -> Result<Vec<ConstituitivePair>>
{
    let mut reannotated = Vec::<ConstituitivePair>::new();
    // bigwig histograms
    let mut start_plus_bw_histo = HashMap::<String,Vec<u32>>::new();
    let mut start_minus_bw_histo = HashMap::<String,Vec<u32>>::new();
    let mut end_plus_bw_histo = HashMap::<String,Vec<u32>>::new();
    let mut end_minus_bw_histo = HashMap::<String,Vec<u32>>::new();
    // allocate the bigwig histograms
    for (chr, &(_, length)) in refs {
        if options.debug_bigwig.is_some() {
            start_plus_bw_histo.insert(chr.to_string(), vec![0u32; length as usize]);
            start_minus_bw_histo.insert(chr.to_string(), vec![0u32; length as usize]);
            end_plus_bw_histo.insert(chr.to_string(), vec![0u32; length as usize]);
            end_minus_bw_histo.insert(chr.to_string(), vec![0u32; length as usize]);
        }
    }
    // store the refseq names and sizes
    // iterate through each constituitive pair
    for pair in pairs {
        let exon1 = &annot.rows[pair.exon1_row];
        let exon2 = &annot.rows[pair.exon2_row];
        let mut start_histo = vec![0u32; (exon2.end-exon1.start+1) as usize];
        let mut end_histo = vec![0u32; (exon2.end-exon1.start+1) as usize];
        let mut mapped_reads = IntervalTree::<u64, String>::new();
        for (b, bamfile) in bamfiles.iter().enumerate() {
            let read1strand = bamstrand[b];
            let mut bam = IndexedReader::from_path(bamfile)?;
            let chr = &annot.rows[pair.exon1_row].seqname;
            if let Some(tid) = bam.header.tid(&chr.clone().into_bytes()) {
                let start = &annot.rows[pair.exon1_row].start-1;
                let end = &annot.rows[pair.exon2_row].end;
                let strand = &annot.rows[pair.exon1_row].strand;
                let strand_is_plus = strand == "+";
                let mut start_bw_histogram = 
                    if strand_is_plus { start_plus_bw_histo.get_mut(chr).r()? } 
                    else { start_minus_bw_histo.get_mut(chr).r()? };
                let mut end_bw_histogram = 
                    if strand_is_plus { end_plus_bw_histo.get_mut(chr).r()? }
                    else { end_minus_bw_histo.get_mut(chr).r()? };
                
                bam.seek(tid, start as u32, *end as u32)?;
                for read in bam.records() {
                    let read = read?;
                    // make sure the read's strand matches
                    if let Some(is_read1strand) = read1strand {
                        if !(((is_read1strand == read.is_first_in_template()) == !read.is_reverse()) == strand_is_plus) {
                            continue;
                        }
                    }
                    let mut exons = Vec::<(u64, u64)>::new();
                    cigar2exons(&mut exons, &read.cigar(), read.pos() as u64)?;
                    // calculate introns
                    let mut introns = Vec::<(u64, u64)>::new();
                    for (i, exon) in exons.iter().enumerate() {
                        mapped_reads.insert(Interval::new(exon.0..exon.1)?, String::from_utf8(read.qname().to_vec())?);
                        if i < exons.len()-1 {
                            introns.push((exon.1, exons[i+1].0));
                        }
                    }
                    let overlaps_constituitives = if options.splice_must_match {
                            introns.iter().any(|intron| intron.0 == exon1.start-1)
                        }
                        else {
                            exons.iter().any(|exon| 
                                (exon1.start-1 < exon.1 && exon.0 < exon1.end) || 
                                   (exon2.start-1 < exon.1 && exon.0 < exon2.end))
                        };
                    if overlaps_constituitives {
                        for intron in &introns {
                            start_histo[(intron.0-exon1.start+1) as usize] += 1;
                            end_histo[(intron.1-exon1.start+1) as usize] += 1;
                            
                            if options.debug_bigwig.is_some() {
                                start_bw_histogram[intron.0 as usize] += 1;
                                end_bw_histogram[intron.1 as usize] += 1;
                            }
                        }
                    }
                }
            }
        }
        if options.use_annotated_cassettes {
            // merge cassettes together in case any are overlapping
            let mut sorted_cassettes = pair.cassettes.clone();
            sorted_cassettes.sort_by(|a, b| a.range.start.cmp(&b.range.start));
            let mut exons = Vec::<Range<u64>>::new();
            for higher in &sorted_cassettes {
                if exons.is_empty() {
                    exons.push(higher.range.clone());
                }
                else {
                    let lower_start = exons[exons.len()-1].start;
                    let lower_end = exons[exons.len()-1].end;
                    if higher.range.start <= lower_end {
                        let upper_bound = std::cmp::max(lower_end, higher.range.end);
                        let exons_len = exons.len()-1;
                        exons[exons_len] = lower_start..upper_bound;
                    }
                    else {
                        exons.push(higher.range.clone());
                    }
                }
            }
            // find max peak for each exon start/stop
            let mut reannotated_cassettes = Vec::<Cassette>::new();
            for (i, exon) in exons.iter().enumerate() {
                // find the start/stop with the max histogram value
                let start_range_start = if i == 0 { 0usize } else { exons[i-1].end as usize };
                let start_range = start_range_start..exon.end as usize;
                let startidx = start_histo[start_range].iter().enumerate().max_by_key(|&(_, v)| v);
                let start = match startidx {
                    Some((s, _)) => s as u64 + exon1.start,
                    None => exon.start,
                };
                let end_range = if i < exons.len()-1 { exon.start as usize..exons[i+1].start as usize } 
                                else { exon.start as usize..end_histo.len() };
                let endidx = end_histo[end_range].iter().enumerate().max_by_key(|&(_, v)| v);
                let end = match endidx {
                    Some((e, _)) => e as u64 + exon1.start,
                    None => exon.end,
                };
                // add the reannotated cassette
                reannotated_cassettes.push(Cassette {range: start..end, cassette_row: None});
            }
            let reannotated_pair = ConstituitivePair {
                exon1_row: pair.exon1_row,
                exon2_row: pair.exon2_row,
                cassettes: pair.cassettes.clone(),
                mapped_reads: mapped_reads,
            };
            reannotated.push(reannotated_pair);
        }
        else {
          return Err("De-novo transcript reannotation is not yet supported!".into());
        }
    }
    if options.debug_bigwig.is_some() { 
        let start_prefix = format!("{}.start", &options.debug_bigwig.clone().r()?);
        let end_prefix = format!("{}.end", &options.debug_bigwig.clone().r()?);
        write_bigwig(&start_prefix, &start_plus_bw_histo, refs, "+")?;
        write_bigwig(&start_prefix, &start_minus_bw_histo, refs, "-")?;
        write_bigwig(&end_prefix, &end_plus_bw_histo, refs, "+")?;
        write_bigwig(&end_prefix, &end_minus_bw_histo, refs, "-")?;
    }
    Ok(reannotated)
}

fn write_bigwig(
    file: &str, 
    histogram: &HashMap<String,Vec<u32>>, 
    refs: &LinkedHashMap<String,(u32,u32)>,
    strand: &str) 
    -> Result<()> 
{
    // write the bedgraph file
    let bedgraph_file = format!("{}.{}.bedgraph", file, strand);
    {   let mut bw = BufWriter::new(File::create(&bedgraph_file)?);
        let mut start: usize = 0;
        let mut end: usize = 0;
        let mut chrs: Vec<_> = histogram.keys().collect();
        // sort the chrs by ascii values. This is required by bedGraphToBigWig
        chrs.sort_by_key(|a| a.as_bytes());
        for chr in chrs {
            let histo = &histogram[chr];
            let ref_length = refs[chr.into()].1 as usize;
            while start < ref_length {
                while (if end < histo.len() { histo[end] } else { 0 }) ==
                      (if start < histo.len() { histo[start] } else { 0 }) &&
                      end < ref_length {
                    end += 1
                }
                if (if start < histo.len() { histo[start] } else { 0 }) > 0 {
                    if strand == "-" {
                        writeln!(bw, "{}\t{}\t{}\t{}\n",
                                 chr, start, end, -(histo[start] as i64))?;
                    }
                    else {
                        writeln!(bw, "{}\t{}\t{}\t{}\n",
                                 chr, start, end, histo[start])?;
                    }
                }
                start = end;
            }
        }
    }
    
    // write the genome file
    let genome_filename = format!("{}.{}.bw.genome", file, strand);
    {   let mut genome_fh = File::create(&genome_filename)?;
        for (chr, &(_, length)) in refs {
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
    total_reads: u64) 
    -> Result<Vec<RpkmStats>> 
{
    let mut rpkmstats = Vec::<RpkmStats>::new();
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
        // rpkm = (10^9 * num_mapped_reads_to_target)/(total_mapped_reads * mappable_target_size_bp)
        for gene_row in &gene_parents {
            let mapped_reads = &pair.mapped_reads;
            let constituitive_bases = annot.rows[pair.exon1_row].end-annot.rows[pair.exon1_row].start+1;
            let mut constituitive_reads = HashSet::<String>::new();
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
            
            let mut cassette_bases = 0u64;
            let mut cassette_reads = HashSet::<String>::new();
            let mut intron_bases = 0u64;
            let mut intron_reads = HashSet::<String>::new();
            let mut max_exon_rpkm = 0f64;
            if pair.cassettes.is_empty() {
                intron_bases += (annot.rows[pair.exon2_row].start-1) - annot.rows[pair.exon1_row].end;
                for read in mapped_reads.find(annot.rows[pair.exon2_row].start-1..annot.rows[pair.exon1_row].end) {
                    intron_reads.insert(read.data().clone());
                }
            }
            else {
                intron_bases += pair.cassettes[0].range.start - annot.rows[pair.exon1_row].end;
                for read in mapped_reads.find(pair.cassettes[0].range.start..annot.rows[pair.exon1_row].end) {
                    intron_reads.insert(read.data().clone());
                }
                for (i, cassette) in pair.cassettes.iter().enumerate() {
                    cassette_bases += cassette.range.end-cassette.range.start;
                    let reads: Vec<_> = mapped_reads.find(&cassette.range).collect();
                    for read in &reads {
                        cassette_reads.insert(read.data().clone());
                    }
                    let exon_rpkm = (1e10f64 * reads.len() as f64) / (total_reads as f64 * (cassette.range.end-cassette.range.start) as f64);
                    if max_exon_rpkm < exon_rpkm {
                        max_exon_rpkm = exon_rpkm;
                    }
                    if i < pair.cassettes.len()-1 {
                        intron_bases += pair.cassettes[i+1].range.start - cassette.range.end;
                        for read in mapped_reads.find(cassette.range.end..pair.cassettes[i+1].range.start) {
                            intron_reads.insert(read.data().clone());
                        }
                    }
                }
                intron_bases += (annot.rows[pair.exon2_row].start-1) - pair.cassettes[pair.cassettes.len()-1].range.end;
                for read in mapped_reads.find(pair.cassettes[0].range.start..annot.rows[pair.exon1_row].end) {
                    intron_reads.insert(read.data().clone());
                }
            }
            let total_cassette_rpkm = (1e10f64 * cassette_reads.len() as f64) / (total_reads as f64 * cassette_bases as f64);
            let intron_rpkm = (1e10f64 * intron_reads.len() as f64) / (total_reads as f64 * intron_bases as f64);
            rpkmstats.push(RpkmStats {
                gene_row: *gene_row,
                intron_rpkm: intron_rpkm,
                max_cassette_rpkm: max_exon_rpkm,
                total_constituitive_rpkm: total_constituitive_rpkm,
                total_cassette_rpkm: total_cassette_rpkm,
            });
            break;
        }
    }
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

fn write_exon_bigbed(
    pairs: &[ConstituitivePair], 
    annot: &IndexedAnnotation,
    file: &str, 
    refs: &LinkedHashMap<String,(u32,u32)>) 
    -> Result<()> 
{
    // bigbed file is already sorted
    
    // write the bed file
    let bed_file = format!("{}.bed", file);
    
    {   let mut bw = BufWriter::new(File::create(&bed_file)?);
        for pair in pairs {
            // find the gene ID
            let mut gene_id = None;
            'FIND_GENE_ID:
            for transcript_row in &annot.row2parents[&pair.exon1_row] {
                for gene_row in &annot.row2parents[transcript_row] {
                    if annot.rows[*gene_row].feature_type == "gene" {
                        if let Some(gene_name_attr) = annot.rows[*gene_row].attributes.get("Name") {
                            gene_id = Some(gene_name_attr);
                            break 'FIND_GENE_ID;
                        } else if let Some(gene_id_attr) = annot.rows[*gene_row].attributes.get("ID") {
                            gene_id = Some(gene_id_attr);
                            break 'FIND_GENE_ID;
                        }
                    }
                }
                if let Some(transcript_name_attr) = annot.rows[*transcript_row].attributes.get("Name") {
                    gene_id = Some(transcript_name_attr);
                    break 'FIND_GENE_ID;
                } else if let Some(transcript_id_attr) = annot.rows[*transcript_row].attributes.get("ID") {
                    gene_id = Some(transcript_id_attr);
                    break 'FIND_GENE_ID;
                }
            }
            gene_id = gene_id.or_else(||
                annot.rows[pair.exon1_row].attributes.get("gene_id").or_else(||
                annot.rows[pair.exon1_row].attributes.get("transcript_id").or_else(||
                annot.rows[pair.exon1_row].attributes.get("Name").or_else(|| 
                annot.rows[pair.exon1_row].attributes.get("ID")))));
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
                
            // write the bed record
            let score = 0;
            let item_rgb = &[0,0,0].iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
            let mut block_sizes = Vec::<u64>::new();
            let mut block_starts = Vec::<u64>::new();
            block_sizes.push(annot.rows[pair.exon1_row].end-annot.rows[pair.exon1_row].start+1);
            for cassette in &pair.cassettes {
                block_sizes.push(cassette.range.end-cassette.range.start);
                block_starts.push(cassette.range.start-annot.rows[pair.exon1_row].start+1);
            }
            block_sizes.push(annot.rows[pair.exon2_row].end-annot.rows[pair.exon2_row].start+1);
            let line = &[
                annot.rows[pair.exon1_row].seqname.clone(),
                (annot.rows[pair.exon1_row].start-1).to_string(),
                annot.rows[pair.exon2_row].end.to_string(),
                name,
                score.to_string(),
                annot.rows[pair.exon2_row].strand.clone(),
                if pair.cassettes.is_empty() { (annot.rows[pair.exon1_row].start-1).to_string() } 
                else { pair.cassettes[0].range.start.to_string() },
                if pair.cassettes.is_empty() { (annot.rows[pair.exon1_row].start-1).to_string() } 
                else { pair.cassettes[pair.cassettes.len()-1].range.end.to_string() },
                item_rgb.clone(),
                (pair.cassettes.len()+2).to_string(),
                block_sizes.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(","),
                block_starts.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(","),
            ].join("\t");
            writeln!(bw, "{}", line)?;
        }
    }
    
    // write the genome file
    let genome_filename = format!("{}.bb.genome", file);
    {   let mut genome_fh = File::create(&genome_filename)?;
        for (chr, &(_, length)) in refs {
            writeln!(genome_fh, "{}\t{}", chr, length)?;
        }
    }
    
    // call bedToBigBed
    let bigbed_file = format!("{}.bb", file);
    let mut command = Command::new("bedToBigBed");
    command.args(&["-type=bed12", "-tab", "extraIndex=name", &bed_file, &genome_filename, &bigbed_file]);
    let status = command.spawn()?.wait()?;
    if !status.success() {
        return Err(format!(
            "Exit code {:?} returned from command {:?}", status.code(), command).into());
    }
    // remove the bed file
    std::fs::remove_file(&bed_file)?;
    // remove the genome file
    std::fs::remove_file(&genome_filename)?;
    
    Ok(())
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    // set defaults for feature types
    options.gene_type =
        (if options.gene_type.is_empty() { vec!["gene".to_string()] } 
         else { options.gene_type.clone() }).into_iter().collect();
    options.transcript_type = 
        (if options.transcript_type.is_empty() { vec!["mRNA".to_string()] } 
         else { options.transcript_type.clone() }).into_iter().collect();
    options.exon_type =
        (if options.exon_type.is_empty() { vec!["exon".to_string()] }
         else { options.exon_type.clone() }).into_iter().collect();
    options.cds_type =
        (if options.cds_type.is_empty() { vec!["CDS".to_string()] }
         else { options.cds_type.clone() }).into_iter().collect();
    if !options.use_annotated_cassettes {
        return Err("De-novo transcript reannotation is not yet supported!".into());
    }
    
    // set up the trackdb writer
    let trackdb: BufWriter<Box<Write>> = BufWriter::new(
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
        return Err("No bam files were passed in!".into());
    }
    let refs = get_bam_refs(&bamfiles[0])?;
    
    let annot = IndexedAnnotation::from_file(&options.annotfile, 
        options.gene_type.get(0).r()?, 
        options.transcript_type.get(0).r()?)?;
    if let Some(debug_annot) = options.debug_annot.clone() {
        annot.to_file(&debug_annot)?; 
    }
    if let Some(debug_annot_bigbed) = options.debug_annot_bigbed.clone() {
        annot.to_bigbed(&debug_annot_bigbed, &refs, &options.exon_type, &options.cds_type, &options.transcript_type)?; 
    }
    if let Some(debug_annot_json) = options.debug_annot_json.clone() {
        annot.to_json(&debug_annot_json)?; 
    }
    
    // get the total bam reads
    let total_reads = get_bam_total_reads(&bamfiles)?;
    
    // find the constituitive exons
    let exonpairs = find_constituitive_exons(&annot, &options)?;
    if let Some(debug_exon_bigbed) = options.debug_exon_bigbed.clone() {
        write_exon_bigbed(&exonpairs, &annot, &debug_exon_bigbed, &refs)?;
    }
        
    let reannotated_pairs = reannotate_regions(
        &annot,
        &exonpairs, 
        &bamfiles, 
        &bamstrand,
        &refs,
        &options)?;
    if let Some(debug_reannot_bigbed) = options.debug_reannot_bigbed.clone() {
        write_exon_bigbed(&reannotated_pairs, &annot, &debug_reannot_bigbed, &refs)?;
    }
    
    let mut rpkmstats = compute_rpkm(&annot, &reannotated_pairs, total_reads)?;
    if let Some(debug_rpkmstats_json) = options.debug_rpkmstats_json.clone() {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if debug_rpkmstats_json == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(debug_rpkmstats_json)?) });
            serde_json::to_writer_pretty(&mut output, &rpkmstats)?;
    }
    write_rpkm_stats(&options.outfile, &annot, &mut rpkmstats)?;
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