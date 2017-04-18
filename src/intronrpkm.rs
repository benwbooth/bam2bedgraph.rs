#![recursion_limit="128"]
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
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter, BufRead, Write};
use std::io::{stdout, stderr, sink};
use std::path::Path;
use std::sync::Arc;

extern crate bam2bedgraph;
use bam2bedgraph::errors::*;
use bam2bedgraph::cigar2exons;

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

extern crate linked_hash_map;
use linked_hash_map::LinkedHashMap;

extern crate serde;
use serde::ser::SerializeStruct;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

extern crate itertools;
use itertools::Itertools;

extern crate concurrent_hashmap;
use concurrent_hashmap::*;

extern crate futures;
use futures::Future;
extern crate futures_cpupool;
use futures_cpupool::CpuPool;

extern crate num_cpus;

extern crate ordered_float;
use ordered_float::OrderedFloat;

#[macro_use]
extern crate duct;

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
    #[structopt(long="outannot", help = "Output Annotation file", name="OUT_ANNOT_FILE")]
    outannot: Option<String>,
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
    #[structopt(long="max_iterations", help = "How many start/stop combinations before we skip this one?", name="MAX_ITERATIONS", default_value="200")]
    max_iterations: usize,
    #[structopt(long="fill_incomplete_exons", help = "Should we try to fill in exons which do not have two splice junctions?")]
    fill_incomplete_exons: bool,
    
    #[structopt(long="cpu_threads", short="t", help = "How many threads to use for processing", default_value="0")]
    cpu_threads: usize,
    
    // debug output files
    #[structopt(long="debug", help = "Output all debug files?")]
    debug: bool,
    #[structopt(long="debug_prefix", help = "Prefix to prepend to debug output files", default_value="debug")]
    debug_prefix: String,
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
    #[structopt(long="debug_rpkmstats_json", help = "Dump the RpkmStats to a JSON file", name="DEBUG_RPKMSTATS_JSON")]
    debug_rpkmstats_json: Option<String>,
    #[structopt(long="debug_trackdb", help = "Write a UCSC trackDb.txt file with all the bigwigs/bigbeds", name="DEBUG_TRACKDB_FILE")]
    debug_trackdb: Option<String>,
    #[structopt(long="debug_outannot_bigbed", help = "Write the output annotation as a bigBed file", name="DEBUG_OUTANNOT_BIGBED_FILE")]
    debug_outannot_bigbed: Option<String>,
}

#[derive(Default, Serialize, Deserialize,Clone)]
struct Record {
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
    fn from_row(row: usize, line: &str, filetype: &str, chrmap: &HashMap<String,String>) -> Result<Record> {
        lazy_static! {
            static ref COMMENT: Regex = Regex::new(r"^#").unwrap();
            static ref GTF_ATTR: Regex = Regex::new(r#"^(?P<key>\S+)\s+(?:"(?P<qval>[^"]*)"|(?P<val>\S+));\s*"#).unwrap();
        }
        if COMMENT.is_match(line) {
            return Err("Comment".into());
        }
        let fields: Vec<_> = line.split('\t').collect();
        let seqname = String::from(*fields.get(0).unwrap_or(&""));
        let seqname = chrmap.get(&seqname).unwrap_or(&seqname);
        Ok(Record {
            row: row,
            seqname: seqname.clone(),
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
                    (String::from(caps["key"].to_string()), 
                     String::from(caps.name("qval").unwrap_or_else(|| caps.name("val").unwrap()).as_str().to_string()))
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
    id2row: HashMap<String, usize>,
    row2parents: HashMap<usize, Vec<usize>>,
    row2children: HashMap<usize, Vec<usize>>,
    tree: HashMap<String, IntervalTree<u64, usize>>,
    chrmap: HashMap<String, String>,
    vizchrmap: HashMap<String, String>,
    refs: LinkedHashMap<String, u64>,
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
        let mut chrmap = HashMap::<String,String>::new();
        if let Some(charmap_file) = chrmap_file.clone() {
            let f = File::open(&charmap_file)?;
            let file = BufReader::new(&f);
            for line in file.lines() {
                let line = line?;
                let cols: Vec<&str> = line.split('\t').collect();
                if let Some(key) = cols.get(0) {
                    if let Some(value) = cols.get(1) {
                        chrmap.insert(String::from(*key), String::from(*value));
                    }
                }
            }
        }
        // read in the optional visualization chr map
        let mut vizchrmap = HashMap::<String,String>::new();
        if let Some(vizcharmap_file) = vizchrmap_file.clone() {
            let f = File::open(&vizcharmap_file)?;
            let file = BufReader::new(&f);
            for line in file.lines() {
                let line = line?;
                let cols: Vec<&str> = line.split('\t').collect();
                if let Some(key) = cols.get(0) {
                    if let Some(value) = cols.get(1) {
                        vizchrmap.insert(String::from(*key), String::from(*value));
                    }
                }
            }
        }
        
        let mut id2row = HashMap::<String, usize>::new();
        let mut rows = Vec::<Record>::new();
        let f = File::open(&annotfile)?;
        let file = BufReader::new(&f);
        let mut refs = HashMap::<String,u64>::new();
        for line in file.lines() {
            let row = rows.len();
            if let Ok(record) = Record::from_row(row, &line?, filetype, &chrmap) {
                if let Some(id) = record.attributes.get("ID") {
                    id2row.insert(id.clone(), row);
                }
                // get the max ref lengths
                let mut reflength = refs.entry(record.seqname.clone()).or_insert(record.end);
                if *reflength < record.end { *reflength = record.end }
                rows.push(record);
            }
        }
        // sort the refs
        let mut sorted_refs = LinkedHashMap::<String,u64>::new();
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
        
        match filetype {
            "gtf" => {
                let mut firstrow = HashMap::<&Record, usize>::new();
                for (row, record) in rows.iter().enumerate() {
                    // in GTF, if two rows both have identical fields (except attributes), they
                    // are the same feature.
                    let row = *firstrow.entry(record).or_insert(row);
                    if let Some(transcript_id) = record.attributes.get("transcript_id") {
                        let mut transcript_record = Record::new();
                        transcript_record.feature_type = String::from(transcript_type);
                        transcript_record.attributes
                            .insert(String::from("ID"), transcript_id.clone());
                        if let Some(transcript_name) = record.attributes.get("transcript_name") {
                            transcript_record.attributes
                                .insert(String::from("Name"), transcript_name.clone());
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
                            gene_record.feature_type = String::from(gene_type);
                            gene_record.attributes.insert(String::from("ID"), gene_id.clone());
                            if let Some(gene_name) = record.attributes.get("gene_name") {
                                gene_record.attributes
                                    .insert(String::from("Name"), gene_name.clone());
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
                            if let Some(parentrow) = id2row.get(&String::from(p)) {
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
                let mut seqname = HashMap::<String, u64>::new();
                let mut source = HashMap::<String, u64>::new();
                let mut start: Option<u64> = Option::None;
                let mut end: Option<u64> = Option::None;
                let mut strand = HashMap::<String, u64>::new();
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
                    .seqname.as_bytes()
                    .cmp(&rows[*b].seqname.as_bytes())
                    .then_with(|| rows[*a].start.cmp(&rows[*b].start))
            });
            *up = vec;
        }
        // sort row2children by coordinates
        let mut row2children_update = HashMap::<usize, Vec<usize>>::new();
        for (row, children) in row2children {
            let mut up = row2children_update.entry(row).or_insert_with(Vec::new);
            let mut vec: Vec<_> = children.into_iter().collect();
            vec.sort_by(|a, b| {
                rows[*a]
                    .seqname.as_bytes()
                    .cmp(&rows[*b].seqname.as_bytes())
                    .then_with(|| rows[*a].start.cmp(&rows[*b].start))
            });
            *up = vec;
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
            chrmap: chrmap,
            vizchrmap: vizchrmap,
            refs: sorted_refs,
        })
    }
    
    fn to_gtf(&self, filename: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if filename == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(filename)?) });
            
        for (row, record) in self.rows.iter().enumerate() {
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
                                    rec.attributes.insert(String::from("gene_id"), gene_id.clone());
                                    rec.attributes.insert(String::from("transcript_id"), transcript_id.clone());
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
                    rec.attributes.insert(k.clone(), String::from(parents.join(",")));
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
        transcript_types: &[String],
        gene_types: &[String])
        -> Result<()>
    {
        let exon_types = exon_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let mut cds_types = cds_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        if cds_types.is_empty() { cds_types.insert(String::from("CDS")); }
        let transcript_types = transcript_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let gene_types = gene_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let mut seen_transcript = HashSet::<usize>::new();
        let mut transcript_names = HashSet::<String>::new();
        
        let unsorted_bed_file = format!("{}.unsorted.bed", bed_file);
        {   let mut bw = BufWriter::new(File::create(&unsorted_bed_file)?);
            let mut chrs = self.tree.keys().collect::<Vec<_>>();
            chrs.sort_by_key(|a| a.as_bytes());
            for chr in chrs {
                for node in self.tree[chr].find(0..std::u64::MAX) {
                    let gene_row = node.data();
                    let gene = &self.rows[*gene_row];
                    if gene_types.is_empty() || gene_types.contains(&gene.feature_type) {
                        if let Some(ref gene_children) = self.row2children.get(gene_row) {
                            for transcript_row in gene_children.iter() {
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
                                                cdss.push((cds_feature.start-1)..cds_feature.end);
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
                                            let transcript_name = 
                                                transcript.attributes.get("transcript_name").or_else(||
                                                transcript.attributes.get("Name").or_else(||
                                                transcript.attributes.get("ID")));
                                            let mut transcript_name = match transcript_name {
                                                Some(t) => t.clone(),
                                                None => String::from(format!("{}:{}..{}:{}", 
                                                    transcript.seqname, transcript.start-1, transcript.end, transcript.strand)),
                                            };
                                            while transcript_names.contains(&transcript_name) {
                                                lazy_static! {
                                                    static ref TRANSCRIPT_RENAME: Regex = Regex::new(r"(?:\.([0-9]+))?$").unwrap();
                                                }
                                                transcript_name = String::from(TRANSCRIPT_RENAME.replace(&transcript_name.as_ref(), |caps: &Captures| {
                                                    format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                                            }
                                            transcript_names.insert(transcript_name.clone());
                                            
                                            let chr = self.vizchrmap.get(&transcript.seqname).unwrap_or(&transcript.seqname);
                                            // write the bed line
                                            let line = &[
                                                &chr,
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
                                            ].iter().join("\t");
                                            writeln!(bw, "{}", line)?;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // sort the bed file
        cmd!("sort","-k1,1","-k2,2n","-o",&bed_file,&unsorted_bed_file).env("LC_COLLATE","C").run()?;
        // remove the unsorted bed file
        std::fs::remove_file(&unsorted_bed_file)?;
        
        Ok(())
    }
    
    fn to_bigbed(&self, 
        file: &str, 
        exon_types: &[String],
        cds_types: &[String],
        transcript_types: &[String],
        gene_types: &[String],
        trackdb: &mut BufWriter<Box<Write>>)
        -> Result<()> 
    {
        // write the bed file
        let bed_file = format!("{}.bed", file);
        self.to_bed(&bed_file, exon_types, cds_types, transcript_types, gene_types)?;
        
        // write the genome file
        let genome_filename = format!("{}.genome", file);
        {   let mut genome_fh = File::create(&genome_filename)?;
            for (chr, length) in &self.refs {
                let chr = self.vizchrmap.get(chr).unwrap_or(chr);
                writeln!(genome_fh, "{}\t{}", chr, length)?;
            }
        }
        
        // call bedToBigBed
        cmd!("bedToBigBed","-type=bed12","-tab","-extraIndex=name", &bed_file, &genome_filename, &file).run()?;
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
        trackdb.flush()?;
        
        Ok(())
    }
}

#[derive(Clone, Serialize)]
struct Cassette {
    range: Range<u64>,
    cassette_row: Option<usize>,
}
impl std::fmt::Debug for Cassette {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{{range: {:?}, cassette_row: {:?}}}", 
            self.range, self.cassette_row)
    }    
}
#[derive(Clone, Serialize)]
struct ConstituitivePair {
    exon1_row: usize,
    exon2_row: usize,
    // start..end, optional cassette_row
    cassettes: Vec<Cassette>,
}
impl std::fmt::Debug for ConstituitivePair {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "exon1_row: {}, exon2_row: {}, cassettes: {:?}", 
            self.exon1_row, self.exon2_row, self.cassettes)
    }    
}

fn find_constituitive_exons(annot: &IndexedAnnotation,
                            options: &Options)
                            -> Result<Vec<ConstituitivePair>> {
    let gene_types: HashSet<_> = options.gene_type.iter().map(|t| String::from(t.as_ref())).collect();
    let transcript_types: HashSet<_> = options.transcript_type.iter().map(|t| String::from(t.as_ref())).collect();
    let exon_types: HashSet<_> = options.exon_type.iter().map(|t| String::from(t.as_ref())).collect();
    
    // set default feature types
    let mut exonpairs = Vec::<ConstituitivePair>::new();
    for (gene_row, gene) in annot.rows.iter().enumerate() {
        if gene_types.is_empty() || gene_types.contains(&gene.feature_type) {
            if let Some(feature_rows) = annot.row2children.get(&gene_row) {
                // get the transcript rows for this gene
                let mut transcript_rows = Vec::<usize>::new();
                'transcript_row:
                for transcript_row in feature_rows {
                    let transcript = &annot.rows[*transcript_row];
                    // make sure this is a transcript type
                    // make sure transcript coordinates are consistent with the gene coordinates
                    if (transcript_types.is_empty() || 
                            transcript_types.contains(&transcript.feature_type)) &&
                        transcript.seqname == gene.seqname &&
                        transcript.strand == gene.strand 
                    {
                        // make sure the transcript has at least 1 exon
                        if let Some(exon_rows) = annot.row2children.get(transcript_row) {
                            for exon_row in exon_rows {
                                let exon = &annot.rows[*exon_row];
                                if (exon_types.is_empty() || 
                                        exon_types.contains(&exon.feature_type)) &&
                                    exon.seqname == transcript.seqname && 
                                    exon.strand == transcript.strand 
                                {
                                    transcript_rows.push(*transcript_row);
                                    break 'transcript_row;
                                }
                            }
                        }
                    }
                }
                // (exon1_row, exon2_row) -> (cassette_rows, transcript_rows)
                let mut constituitive_pairs = HashMap::<(usize,usize),(HashSet<usize>,HashSet<usize>)>::new();
                // look in each transcript
                'transcript: 
                for transcript_row in &transcript_rows {
                    let transcript = &annot.rows[*transcript_row];
                    if let Some(feature_rows) = annot.row2children.get(transcript_row) {
                        // get the set of exon features we care about
                        let mut exon_rows = Vec::<usize>::new();
                        for feature_row in feature_rows {
                            let feature = &annot.rows[*feature_row];
                            if exon_types.contains(&feature.feature_type) {
                                // skip any transcripts with weird exons
                                if feature.seqname != transcript.seqname || feature.strand != transcript.strand {
                                    continue 'transcript;
                                }
                                exon_rows.push(*feature_row);
                            }
                        }
                        // get the constituitive exons
                        let mut constituitive = HashSet::<usize>::new();
                        'exon: 
                        for exon_row in &exon_rows {
                            if let Some(parent_rows) = annot.row2parents.get(exon_row) {
                                let parents: HashSet<_> = parent_rows.iter().collect();
                                // are there any transcripts in this gene that this exon
                                // does NOT have as a parent?
                                for transcript_row in &transcript_rows {
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
                // filter out constituitive pairs that do not cover every transcript
                let mut valid_constituitive_pairs: Vec<_> = constituitive_pairs.iter().
                    filter(|p| (p.1).1.len() == transcript_rows.len()).
                    map(|p| {
                        // sort the cassettes
                        let mut cassettes = 
                            (p.1).0.iter().map(|c| Cassette {
                                    range: (annot.rows[*c].start-1)..annot.rows[*c].end, 
                                    cassette_row: Some(*c),
                                }).collect::<Vec<_>>();
                        cassettes.sort_by_key(|a| a.range.start);
                        
                        ConstituitivePair {
                            exon1_row: (p.0).0,
                            exon2_row: (p.0).1,
                            cassettes: cassettes,
                        }}).
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
    refs: &LinkedHashMap<String,u64>, 
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
        cmd!("sort","-k1,1","-k2,2n","-o",&sorted_bed_file,&bed_file).env("LC_COLLATE","C").run()?;
    }
    
    // call bedToBigBed
    cmd!("bedToBigBed","-type=bed12","-tab","-extraIndex=name",
        if sort_bed {&sorted_bed_file} else {bed_file},
        &genome_filename,
        &bigbed_file).run()?;
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
    trackdb.flush()?;
        
    Ok(())
}

fn read_sizes_file(sizes_file: &str, chrmap: &HashMap<String,String>) -> Result<LinkedHashMap<String,u64>> {
    let mut refs = HashMap::<String,u64>::new();
    let f = File::open(&sizes_file)?;
    let file = BufReader::new(&f);
    for line in file.lines() {
        let line = line?;
        let cols: Vec<&str> = line.split('\t').collect();
        if let Some(chr) = cols.get(0) {
            let chr = String::from(*chr);
            let chr = chrmap.get(&chr).unwrap_or(&chr);
            if let Some(size) = cols.get(1) {
                let size = size.parse::<u64>()?;
                refs.insert(chr.clone(), size);
            }
        }
    }
    let mut chrs = refs.keys().collect::<Vec<_>>();
    chrs.sort_by_key(|a| a.as_bytes());
    let mut sorted_refs = LinkedHashMap::<String,u64>::new();
    for chr in chrs {
        let size = refs[chr];
        sorted_refs.insert(chr.clone(), size);
    }
    Ok(sorted_refs)
}


fn get_bam_refs(bamfile: &str, chrmap: &HashMap<String,String>) -> Result<LinkedHashMap<String,u64>> {
    let mut refs = LinkedHashMap::<String,u64>::new();
    let bam = IndexedReader::from_path(bamfile)?;
    let header = bam.header();
    let target_names = header.target_names();
    for target_name in target_names {
        let tid = header.tid(target_name).r()?;
        let target_len = header.target_len(tid).r()? as u64;
        let target_name = String::from(std::str::from_utf8(target_name)?);
        let chr = chrmap.get(&target_name).unwrap_or(&target_name);
        refs.insert(chr.clone(), target_len);
    }
    Ok(refs)
}

fn reannotate_pair(
    pair_name: &str,
    exon1: &Record,
    exon2: &Record,
    read_pairs: &HashMap<(String,String),Vec<Vec<Range<u64>>>>,
    debug_bigwig: &Option<String>,
    fill_incomplete_exons: bool,
    max_iterations: usize,
    bw_histogram: Arc<ConcHashMap<usize,i32>>,
    start_bw_histogram: Arc<ConcHashMap<usize,i32>>,
    end_bw_histogram: Arc<ConcHashMap<usize,i32>>) 
    -> Result<(ConstituitivePair, IntervalTree<u64,String>)> 
{
    let start = exon1.end as usize;
    let end = (exon2.start-1) as usize;
    let region_size = end-start;
    let mut start_histo = vec![0i32; region_size];
    let mut end_histo = vec![0i32; region_size];
    let mut histo = vec![0i32; region_size];
    let mut mapped_reads = IntervalTree::<u64, String>::new();
    let mut cassettes = Vec::<Cassette>::new();
    for (&(_,ref read_name),read_pair) in read_pairs {
        // at least one of the read pairs must match a constituitive splice junction
        let mut matches_splice = false;
        for exons in read_pair {
            matches_splice = exons.iter().enumerate().any(|(i,exon)| 
                (i < exons.len()-1 && exon.end == exon1.end) || 
                (i > 0 && exon.start == exon2.start-1));
            if matches_splice { break }
        }
        if !matches_splice { continue }
        
        for exons in read_pair {
            // write the internal reads histogram and mapped_reads
            for exon in exons {
                mapped_reads.insert(Interval::new(exon.clone())?, read_name.clone());
                for pos in std::cmp::max(start as u64, exon.start)..std::cmp::min(end as u64, exon.end) {
                    if debug_bigwig.is_some() {
                        bw_histogram.upsert(pos as usize, 1, &|v| *v += 1);
                    }
                }
            }
            
            // write the start/stop histograms
            for (i, exon) in exons.iter().enumerate() {
                // write start/stop histograms
                if i > 0 {
                    if start <= (exon.start as usize) && (exon.start as usize) < end {
                        start_histo[exon.start as usize - start] += 1;
                    }
                }
                if i < exons.len()-1 {
                    if start <= (exon.end as usize) && (exon.end as usize) < end {
                        end_histo[std::cmp::max(0, exon.end as i64 - start as i64 - 1) as usize] += 1;
                    }    
                }
                if debug_bigwig.is_some() {
                    if i > 0 {
                        if start <= (exon.start as usize) && (exon.start as usize) < end {
                            start_bw_histogram.upsert(exon.start as usize, 1, &|v| *v += 1);
                        }
                    }
                    if i < exons.len()-1 {
                        if start <= (exon.end as usize) && (exon.end as usize) < end {
                            end_bw_histogram.upsert(std::cmp::max(0, exon.end as i64-1) as usize, 1, &|v| *v += 1);
                        }
                    }
                }
                
                // write the exon_regions histogram
                if start < (exon.end as usize) && (exon.start as usize) < end {
                    for pos in std::cmp::max(start as u64, exon.start)..std::cmp::min(end as u64, exon.end) {
                        histo[pos as usize - start] += 1;
                    }
                }
            }
        }
    }
    // iterate through exon regions
    let mut exon_start = start;
    let mut exon_value = histo[exon_start-start];
    let mut exon_regions = Vec::<Range<usize>>::new();
    for (i, value) in histo.iter().enumerate() {
        if (*value == 0) != (exon_value == 0) {
            if exon_value > 0 {
                let exon_end = i+start;
                exon_regions.push(exon_start..exon_end);
            }
            exon_start = i+start;
            exon_value = *value;
        }
    }
    'EXON_REGION:
    for exon_region in exon_regions {
        let mut starts = Vec::<(usize,i32)>::new();
        for (i, value) in start_histo[(exon_region.start-start)..(exon_region.end-start)].iter().enumerate() {
            if *value > 0 {
                starts.push((i+exon_region.start,*value));
            }
        }
        let mut ends = Vec::<(usize,i32)>::new();
        for (i, value) in end_histo[(exon_region.start-start)..(exon_region.end-start)].iter().enumerate() {
            if *value > 0 {
                ends.push((i+exon_region.start,*value));
            }
        }
        let mut iterations = 0;
        let mut start_score = 0;
        let mut end_score = 0;
        let mut best_start = exon_region.start;
        let mut best_end = exon_region.end;
        for &(s, sscore) in &starts {
            for &(e, escore) in &ends {
                if s < e && 
                    (start_score+end_score < sscore+escore 
                        || (start_score+end_score == sscore+escore && best_end-best_start < e-s)) 
                {
                    start_score = sscore;
                    end_score = escore;
                    best_start = s;
                    best_end = e;
                }
                iterations += 1;
                if iterations > max_iterations {
                    writeln!(stderr(), "More than {} iterations on pair {}", max_iterations, pair_name)?;
                    continue 'EXON_REGION;
                }
            }
        }
        // store the reannotated cassette exon
        if fill_incomplete_exons || (start_score > 0  && end_score > 0) {
            cassettes.push(Cassette {
                range: best_start as u64..best_end as u64,
                cassette_row: None,
            });
        }
    }
    let reannotpair = ConstituitivePair {
        exon1_row: exon1.row,
        exon2_row: exon2.row,
        cassettes: cassettes,
    };
    //writeln!(stderr(), "Writing reannotated pair: {:?}", reannotpair)?;
    Ok((reannotpair, mapped_reads))
}

fn reannotate_regions(
    annot: &Arc<IndexedAnnotation>,
    pairs: &[ConstituitivePair], 
    bamfiles: &[String], 
    bamstrand: &[Option<bool>], 
    total_reads: u64,
    options: &Options,
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<(Vec<ConstituitivePair>,Vec<RpkmStats>)>
{
    let num_cpus = num_cpus::get();
    let pool_cpus = if options.cpu_threads==0 {num_cpus} else {options.cpu_threads};
    
    // bigwig histograms
    let mut plus_bw_histo = HashMap::<String,Arc<ConcHashMap<usize,i32>>>::new();
    let mut minus_bw_histo = HashMap::<String,Arc<ConcHashMap<usize,i32>>>::new();
    let mut start_plus_bw_histo = HashMap::<String,Arc<ConcHashMap<usize,i32>>>::new();
    let mut start_minus_bw_histo = HashMap::<String,Arc<ConcHashMap<usize,i32>>>::new();
    let mut end_plus_bw_histo = HashMap::<String,Arc<ConcHashMap<usize,i32>>>::new();
    let mut end_minus_bw_histo = HashMap::<String,Arc<ConcHashMap<usize,i32>>>::new();
    // allocate the bigwig histograms
    {   let opts = || {
            concurrent_hashmap::Options {concurrency: pool_cpus as u16, ..Default::default() }
        };
        for (chr, _) in &annot.refs {
            plus_bw_histo.insert(chr.clone(),Arc::new(ConcHashMap::<usize,i32>::with_options(opts())));
            minus_bw_histo.insert(chr.clone(),Arc::new(ConcHashMap::<usize,i32>::with_options(opts())));
            start_plus_bw_histo.insert(chr.clone(),Arc::new(ConcHashMap::<usize,i32>::with_options(opts())));
            start_minus_bw_histo.insert(chr.clone(),Arc::new(ConcHashMap::<usize,i32>::with_options(opts())));
            end_plus_bw_histo.insert(chr.clone(),Arc::new(ConcHashMap::<usize,i32>::with_options(opts())));
            end_minus_bw_histo.insert(chr.clone(),Arc::new(ConcHashMap::<usize,i32>::with_options(opts())));
        }
    }
    // wrap them in an Arc
    let plus_bw_histo = Arc::new(plus_bw_histo);
    let minus_bw_histo = Arc::new(minus_bw_histo);
    let start_plus_bw_histo = Arc::new(start_plus_bw_histo);
    let start_minus_bw_histo = Arc::new(start_minus_bw_histo);
    let end_plus_bw_histo = Arc::new(end_plus_bw_histo);
    let end_minus_bw_histo = Arc::new(end_minus_bw_histo);
    
    let mut tidmaps = HashMap::<String,HashMap<String,u32>>::new();
    for bamfile in bamfiles {
        let bam = IndexedReader::from_path(bamfile)?;
        // build the tid map for this bam file
        let mut tidmap = HashMap::<String,u32>::new();
        {   let header = bam.header();
            for target_name in header.target_names() {
                let tid = header.tid(target_name).r()?;
                let target_name = String::from(std::str::from_utf8(target_name)?);
                let chr = annot.chrmap.get(&target_name).unwrap_or(&target_name);
                tidmap.insert(chr.clone(), tid);
            }
        }
        tidmaps.insert(bamfile.clone(), tidmap);
    }
    let tidmaps = Arc::new(tidmaps);
    
    let mut reannotated = Vec::new();
    let mut rpkmstats = Vec::new();
    let num_cpus = num_cpus::get();
    let pool = Arc::new(CpuPool::new(if options.cpu_threads==0 {num_cpus} else {options.cpu_threads}));
    let mut pair_futures = Vec::new();
    for pair in pairs {
        let exon1 = &annot.rows[pair.exon1_row];
        let exon2 = &annot.rows[pair.exon2_row];
        let start = exon1.end;
        let end = (exon2.start-1) as usize;
        let strand = &exon1.strand;
        let strand_is_plus = strand == "+";
        
        let chr = exon1.seqname.clone();
        let bw_histogram = 
            if strand_is_plus { plus_bw_histo[&chr].clone() } 
            else { minus_bw_histo[&chr].clone() };
        let start_bw_histogram = 
            if strand_is_plus { start_plus_bw_histo[&chr].clone() } 
            else { start_minus_bw_histo[&chr].clone() };
        let end_bw_histogram = 
            if strand_is_plus { end_plus_bw_histo[&chr].clone() }
            else { end_minus_bw_histo[&chr].clone() };
            
        // process the constituitivepair in parallel
        let exon1 = exon1.clone();
        let exon2 = exon2.clone();
        let debug_bigwig = options.debug_bigwig.clone().map(String::from);
        let fill_incomplete_exons = options.fill_incomplete_exons;
        let max_iterations = options.max_iterations;
        let bamfiles = Arc::new(bamfiles.to_vec());
        let bamstrand = Arc::new(bamstrand.to_vec());
        let tidmaps = tidmaps.clone();
        let annot = annot.clone();
        let pair_name = get_pair_name(pair, &annot);
        let pair_future = pool.spawn_fn(move ||->Result<(ConstituitivePair,RpkmStats)> {
            //get all the bam reads in parallel
            let mut read_pairs = HashMap::<(String,String),Vec<Vec<Range<u64>>>>::new();
            for (i,bamfile) in bamfiles.iter().enumerate() {
                let read1strand = bamstrand[i];
                let tidmap = &tidmaps[bamfile];
                if let Some(tid) = tidmap.get(&chr) {
                    let mut bam = IndexedReader::from_path(bamfile)?;
                    bam.seek(*tid, start as u32, end as u32)?;
                    for read in bam.records() {
                        let read = read?;
                        // make sure the read's strand matches
                        if let Some(is_read1strand) = read1strand {
                            if !(((is_read1strand == read.is_first_in_template()) == !read.is_reverse()) == strand_is_plus) {
                                continue;
                            }
                        }
                        let read_name = String::from(str::from_utf8(read.qname())?);
                        //writeln!(stderr(), "Looking at read: {}", read_name)?;
                        let exons = cigar2exons(&read.cigar(), read.pos() as u64)?;
                        read_pairs.entry((bamfile.clone(),read_name)).or_insert_with(Vec::new).push(exons);
                    }
                }
            }
            let (pair,mapped_reads) = reannotate_pair(
                &pair_name,
                &exon1,
                &exon2,
                &read_pairs,
                &debug_bigwig,
                fill_incomplete_exons,
                max_iterations,
                bw_histogram,
                start_bw_histogram,
                end_bw_histogram)?;
            let rpkmstats = compute_rpkm( 
                    &annot,
                    &pair,
                    &mapped_reads,
                    total_reads)?;
            Ok((pair, rpkmstats))
        });
        pair_futures.push(pair_future);
    }
    for future in pair_futures {
        match future.wait() {
            Ok((pair, stats)) => {
                reannotated.push(pair);
                rpkmstats.push(stats);
            }
            Err(ref e) => {
                writeln!(stderr(), "Got Err in reannotation: {:?}", e)?;
            }
        };
    }
    
    if options.debug_bigwig.is_some() { 
        let plus_parent = format!("{}_+", options.debug_prefix);
        let minus_parent = format!("{}_-", options.debug_prefix);
        let start_prefix = format!("{}.start", &options.debug_bigwig.clone().r()?);
        let end_prefix = format!("{}.end", &options.debug_bigwig.clone().r()?);
        write_bigwig(&options.debug_bigwig.clone().r()?, &plus_bw_histo, &annot.refs, &annot.vizchrmap, "+", trackdb, &plus_parent, true)?;
        write_bigwig(&start_prefix, &start_plus_bw_histo, &annot.refs, &annot.vizchrmap, "+", trackdb, &plus_parent, false)?;
        write_bigwig(&end_prefix, &end_plus_bw_histo, &annot.refs, &annot.vizchrmap, "+", trackdb, &plus_parent, false)?;
        write_bigwig(&options.debug_bigwig.clone().r()?, &minus_bw_histo, &annot.refs, &annot.vizchrmap, "-", trackdb, &minus_parent, true)?;
        write_bigwig(&start_prefix, &start_minus_bw_histo, &annot.refs, &annot.vizchrmap, "-", trackdb, &minus_parent, false)?;
        write_bigwig(&end_prefix, &end_minus_bw_histo, &annot.refs, &annot.vizchrmap, "-", trackdb, &minus_parent, false)?;
    }
    Ok((reannotated,rpkmstats))
}

fn write_bigwig(
    file: &str, 
    histogram: &HashMap<String,Arc<ConcHashMap<usize,i32>>>, 
    refs: &LinkedHashMap<String,u64>,
    vizchrmap: &HashMap<String,String>,
    strand: &str,
    trackdb: &mut BufWriter<Box<Write>>,
    parent: &str,
    write_parent: bool) 
    -> Result<()> 
{
    // write the bedgraph file
    let bedgraph_file = format!("{}.{}.bedgraph", file, strand);
    {   let mut bw = BufWriter::new(File::create(&bedgraph_file)?);
        // sort the chrs by ascii values. This is required by bedGraphToBigWig
        let mut chrs = histogram.keys().map(|c| (c, vizchrmap.get(c).unwrap_or(c))).collect::<Vec<_>>();
        chrs.sort_by_key(|a| a.1.as_bytes());
        for (chr, vizchr) in chrs {
            let histo = &histogram[chr];
            let reflength = refs[chr] as usize;
            let mut start = 0usize;
            let mut start_value = 0i32;
            if let Some(start_val) = histo.find(&start) {
                start_value = *start_val.get();
            }
            for i in 0usize..reflength {
                let mut value = 0i32;
                if let Some(v) = histo.find(&i) {
                    value = *v.get();
                }
                if value != start_value {
                    if start_value > 0i32 {
                        if strand == "-" {
                            writeln!(bw, "{}\t{}\t{}\t{}\n",
                                     vizchr, start, i, -(start_value as i64))?;
                        }
                        else {
                            writeln!(bw, "{}\t{}\t{}\t{}\n",
                                     vizchr, start, i, start_value)?;
                        }
                    }
                    start = i;
                    start_value = value;
                }
            }
        }
    }
    
    let vizrefs = refs.iter().
        map(|(k,v)| (vizchrmap.get(k).unwrap_or(k).clone(), *v)).
        collect::<LinkedHashMap<String,u64>>();
    // write the genome file
    let genome_filename = format!("{}.{}.bw.genome", file, strand);
    {   let mut genome_fh = File::create(&genome_filename)?;
        for (chr, length) in vizrefs {
            writeln!(genome_fh, "{}\t{}", chr, length)?;
        }
    }
    
    // call bedGraphToBigWig
    let bigwig_file = format!("{}.{}.bw", file, strand);
    cmd!("bedGraphToBigWig",&bedgraph_file,&genome_filename,&bigwig_file).run()?;
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
        trackdb.flush()?;
    }
    trackdb.write_fmt(format_args!(r##"
    track {}
    shortLabel {}
    longLabel {}
    type bigWig
    parent {}
    bigDataUrl {}
    "##, track_name, track_name, track_name, parent, url))?;
    trackdb.flush()?;
    
    Ok(())
}

fn get_bam_total_reads(bamfiles: &[String]) -> Result<u64> {
    let mut total_reads = 0u64;
    for bamfile in bamfiles {
        let stdout = cmd!("samtools","idxstats",bamfile).read()?;
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
    pair: &ConstituitivePair,
    mapped_reads: &IntervalTree<u64,String>,
    total_reads: u64) 
    -> Result<RpkmStats>
{
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
        
        let mut cassette_reads = HashMap::<String,Vec<Range<u64>>>::new();
        let mut cassette_features = Vec::<(Range<u64>,f64)>::new();
        let mut intron_reads = HashMap::<String,Vec<Range<u64>>>::new();
        let mut intron_features = Vec::<Range<u64>>::new();
        if pair.cassettes.is_empty() {
            let intron_range = annot.rows[pair.exon1_row].end..(annot.rows[pair.exon2_row].start-1);
            for read in mapped_reads.find(&intron_range) {
                let mut intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                intron_read.push(read.interval().start..read.interval().end);
            }
            intron_features.push(intron_range);
        }
        else {
            let intron_range = annot.rows[pair.exon1_row].end..pair.cassettes[0].range.start;
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
        let rpkmstats = RpkmStats {
            gene_row: *gene_row,
            intron_rpkm: intron_rpkm,
            max_cassette_rpkm: max_exon_rpkm,
            total_constituitive_rpkm: total_constituitive_rpkm,
            total_cassette_rpkm: total_cassette_rpkm,
        };
        return Ok(rpkmstats);
    }
    Err(format!("Could not find gene parents for {:?}", pair).into())
}

fn write_rpkm_stats(
    outfile: &str, 
    annot: &Arc<IndexedAnnotation>, 
    rpkmstats: &mut[RpkmStats]) 
    -> Result<()> 
{
    let mut output: BufWriter<Box<Write>> = BufWriter::new(
        if outfile == "-" { Box::new(stdout()) } 
        else { Box::new(File::create(outfile)?) });
        
    // sort by intron_rpkm / max_exon_rpkm, descending
    let unknown = String::from("unknown");
    rpkmstats.sort_by(|a, b| 
        ((b.intron_rpkm/b.max_cassette_rpkm).is_finite()).
            cmp(&((a.intron_rpkm/a.max_cassette_rpkm).is_finite())).
        then_with(|| OrderedFloat(b.intron_rpkm/b.max_cassette_rpkm).
            cmp(&OrderedFloat(a.intron_rpkm/a.max_cassette_rpkm))).
        then_with(|| {
            let gene1 = &annot.rows[a.gene_row];
            let gene_name1 = 
                gene1.attributes.get("gene_name").or_else(||
                gene1.attributes.get("Name").or_else(||
                gene1.attributes.get("ID"))).unwrap_or(&unknown);
            let pair_name1 = format!("{}:{}:{}..{}:{}", 
                    gene_name1, gene1.seqname, gene1.start-1, gene1.end, gene1.strand);
            let gene2 = &annot.rows[b.gene_row];
            let gene_name2 = 
                gene2.attributes.get("gene_name").or_else(||
                gene2.attributes.get("Name").or_else(||
                gene2.attributes.get("ID"))).unwrap_or(&unknown);
            let pair_name2 = format!("{}:{}:{}..{}:{}", 
                    gene_name2, gene2.seqname, gene2.start-1, gene2.end, gene2.strand);
            pair_name1.cmp(&pair_name2)
        }));
        
    // write the header
    output.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\n", 
        "constituitive_pair_name",
        "intron_rpkm/max_cassette_rpkm",
        "intron_rpkm",
        "max_cassette_rpkm",
        "total_constituitive_rpkm",
        "total_cassette_rpkm",
    ))?;
    for rpkm in rpkmstats {
        let gene = &annot.rows[rpkm.gene_row];
        
        let gene_name = 
            gene.attributes.get("gene_name").or_else(||
            gene.attributes.get("Name").or_else(||
            gene.attributes.get("ID"))).unwrap_or(&unknown);
        let pair_name = format!("{}:{}:{}..{}:{}", 
                gene_name, gene.seqname, gene.start-1, gene.end, gene.strand);
        let ratio = rpkm.intron_rpkm / rpkm.max_cassette_rpkm;
        if !ratio.is_finite() { continue }
        
        output.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\n", 
            pair_name, 
            ratio,
            rpkm.intron_rpkm, 
            rpkm.max_cassette_rpkm, 
            rpkm.total_constituitive_rpkm, 
            rpkm.total_cassette_rpkm))?;
    }
    Ok(())
}

fn get_pair_name(pair: &ConstituitivePair, annot: &IndexedAnnotation) -> String {
    let mut gene_id = None;
    'FIND_GENE_ID:
    for transcript_row in &annot.row2parents[&pair.exon1_row] {
        for gene_row in &annot.row2parents[transcript_row] {
            if annot.rows[*gene_row].feature_type == *"gene" {
                if let Some(gene_name_attr) = annot.rows[*gene_row].attributes.get("gene_name") {
                    gene_id = Some(gene_name_attr);
                    break 'FIND_GENE_ID;
                } else if let Some(gene_name_attr) = annot.rows[*gene_row].attributes.get("Name") {
                    gene_id = Some(gene_name_attr);
                    break 'FIND_GENE_ID;
                } else if let Some(gene_id_attr) = annot.rows[*gene_row].attributes.get("ID") {
                    gene_id = Some(gene_id_attr);
                    break 'FIND_GENE_ID;
                }
            }
        }
        if let Some(transcript_name_attr) = annot.rows[*transcript_row].attributes.get("transcript_name") {
            gene_id = Some(transcript_name_attr);
            break 'FIND_GENE_ID;
        } else if let Some(transcript_name_attr) = annot.rows[*transcript_row].attributes.get("Name") {
            gene_id = Some(transcript_name_attr);
            break 'FIND_GENE_ID;
        } else if let Some(transcript_id_attr) = annot.rows[*transcript_row].attributes.get("ID") {
            gene_id = Some(transcript_id_attr);
            break 'FIND_GENE_ID;
        }
    }
    gene_id = gene_id.or_else(||
        annot.rows[pair.exon1_row].attributes.get("gene_name").or_else(||
        annot.rows[pair.exon1_row].attributes.get("gene_id").or_else(||
        annot.rows[pair.exon1_row].attributes.get("transcript_name").or_else(||
        annot.rows[pair.exon1_row].attributes.get("transcript_id").or_else(||
        annot.rows[pair.exon1_row].attributes.get("Name").or_else(|| 
        annot.rows[pair.exon1_row].attributes.get("ID")))))));
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
    annot: Arc<IndexedAnnotation>,
    file: &str, 
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<()> 
{
    // write the bed file
    let bed_file = format!("{}.bed", file);
    
    {   let mut bw = BufWriter::new(File::create(&bed_file)?);
        for pair in pairs {
            // find the gene ID
            let name = get_pair_name(pair, &annot);
            // write the bed record
            let score = 0;
            let item_rgb = [0,0,0].iter().map(|v| v.to_string()).join(",");
            let mut block_sizes = Vec::<u64>::new();
            let mut block_starts = Vec::<u64>::new();
            block_sizes.push(annot.rows[pair.exon1_row].end-annot.rows[pair.exon1_row].start+1);
            block_starts.push(0u64);
            for cassette in &pair.cassettes {
                block_sizes.push(cassette.range.end-cassette.range.start);
                block_starts.push(cassette.range.start-annot.rows[pair.exon1_row].start+1);
            }
            block_sizes.push(annot.rows[pair.exon2_row].end-annot.rows[pair.exon2_row].start+1);
            block_starts.push(annot.rows[pair.exon2_row].start-annot.rows[pair.exon1_row].start);
            let chr = &annot.rows[pair.exon1_row].seqname;
            let chr = annot.vizchrmap.get(chr).unwrap_or(&chr);
            let thick_start = if pair.cassettes.is_empty() { (annot.rows[pair.exon1_row].start-1).to_string() } 
                else { pair.cassettes[0].range.start.to_string() };
            let thick_end = if pair.cassettes.is_empty() { (annot.rows[pair.exon1_row].start-1).to_string() } 
                else { pair.cassettes[pair.cassettes.len()-1].range.end.to_string() };
            let line = &[
                &chr,
                &(annot.rows[pair.exon1_row].start-1).to_string(),
                &annot.rows[pair.exon2_row].end.to_string(),
                &name,
                &score.to_string(),
                &annot.rows[pair.exon2_row].strand,
                &thick_start,
                &thick_end,
                &item_rgb,
                &(pair.cassettes.len()+2).to_string(),
                &block_sizes.iter().map(|v| v.to_string()).join(","),
                &block_starts.iter().map(|v| v.to_string()).join(","),
            ].iter().join("\t");
            writeln!(bw, "{}", line)?;
        }
    }
    let vizrefs = annot.refs.iter().
        map(|(k,v)| (annot.vizchrmap.get(k).unwrap_or(k).clone(), *v)).
        collect::<LinkedHashMap<String,u64>>();
    bed2bigbed(&bed_file,file,&vizrefs,true,true,trackdb)?;
    Ok(())
}

fn write_enriched_annotation(
    annot: &Arc<IndexedAnnotation>, 
    reannotated_pairs: &Vec<ConstituitivePair>,
    outannot: &str,
    options: &Options,
    trackdb: &mut BufWriter<Box<Write>>) 
    -> Result<()> 
{
    annot.to_gff(&outannot)?; 
    {
        let mut output = OpenOptions::new().append(true).open(outannot)?;
        // get the set of transcript -> pair associations
        let mut transcript2pair = HashMap::<usize,HashSet<usize>>::new();
        for (i, pair) in reannotated_pairs.iter().enumerate() {
            if !pair.cassettes.is_empty() {
                if let Some(ref exon1parents) = annot.row2parents.get(&pair.exon1_row) {
                    if let Some(ref exon2parents) = annot.row2parents.get(&pair.exon2_row) {
                        let exon2parents = exon2parents.iter().collect::<HashSet<_>>();
                        for transcript_row in exon1parents.iter() {
                            if exon2parents.contains(transcript_row) {
                                transcript2pair.entry(*transcript_row).or_insert_with(HashSet::new).insert(i);
                            }
                        }
                    }
                }
            }
        }
        
        // keep track of unique transcript identifiers
        lazy_static! {
            static ref TRANSCRIPT_RENAME: Regex = Regex::new(r"(?:\.([0-9]+))?$").unwrap();
        }
        let mut ids = HashSet::<String>::new();
        let mut names = HashSet::<String>::new();
        // sort transcripts by largest number of associated constituitive pairs
        let mut transcript_order = transcript2pair.keys().collect::<Vec<_>>();
        transcript_order.sort_by(|a,b| transcript2pair[b].len().cmp(&transcript2pair[a].len()));
        let mut seen_pair = HashSet::<usize>::new();
        for transcript_row in transcript_order {
            let transcript = &annot.rows[*transcript_row];
            let mut records = Vec::<Record>::new();
            
            // make sure there is at least one never-before-seen pair that this transcript owns
            let mut any_new_pairs = false;
            for pair_row in &transcript2pair[transcript_row] {
                if !seen_pair.contains(pair_row) {
                    any_new_pairs = true;
                }
            }
            if !any_new_pairs { continue }
                
            // find unique transcript IDs and names
            let mut transcript_id = transcript.attributes.get("transcript_id").map(|s| format!("{}.reannot", s));
            let mut transcript_name = transcript.attributes.get("transcript_name").map(|s| format!("{}.reannot", s));
            transcript_id = if let Some(mut transcript_id) = transcript_id {
                while ids.contains(&transcript_id) {
                    transcript_id = String::from(TRANSCRIPT_RENAME.replace(&transcript_id.as_ref(), |caps: &Captures| {
                        format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                }
                ids.insert(transcript_id.clone());
                Some(transcript_id)
            } else {transcript_id};
            transcript_name = if let Some(mut transcript_name) = transcript_name {
                while ids.contains(&transcript_name) {
                    transcript_name = String::from(TRANSCRIPT_RENAME.replace(&transcript_name.as_ref(), |caps: &Captures| {
                        format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                }
                names.insert(transcript_name.clone());
                Some(transcript_name)
            } else {transcript_name};
            // create a new transcript record
            let mut new_transcript = transcript.clone();
            if let Some(ref transcript_id) = transcript_id {
                {   let mut entry = new_transcript.attributes.entry("ID".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_id.clone());
                }
                {   let mut entry = new_transcript.attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_id.clone());
                }
            }
            if let Some(ref transcript_name) = transcript_name {
                {   let mut entry = new_transcript.attributes.entry("Name".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_name.clone());
                }
                {   let mut entry = new_transcript.attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_name.clone());
                }
            }
            // write new transcript record to file
            records.push(new_transcript);
            
            // get a list of feature starts/stops, and a tree of CDS features
            let mut featurestarts = HashSet::<(String,u64)>::new();
            let mut featurestops = HashSet::<(String,u64)>::new();
            let mut cdstree = IntervalTree::<u64,usize>::new();
            // for each child of the original transcript
            if let Some(ref transcript_row2children) = annot.row2children.get(transcript_row) {
                for child_row in transcript_row2children.iter() {
                    // make a copy
                    let mut child = annot.rows[*child_row].clone();
                    // populate featurestarts/stops and cdstree
                    featurestarts.insert((child.feature_type.clone(), child.start-1));
                    featurestops.insert((child.feature_type.clone(), child.end));
                    if child.feature_type == "CDS" {
                        cdstree.insert(Interval::new((child.start - 1)..(child.end))?, *child_row);
                    }
                    // update the parents list with the new transcript ID
                    let mut newparents = Vec::<String>::new();
                    if let Some(parent) = child.attributes.get("Parent") {
                        if let Some(oldid) = transcript.attributes.get("ID") {
                            if let Some(ref id) = transcript_id {
                                for p in parent.split(",") {
                                    newparents.push(if p == oldid {id.clone()} else {p.to_string()});
                                }
                            }
                        }
                    }
                    
                    {   let mut entry = child.attributes.entry("Parent".to_string()).or_insert_with(|| "".to_string());
                        entry.clear();
                        entry.push_str(&mut newparents.join(","));
                    }
                    // update transcript_id and transcript_name attributes
                    if child.attributes.contains_key("transcript_id") {
                        if let Some(ref transcript_id) = transcript_id {
                            let mut entry = child.attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_id.clone());
                        }
                    }
                    if child.attributes.contains_key("transcript_name") {
                        if let Some(ref transcript_name) = transcript_name {
                            let mut entry = child.attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_name.clone());
                        }
                    }
                    // store child record
                    records.push(child);
                    
                    // update any grandchildren with correct transcript_ids and transcript_names
                    let mut seen = HashSet::<usize>::new();
                    let mut children = Vec::<usize>::new();
                    if let Some(ref mut cs) = annot.row2children.get(child_row) {
                        children.append(&mut cs.clone());
                    }
                    while !children.is_empty() {
                        let mut newchildren = Vec::<usize>::new();
                        for c in children {
                            if !seen.contains(&c)  {
                                seen.insert(c);
                                let mut cc = annot.rows[c].clone();
                                featurestarts.insert((cc.feature_type.clone(), cc.start-1));
                                featurestops.insert((cc.feature_type.clone(), cc.end));
                                if cc.feature_type == "CDS" {
                                    cdstree.insert(Interval::new((cc.start - 1)..(cc.end))?, c);
                                }
                                // update transcript_id and transcript_name attributes
                                if cc.attributes.contains_key("transcript_id") {
                                    if let Some(ref transcript_id) = transcript_id {
                                        let mut entry = cc.attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                                        entry.clear();
                                        entry.push_str(&mut transcript_id.clone());
                                    }
                                }
                                if cc.attributes.contains_key("transcript_name") {
                                    if let Some(ref transcript_name) = transcript_name {
                                        let mut entry = cc.attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                                        entry.clear();
                                        entry.push_str(&mut transcript_name.clone());
                                    }
                                }
                                // write child record
                                records.push(cc);
                                // add more children
                                if let Some(ref mut cs) = annot.row2children.get(&c) {
                                    newchildren.append(&mut cs.clone());
                                }
                            }
                        }
                        children = newchildren;
                    }
                }
            }
            
            for pair_row in &transcript2pair[transcript_row] {
                if !seen_pair.contains(pair_row) {
                    seen_pair.insert(*pair_row);
                    let pair = &reannotated_pairs[*pair_row];
                    let exon1 = &annot.rows[pair.exon1_row];
                    let exon2 = &annot.rows[pair.exon2_row];
                        
                    let cdstype = 
                        if featurestops.contains(&("CDS".to_string(), exon1.end)) && 
                            featurestarts.contains(&("CDS".to_string(), exon2.start-1))
                        { Some(("CDS","CDS")) }
                        else if featurestops.contains(&("five_prime_UTR".to_string(), exon1.end)) && 
                            featurestarts.contains(&("five_prime_UTR".to_string(), exon2.start-1))
                        { Some(("five_prime_UTR","UTR5")) }
                        else if featurestops.contains(&("three_prime_UTR".to_string(), exon1.end)) && 
                            featurestarts.contains(&("three_prime_UTR".to_string(), exon2.start-1))
                        { Some(("three_prime_UTR","UTR3")) }
                        else if featurestops.contains(&("UTR".to_string(), exon1.end)) && 
                            featurestarts.contains(&("UTR".to_string(), exon2.start-1))
                        { Some(("UTR","UTR")) }
                        else { None };
                    
                    // finally, add the reannotated cassettes
                    lazy_static! {
                        static ref EXON_RENAME: Regex = Regex::new(r"(?:[:]([0-9]+))?$").unwrap();
                    }
                    for cassette in &pair.cassettes {
                        // get a unique ID for the cassette
                        let mut cassette_id = transcript.attributes.get("ID").map(|s| format!("exon:{}:1", s));
                        cassette_id = if let Some(mut cassette_id) = cassette_id {
                            while ids.contains(&cassette_id) {
                                cassette_id = String::from(EXON_RENAME.replace(&cassette_id.as_ref(), |caps: &Captures| {
                                    format!(":{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                            }
                            ids.insert(cassette_id.clone());
                            Some(cassette_id)
                        } else {cassette_id};
                        
                        let mut attributes = exon1.attributes.clone();
                        if let Some(cassette_id) = cassette_id {
                            let mut entry = attributes.entry("ID".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut cassette_id.clone());
                        }
                        if let Some(ref transcript_id) = transcript_id {
                            let mut entry = attributes.entry("Parent".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_id.clone());
                        }
                        if let Some(ref transcript_id) = transcript_id {
                            let mut entry = attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_id.clone());
                        }
                        if let Some(ref transcript_name) = transcript_name {
                            let mut entry = attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_name.clone());
                        }
                        
                        let record = Record {
                            row: 0,
                            seqname: exon1.seqname.clone(),
                            source: exon1.source.clone(),
                            feature_type: "exon".to_string(),
                            start: cassette.range.start+1,
                            end: cassette.range.end,
                            score: ".".to_string(),
                            strand: exon1.strand.clone(),
                            frame: ".".to_string(),
                            attributes: attributes.clone(),
                        };
                        records.push(record);
                        if let Some(cdstype) = cdstype {
                            let mut attributes = attributes.clone();
                            // get a unique ID for the cassette CDS feature
                            let mut cds_id = transcript.attributes.get("ID").map(|s| format!("{}:{}:1", cdstype.1, s));
                            cds_id = if let Some(mut cds_id) = cds_id {
                                while ids.contains(&cds_id) {
                                    cds_id = String::from(EXON_RENAME.replace(&cds_id.as_ref(), |caps: &Captures| {
                                        format!(":{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                                }
                                ids.insert(cds_id.clone());
                                Some(cds_id)
                            } else {cds_id};
                            if let Some(cds_id) = cds_id {
                                let mut entry = attributes.entry("ID".to_string()).or_insert_with(|| "".to_string());
                                entry.clear();
                                entry.push_str(&mut cds_id.clone());
                            }
                            let record = Record {
                                row: 0,
                                seqname: exon1.seqname.clone(),
                                source: exon1.source.clone(),
                                feature_type: cdstype.0.to_string(),
                                start: cassette.range.start+1,
                                end: cassette.range.end,
                                score: ".".to_string(),
                                strand: exon1.strand.clone(),
                                frame: ".".to_string(),
                                attributes: attributes,
                            };
                            records.push(record);
                        }
                    }
                }
            }
            // compute frame field for all the CDS features
            // write features to output file
            for record in records {
                let mut record = record;
                // compute the frame for CDS features
                if record.feature_type == "CDS" {
                    let prevcdslen: u64 = cdstree.find(
                            if transcript.strand == "-" {record.end..std::u64::MAX}
                            else {0..(record.start-1)}).
                        map(|c| annot.rows[*c.data()].end-annot.rows[*c.data()].start+1).sum();
                    let frame = (prevcdslen % 3).to_string();
                    record.frame = frame;
                }
                writeln!(output, "{}", record.to_gff()?)?;
            }
        }
    }
    
    if let Some(ref debug_outannot_bigbed) = options.debug_outannot_bigbed {
        let transcript_type = String::from("transcript");
        let newannot = IndexedAnnotation::from_gff(&outannot, 
            options.gene_type.get(0).r()?, 
            options.transcript_type.get(0).unwrap_or(&transcript_type), 
            &options.chrmap_file,
            &options.vizchrmap_file)?;
        newannot.to_bigbed(
            &debug_outannot_bigbed, 
            &options.exon_type, 
            &options.cds_type, 
            &options.transcript_type,
            &options.gene_type,
            trackdb)?; 
    }
    Ok(())
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    // set defaults for feature types
    options.gene_type =
        (if options.gene_type.is_empty() { vec!["gene".to_string()] } 
         else { options.gene_type.clone() }).into_iter().collect();
    //options.transcript_type = 
    //    (if options.transcript_type.is_empty() { vec!["transcript".to_string()] } 
    //     else { options.transcript_type.clone() }).into_iter().collect();
    options.exon_type =
        (if options.exon_type.is_empty() { vec!["exon".to_string()] }
         else { options.exon_type.clone() }).into_iter().collect();
    options.cds_type =
        (if options.cds_type.is_empty() { vec!["CDS".to_string()] }
         else { options.cds_type.clone() }).into_iter().collect();
    // set debug options if --debug flag is set
    if options.debug {
        if options.debug_annot_gff.is_none() { options.debug_annot_gff = Some(format!("{}.annot.gff", options.debug_prefix)) }
        // if options.debug_annot_gtf.is_none() { options.debug_annot_gtf = Some(format!("{}.annot.gtf", options.debug_prefix)) }
        if options.debug_annot_bigbed.is_none() { options.debug_annot_bigbed = Some(format!("{}.annot.bb", options.debug_prefix)) }
        // if options.debug_annot_json.is_none() { options.debug_annot_json = Some(format!("{}.annot.json", options.debug_prefix)) }
        if options.debug_exon_bigbed.is_none() { options.debug_exon_bigbed = Some(format!("{}.exon.bb", options.debug_prefix)) }
        if options.debug_bigwig.is_none() { options.debug_bigwig = Some(format!("{}.intronrpkm", options.debug_prefix)) }
        if options.debug_reannot_bigbed.is_none() { options.debug_reannot_bigbed = Some(format!("{}.reannot.bb", options.debug_prefix)) }
        if options.debug_rpkm_region_bigbed.is_none() { options.debug_rpkm_region_bigbed = Some(format!("{}.rpkm.bb", options.debug_prefix)) }
        if options.debug_rpkmstats_json.is_none() { options.debug_rpkmstats_json = Some(format!("{}.rpkmstats.json", options.debug_prefix)) }
        if options.debug_trackdb.is_none() { options.debug_trackdb = Some(format!("{}.trackDb.txt", options.debug_prefix)) }
        if options.debug_outannot_bigbed.is_none() { options.debug_outannot_bigbed = Some(format!("{}.outannot.bb", options.debug_prefix)) }
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
    let transcript_type = String::from("transcript");
    let mut annot = if let Some(annotfile_gff) = options.annotfile_gff.clone() {
        writeln!(stderr(), "Reading annotation file {:?}", &annotfile_gff)?;
        IndexedAnnotation::from_gff(
            &annotfile_gff, 
            options.gene_type.get(0).r()?, 
            options.transcript_type.get(0).unwrap_or(&transcript_type),
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else if let Some(annotfile_gtf) = options.annotfile_gtf.clone() {
        writeln!(stderr(), "Reading annotation file {:?}", &annotfile_gtf)?;
        IndexedAnnotation::from_gtf(&annotfile_gtf, 
            options.gene_type.get(0).r()?, 
            options.transcript_type.get(0).unwrap_or(&transcript_type),
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
    
    if let Some(ref debug_annot_gff) = options.debug_annot_gff {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_gff)?;
        annot.to_gff(&debug_annot_gff)?; 
    }
    if let Some(ref debug_annot_gtf) = options.debug_annot_gtf {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_gtf)?;
        annot.to_gtf(&debug_annot_gtf)?; 
    }
    if let Some(ref debug_annot_bigbed) = options.debug_annot_bigbed {
        writeln!(stderr(), "Writing annotation file to {:?}", &debug_annot_bigbed)?;
        annot.to_bigbed(
            &debug_annot_bigbed, 
            &options.exon_type, 
            &options.cds_type, 
            &options.transcript_type,
            &options.gene_type,
            &mut trackdb)?; 
    }
    if let Some(ref debug_annot_json) = options.debug_annot_json {
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
    
    let annot = Arc::new(annot);
    if let Some(ref debug_exon_bigbed) = options.debug_exon_bigbed {
        writeln!(stderr(), "Writing constituitive pairs to a bigbed file")?;
        write_exon_bigbed(&exonpairs, annot.clone(), &debug_exon_bigbed, &mut trackdb)?;
    }
        
    writeln!(stderr(), "Reannotating cassette regions")?;
    let (reannotated_pairs, mut rpkmstats) = reannotate_regions(
        &annot,
        &exonpairs, 
        &bamfiles, 
        &bamstrand,
        total_reads,
        &options,
        &mut trackdb)?;
    if let Some(ref debug_reannot_bigbed) = options.debug_reannot_bigbed {
        writeln!(stderr(), "Writing reannotation to bigbed")?;
        write_exon_bigbed(&reannotated_pairs, annot.clone(), &debug_reannot_bigbed, &mut trackdb)?;
    }
    
    writeln!(stderr(), "Computing RPKM stats")?;
    if let Some(ref debug_rpkmstats_json) = options.debug_rpkmstats_json {
        writeln!(stderr(), "Writing RPKM stats to {:?}", &debug_rpkmstats_json)?;
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if debug_rpkmstats_json == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(debug_rpkmstats_json)?) });
            serde_json::to_writer_pretty(&mut output, &rpkmstats)?;
    }
    writeln!(stderr(), "Writing RPKM stats to {:?}", &options.outfile)?;
    write_rpkm_stats(&options.outfile, &annot, &mut rpkmstats)?;
    
    if let Some(ref outannot) = options.outannot {
        writeln!(stderr(), "Writing output annotation file")?;
        write_enriched_annotation(&annot, &reannotated_pairs, &outannot, &options, &mut trackdb)?;
    }
    Ok(())
}

fn main() {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "full");

    if let Err(ref e) = run() {
        writeln!(stderr(), "error: {}", e).unwrap();
        for e in e.iter().skip(1) {
            writeln!(stderr(), "caused by: {}", e).unwrap();
        }
        if let Some(backtrace) = e.backtrace() {
            writeln!(stderr(), "backtrace: {:?}", backtrace).unwrap();
        }
        ::std::process::exit(1);
    }
}
