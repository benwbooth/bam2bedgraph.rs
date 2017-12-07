#![recursion_limit="128"]
#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::str;
use std::vec::Vec;
use std::collections::HashMap;
use std::collections::HashSet;
use std::ops::Range;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::io::{stdout, sink};
use std::path::Path;
use std::sync::Arc;
use std::fs::File;

#[macro_use] 
extern crate failure;

extern crate bam2bedgraph;
use bam2bedgraph::*;
use bam2bedgraph::error::*;
use bam2bedgraph::power_set::*;
use bam2bedgraph::indexed_annotation::*;

#[macro_use] 
extern crate lazy_static;

#[macro_use]
extern crate url;
use url::percent_encoding::{utf8_percent_encode, SIMPLE_ENCODE_SET};

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

extern crate serde;
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

extern crate unindent;
use unindent::unindent;

extern crate linked_hash_map;
use linked_hash_map::LinkedHashMap;

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
    #[structopt(long="genome", help = "The input genome FASTA file", name="GENOME_FASTA_FILE")]
    genome_file: Option<String>,
    
    // flags
    #[structopt(long="max_iterations", help = "How many start/stop combinations before we skip this one?", name="MAX_ITERATIONS", default_value="1000000")]
    max_iterations: usize,
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
    #[structopt(long="debug_outannot_gtf", help = "Write the output annotation as a gtf file", name="DEBUG_OUTANNOT_GTF_FILE")]
    debug_outannot_gtf: Option<String>,
    #[structopt(long="debug_outannot_fasta", help = "Write the output annotation as a fasta file", name="DEBUG_OUTANNOT_FASTA_FILE")]
    debug_outannot_fasta: Option<String>,
    #[structopt(long="debug_retained_introns", help = "Write the retained introns to a file", name="DEBUG_RETAINED_INTRONS")]
    debug_retained_introns: Option<String>,
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
    is_retained_intron: bool,
}
impl std::fmt::Debug for ConstituitivePair {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "exon1_row: {}, exon2_row: {}, cassettes: {:?}, is_retained_intron: {}", 
            self.exon1_row, self.exon2_row, self.cassettes, self.is_retained_intron)
    }    
}

fn find_constituitive_splice_pairs(annot: &IndexedAnnotation,
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
                let mut transcript_rows = HashSet::<usize>::new();
                let mut transcript2exon = HashMap::<usize,HashSet<usize>>::new();
                let mut start2transcript = HashMap::<u64,HashSet<usize>>::new();
                let mut end2transcript = HashMap::<u64,HashSet<usize>>::new();
                let mut transcript_tree = IntervalTree::<u64, usize>::new();
                // splice start..splice end -> vec![(transcript_row, exon1_row, exon2_row)]
                let mut splices = HashMap::<Range<u64>,Vec<(usize,usize,usize)>>::new();
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
                        if let Some(exon_rows) = annot.row2children.get(transcript_row) {
                            // first get exon start/stop -> transcript associations
                            // and splice start/stop -> transcript associations
                            for exon_row in exon_rows {
                                let exon = &annot.rows[*exon_row];
                                if exon_types.is_empty() || exon_types.contains(&exon.feature_type) {
                                    // skip any transcripts with weird exons
                                    if exon.seqname != transcript.seqname || exon.strand != transcript.strand {
                                        continue 'transcript_row;
                                    }
                                    // make sure the transcript has at least 1 exon
                                    if !transcript_rows.contains(transcript_row) {
                                        transcript_rows.insert(*transcript_row);
                                        transcript_tree.insert(
                                            Interval::new(transcript.start-1..transcript.end)?, 
                                            *transcript_row);
                                    }
                                    transcript2exon.entry(*transcript_row).or_insert_with(HashSet::new).
                                        insert(*exon_row);
                                }
                            }
                        }
                    }
                }
                // record each transcript's splice starts/ends
                for transcript_row in &transcript_rows {
                    if let Some(exon_rows) = transcript2exon.get(&transcript_row) {
                        let mut exon_rows = exon_rows.iter().collect::<Vec<_>>();
                        exon_rows.sort_by_key(|a| &annot.rows[**a].start);
                        for exon_row in exon_rows {
                            let exon = &annot.rows[*exon_row];
                            start2transcript.entry(exon.start-1).or_insert_with(HashSet::new).
                                insert(*transcript_row);
                            end2transcript.entry(exon.end).or_insert_with(HashSet::new).
                                insert(*transcript_row);
                        }
                    }
                }
                // find constituitive splice pairs
                for transcript_row in &transcript_rows {
                    if let Some(exon_rows) = transcript2exon.get(&transcript_row) {
                        let mut exon_rows = exon_rows.iter().collect::<Vec<_>>();
                        exon_rows.sort_by_key(|a| &annot.rows[**a].start);
                        let mut exon1_row: Option<usize> = None;
                        for exon_row in exon_rows {
                            let exon = &annot.rows[*exon_row];
                            if let Some(exon1_row_) = exon1_row {
                                let exon1 = &annot.rows[exon1_row_];
                                let containing_trs = transcript_tree.
                                    find((exon.start-1)..exon.start).
                                    map(|t| *t.data()).
                                    collect::<HashSet<_>>();
                                if start2transcript[&(exon.start-1)].is_superset(&containing_trs) {
                                    splices.entry(exon1.end..(exon.start-1)).
                                        or_insert_with(Vec::new).
                                        push((*transcript_row, exon1_row_, *exon_row));
                                    exon1_row = None;
                                    let containing_trs = transcript_tree.
                                        find(exon.end..exon.end+1).
                                        map(|t| *t.data()).
                                        collect::<HashSet<_>>();
                                    if end2transcript[&(exon.end)].is_superset(&containing_trs) {
                                        exon1_row = Some(*exon_row);
                                    }
                                }
                            }
                            else {
                                let containing_trs = transcript_tree.
                                    find(exon.end..exon.end+1).
                                    map(|t| *t.data()).
                                    collect::<HashSet<_>>();
                                if end2transcript[&(exon.end)].is_superset(&containing_trs) {
                                    exon1_row = Some(*exon_row);
                                }
                            }
                        }
                    }
                }
                for (splice_range, splices) in splices {
                    let containing_trs = transcript_tree.
                        find(&splice_range).
                        map(|t| *t.data()).
                        filter(|t| annot.rows[*t].start-1 < splice_range.start && 
                            splice_range.end < annot.rows[*t].end).collect::<HashSet<_>>();
                    let splice_trs = splices.iter().map(|&(transcript_row,_,_)| transcript_row).collect::<HashSet<_>>();
                    // look at constituitive splices
                    if splice_trs.is_superset(&containing_trs) {
                        let mut exon_rows = HashSet::<(usize,usize)>::new();
                        for &(_, exon1_row, exon2_row) in &splices {
                            exon_rows.insert((exon1_row, exon2_row));
                        }
                        // sort by shortest sum of exon lengths
                        let mut exon_rows = exon_rows.iter().collect::<Vec<_>>();
                        exon_rows.sort_by(|a,b| 
                            ((&annot.rows[a.0].end-&annot.rows[a.0].start+1) + (&annot.rows[a.1].end-&annot.rows[a.1].start+1)).
                            cmp(&((&annot.rows[b.0].end-&annot.rows[b.0].start+1) + (&annot.rows[b.1].end-&annot.rows[b.1].start+1))).
                            // then sort by lowest exon1_row 
                            then_with(|| a.0.cmp(&b.0)).
                            // then sort by lowest exon2_row
                            then_with(|| a.1.cmp(&b.1)));
                        if let Some(exon_row) = exon_rows.get(0) {
                            let exon1 = &annot.rows[exon_row.0];
                            let exon2 = &annot.rows[exon_row.1];
                            let region = (exon2.start-1)-exon1.end;
                            if region == 0 {
                                eprintln!("Region between exon1 (row {}) and exon2 (row {}) is zero, skipping constituitive pair", exon_row.0, exon_row.1);
                            }
                            else {
                                exonpairs.push(ConstituitivePair {
                                    exon1_row: exon_row.0,
                                    exon2_row: exon_row.1,
                                    cassettes: Vec::new(),
                                    is_retained_intron: false,
                                });
                            }
                        }
                    }
                }
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
    {   let mut genome_fh = BufWriter::new(File::create(&genome_filename)?);
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
    trackdb.write_fmt(format_args!("{}", unindent(&format!(r##"
    
        track {}
        type bigBed 12
        bigDataUrl {}
        shortLabel {}
        longLabel {}
        itemRgb On
        visibility hide
    "##, track_name, url, track_name, track_name))))?;
    trackdb.flush()?;
        
    Ok(())
}

fn reannotate_pair(
    pair_name: &str,
    exon1: &Record,
    exon2: &Record,
    read_pairs: &HashMap<(String,String),Vec<Vec<Range<u64>>>>,
    debug_bigwig: &Option<String>,
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
    let mut read_coverage = vec![0i32; region_size];
    for (&(_,ref read_name),read_pair) in read_pairs {
        // fill in the read_coverage histogram
        for exons in read_pair {
            for exon in exons {
                if start < (exon.end as usize) && (exon.start as usize) < end {
                    for pos in std::cmp::max(start as u64, exon.start)..std::cmp::min(end as u64, exon.end) {
                        read_coverage[pos as usize - start] += 1;
                    }
                }
            }
        }
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
    const MIN_READ_COUNT_PER_BASE: i32 = 1;
    for (i, value) in histo.iter().enumerate() {
        if (*value == 0) != (exon_value == 0) {
            if exon_value >= MIN_READ_COUNT_PER_BASE {
                let exon_end = i+start;
                exon_regions.push(exon_start..exon_end);
            }
            exon_start = i+start;
            exon_value = *value;
        }
    }
    const MIN_STARTS: i32 = 2;
    const MIN_STOPS: i32 = 2;
    'EXON_REGION:
    for exon_region in exon_regions {
        let mut starts = Vec::<(usize,i32)>::new();
        for (i, value) in start_histo[(exon_region.start-start)..(exon_region.end-start)].iter().enumerate() {
            if *value >= MIN_STARTS {
                starts.push((i+exon_region.start,*value));
            }
        }
        let mut ends = Vec::<(usize,i32)>::new();
        for (i, value) in end_histo[(exon_region.start-start)..(exon_region.end-start)].iter().enumerate() {
            if *value >= MIN_STOPS {
                ends.push((i+exon_region.start,*value));
            }
        }
        // get the set of possible starts/stops
        let mut pairs = HashMap::<Range<usize>,(i32,i32)>::new();
        let mut iterations = 0;
        for &(s, sscore) in &starts {
            for &(e, escore) in &ends {
                if s < e {
                    pairs.insert(s..(e+1), (sscore, escore));
                }
                iterations += 1;
                if iterations > max_iterations {
                    eprintln!("More than {} iterations on pair {}", max_iterations, pair_name);
                    continue 'EXON_REGION;
                }
            }
        }
        if pairs.len() >= 64 {
            eprintln!("Too many possible cassette start/stops for pair {}", pair_name);
            continue 'EXON_REGION;
        }
        // iterate over the powerset of the start/stop set
        let mut best_set_score = 0;
        let mut best_set = None;
        let mut iterations = 0;
        'SET:
        for set in PowerSet::new(&pairs.keys().collect::<Vec<_>>()) {
            let mut set_score = 0;
            let mut overlaps = IntervalTree::<usize,()>::new();
            let mut set = set.clone();
            set.sort_by_key(|a| a.start);
            for pair in &set {
                // make sure there are no overlaps in this set
                for _ in overlaps.find(pair.start..pair.end) { continue 'SET; }
                overlaps.insert(Interval::new(pair.start..pair.end)?, ());
                // add to the set score
                let (sscore, escore) = pairs[pair];
                set_score += sscore + escore;
                // make sure we don't go over the max iterations
                iterations += 1;
                if iterations > max_iterations {
                    eprintln!("More than {} iterations on pair {}", max_iterations, pair_name);
                    continue 'EXON_REGION;
                }
            }
            if best_set.is_none() || set_score > best_set_score {
                best_set = Some(set);
                best_set_score = set_score;
            }
        }
        // store the reannotated cassette exon
        if let Some(best_set) = best_set {
            for pair in best_set {
                cassettes.push(Cassette {
                    range: pair.start as u64..pair.end as u64,
                    cassette_row: None,
                });
            }
        }
    }
    let mut cassette_coverage = vec![0i32; region_size];
    for cassette in &cassettes {
        for pos in cassette.range.start..cassette.range.end {
            cassette_coverage[pos as usize-start] += 1;
        }
    }
    let mut total_bases = 0u64;
    let mut covered_bases = 0u64;
    for (i, read_ct) in read_coverage.iter().enumerate() {
        if cassette_coverage[i] == 0 {
            total_bases += 1;
            if *read_ct > 0 { covered_bases += 1 }
        }
    }
    let base_coverage = covered_bases as f64 / total_bases as f64;
    let is_retained_intron = base_coverage >= 0.9;
    let reannotpair = ConstituitivePair {
        exon1_row: exon1.row,
        exon2_row: exon2.row,
        cassettes: cassettes,
        is_retained_intron: is_retained_intron,
    };
    //eprintln!("Writing reannotated pair: {:?}", reannotpair);
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
                    bam.fetch(*tid, start as u32, end as u32)?;
                    for read in bam.records() {
                        let read = read?;
                        // make sure the read's strand matches
                        if let Some(is_read1strand) = read1strand {
                            if !(((is_read1strand == read.is_first_in_template()) == !read.is_reverse()) == strand_is_plus) {
                                continue;
                            }
                        }
                        let read_name = String::from(str::from_utf8(read.qname())?);
                        //eprintln!("Looking at read: {}", read_name);
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
                eprintln!("Got Err in reannotation: {:?}", e);
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
    {   let mut genome_fh = BufWriter::new(File::create(&genome_filename)?);
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
        trackdb.write_fmt(format_args!("{}", unindent(&format!(r##"

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
            "##, parent, parent, parent, viewlimits))))?;
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

#[derive(Serialize)]
struct RpkmStats {
    pair_name: String,
    intron_rpkm: f64,
    max_cassette_rpkm: f64,
    cassette_cov: f64,
    exon1_rpkm: f64,
    exon2_rpkm: f64,
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
    // rpkm = (10^9 * num_mapped_reads_to_target)/(total_mapped_reads * mappable_target_size_bp)
    let exon1_bases = annot.rows[pair.exon1_row].end-annot.rows[pair.exon1_row].start+1;
    let exon2_bases = annot.rows[pair.exon2_row].end-annot.rows[pair.exon2_row].start+1;
    let constituitive_bases = exon1_bases + exon2_bases;
    let mut exon1_reads = HashSet::<String>::new();
    let mut exon2_reads = HashSet::<String>::new();
    for read in mapped_reads.find(annot.rows[pair.exon1_row].start-1..annot.rows[pair.exon1_row].end) {
        exon1_reads.insert(read.data().clone());
    }
    for read in mapped_reads.find(annot.rows[pair.exon2_row].start-1..annot.rows[pair.exon2_row].end) {
        exon2_reads.insert(read.data().clone());
    }
    let constituitive_reads = exon1_reads.len() + exon2_reads.len();
    let exon1_rpkm = if total_reads == 0 || exon1_bases == 0 { 0f64 }
    else {
        (1e10f64 * exon1_reads.len() as f64) 
            / (total_reads as f64 * exon1_bases as f64)
    };
    let exon2_rpkm = if total_reads == 0 || exon2_bases == 0 { 0f64 }
    else {
        (1e10f64 * exon2_reads.len() as f64) 
            / (total_reads as f64 * exon2_bases as f64)
    };
    let total_constituitive_rpkm = if total_reads == 0 || constituitive_bases == 0 { 0f64 } 
    else {
        (1e10f64 * constituitive_reads as f64) 
            / (total_reads as f64 * constituitive_bases as f64)
    };
    
    let mut cassette_cov = 0f64;
    let mut cassette_reads = HashMap::<String,Vec<Range<u64>>>::new();
    let mut cassette_features = Vec::<(Range<u64>,f64)>::new();
    let mut intron_reads = HashMap::<String,Vec<Range<u64>>>::new();
    let mut intron_features = Vec::<Range<u64>>::new();
    if pair.cassettes.is_empty() {
        let intron_range = annot.rows[pair.exon1_row].end..(annot.rows[pair.exon2_row].start-1);
        for read in mapped_reads.find(&intron_range) {
            let intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
            intron_read.push(read.interval().start..read.interval().end);
        }
        intron_features.push(intron_range);
    }
    else {
        let intron_range = annot.rows[pair.exon1_row].end..pair.cassettes[0].range.start;
        for read in mapped_reads.find(&intron_range) {
            let intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
            intron_read.push(read.interval().start..read.interval().end);
        }
        intron_features.push(intron_range);
        
        for (i, cassette) in pair.cassettes.iter().enumerate() {
            let reads: Vec<_> = mapped_reads.find(&cassette.range).collect();
            for read in &reads {
                let cassette_read = cassette_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                cassette_read.push(read.interval().start..read.interval().end);
                cassette_cov +=
                    std::cmp::min(cassette.range.end, read.interval().end) as f64 -
                    std::cmp::max(cassette.range.start, read.interval().start) as f64;
            }
            let exon_rpkm = (1e10f64 * reads.len() as f64) / 
                            (total_reads as f64 * (cassette.range.end-cassette.range.start) as f64);
            cassette_features.push((cassette.range.clone(), exon_rpkm));

            if i < pair.cassettes.len()-1 {
                let intron_range = cassette.range.end..pair.cassettes[i+1].range.start;
                for read in mapped_reads.find(&intron_range) {
                    let intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                    intron_read.push(read.interval().start..read.interval().end);
                }
                intron_features.push(intron_range);
            }
        }
        
        let intron_range = pair.cassettes[pair.cassettes.len()-1].range.start..annot.rows[pair.exon2_row].start;
        for read in mapped_reads.find(&intron_range) {
                let intron_read = intron_reads.entry(read.data().clone()).or_insert_with(Vec::new);
                intron_read.push(read.interval().start..read.interval().end);
        }
        intron_features.push(intron_range);
    }
    let cassette_bases = cassette_features.iter().map(|f| f.0.end-f.0.start).sum::<u64>();
    let intron_bases = intron_features.iter().map(|f| f.end-f.start).sum::<u64>();
    let max_exon_rpkm = cassette_features.iter().map(|f| f.1).fold(std::f64::NAN, f64::max);
    let total_cassette_rpkm = (1e10f64 * cassette_reads.len() as f64) / (total_reads as f64 * cassette_bases as f64);
    let intron_rpkm = (1e10f64 * intron_reads.len() as f64) / (total_reads as f64 * intron_bases as f64);
    let pair_name = get_pair_name(pair, &annot);
    let rpkmstats = RpkmStats {
        pair_name: pair_name,
        intron_rpkm: intron_rpkm,
        max_cassette_rpkm: max_exon_rpkm,
        cassette_cov: cassette_cov / cassette_bases as f64,
        exon1_rpkm: exon1_rpkm,
        exon2_rpkm: exon2_rpkm,
        total_constituitive_rpkm: total_constituitive_rpkm,
        total_cassette_rpkm: total_cassette_rpkm,
    };
    return Ok(rpkmstats);
}

fn write_rpkm_stats(
    outfile: &str, 
    rpkmstats: &mut[RpkmStats]) 
    -> Result<()> 
{
    let mut output: BufWriter<Box<Write>> = BufWriter::new(
        if outfile == "-" { Box::new(stdout()) } 
        else { Box::new(File::create(outfile)?) });
        
    // sort by intron_rpkm / max_exon_rpkm, descending
    rpkmstats.sort_by(|a, b| 
        ((b.intron_rpkm/b.max_cassette_rpkm).is_finite()).
            cmp(&((a.intron_rpkm/a.max_cassette_rpkm).is_finite())).
        then_with(|| OrderedFloat(b.intron_rpkm/b.max_cassette_rpkm).
            cmp(&OrderedFloat(a.intron_rpkm/a.max_cassette_rpkm))).
        then_with(|| a.pair_name.cmp(&b.pair_name)));
        
    // write the header
    output.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
        "constituitive_pair_name",
        "intron_rpkm/max_cassette_rpkm",
        "intron_rpkm",
        "max_cassette_rpkm",
        "cassette_cov",
        "exon1_rpkm",
        "exon2_rpkm",
        "total_constituitive_rpkm",
        "total_cassette_rpkm",
    ))?;
    for rpkm in rpkmstats {
        let ratio = rpkm.intron_rpkm / rpkm.max_cassette_rpkm;
        if !ratio.is_finite() { continue }
        
        output.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            rpkm.pair_name, 
            ratio,
            rpkm.intron_rpkm, 
            rpkm.max_cassette_rpkm, 
            rpkm.cassette_cov,
            rpkm.exon1_rpkm, 
            rpkm.exon2_rpkm, 
            rpkm.total_constituitive_rpkm, 
            rpkm.total_cassette_rpkm))?;
    }
    Ok(())
}

fn get_name(row: usize, annot: &IndexedAnnotation) -> Option<String> {
    let name = annot.rows[row].attributes.get("transcript_name").or_else(||
        annot.rows[row].attributes.get("transcript").or_else(||
        annot.rows[row].attributes.get("Name").or_else(||
        annot.rows[row].attributes.get("ID").or_else(||
        annot.rows[row].attributes.get("transcript_id").or_else(||
        annot.rows[row].attributes.get("gene_name").or_else(||
        annot.rows[row].attributes.get("gene").or_else(||
        annot.rows[row].attributes.get("gene_id"))))))));
    name.map(|n| n.to_string())
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
    let name = format!("{}{}:{}..{}:{}{}", 
        match gene_id {
            Some(g) => format!("{}:", g),
            None => "".to_string(),
        },
        annot.rows[pair.exon1_row].seqname,
        annot.rows[pair.exon1_row].start-1,
        annot.rows[pair.exon2_row].end,
        annot.rows[pair.exon2_row].strand,
        if pair.is_retained_intron {".retained_intron"} else {""});
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

fn write_retained_introns(
    pairs: &[ConstituitivePair], 
    annot: &IndexedAnnotation,
    file: &str)
    -> Result<()> 
{
    let mut bw = BufWriter::new(File::create(&file)?);
    for pair in pairs {
        if pair.is_retained_intron {
            let pair_name = get_pair_name(pair, annot);
            let line = &[
                &pair_name,
                &annot.rows[pair.exon1_row].seqname,
                &annot.rows[pair.exon1_row].end.to_string(),
                &(annot.rows[pair.exon2_row].start-1).to_string(),
                &annot.rows[pair.exon1_row].strand,
            ].iter().join("\t");
            writeln!(bw, "{}", line)?;
        }
    }
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
            static ref EXON_RENAME: Regex = Regex::new(r"(?:[:]([0-9]+))?$").unwrap();
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
            } else {transcript_id.clone()};
            // create a new transcript record
            let mut new_transcript = transcript.clone();
            if let Some(ref transcript_id) = transcript_id {
                {   let entry = new_transcript.attributes.entry("ID".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_id.clone());
                }
                {   let entry = new_transcript.attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_id.clone());
                }
            }
            if let Some(ref transcript_name) = transcript_name {
                {   let entry = new_transcript.attributes.entry("Name".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_name.clone());
                }
                {   let entry = new_transcript.attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                    entry.clear();
                    entry.push_str(&mut transcript_name.clone());
                }
            }
            // write new transcript record to file
            records.push(new_transcript);
            
            // get a list of feature starts/stops, and a tree of CDS features
            let mut featurestarts = HashSet::<(String,u64)>::new();
            let mut featurestops = HashSet::<(String,u64)>::new();
            // for each child of the original transcript
            if let Some(ref transcript_row2children) = annot.row2children.get(transcript_row) {
                'child:
                for child_row in transcript_row2children.iter() {
                    // make a copy
                    let mut child = annot.rows[*child_row].clone();
                    // if this is a cassette that overlaps a reannotated pair,
                    // do not include it
                    for pair_row in &transcript2pair[transcript_row] {
                        let pair = &reannotated_pairs[*pair_row];
                        if child.start-1 < annot.rows[pair.exon2_row].start-1 &&
                            annot.rows[pair.exon1_row].end < child.end
                        { continue 'child; }
                    }
                    // populate featurestarts/stops
                    featurestarts.insert((child.feature_type.clone(), child.start-1));
                    featurestops.insert((child.feature_type.clone(), child.end));
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
                    
                    {   let entry = child.attributes.entry("Parent".to_string()).or_insert_with(|| "".to_string());
                        entry.clear();
                        entry.push_str(&mut newparents.join(","));
                    }
                    // update transcript_id and transcript_name attributes
                    if child.attributes.contains_key("transcript_id") {
                        if let Some(ref transcript_id) = transcript_id {
                            let entry = child.attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_id.clone());
                        }
                    }
                    if child.attributes.contains_key("transcript_name") {
                        if let Some(ref transcript_name) = transcript_name {
                            let entry = child.attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_name.clone());
                        }
                    }
                    if let Some(id) = child.attributes.get_mut("ID") {
                        let mut feature_id = format!("{}.reannot", id);
                        while ids.contains(&feature_id) {
                            feature_id = String::from(EXON_RENAME.replace(&feature_id.as_ref(), |caps: &Captures| {
                                format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                        }
                        id.clear();
                        id.push_str(&mut feature_id);
                        ids.insert(id.clone());
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
                        'cc:
                        for c in children {
                            if !seen.contains(&c)  {
                                seen.insert(c);
                                let mut cc = annot.rows[c].clone();
                                // if this is a cassette that overlaps a reannotated pair,
                                // do not include it
                                for pair_row in &transcript2pair[transcript_row] {
                                    let pair = &reannotated_pairs[*pair_row];
                                    if cc.start-1 < annot.rows[pair.exon2_row].start-1 &&
                                        annot.rows[pair.exon1_row].end < cc.end
                                    { continue 'cc; }
                                }
                                featurestarts.insert((cc.feature_type.clone(), cc.start-1));
                                featurestops.insert((cc.feature_type.clone(), cc.end));
                                // update transcript_id and transcript_name attributes
                                if cc.attributes.contains_key("transcript_id") {
                                    if let Some(ref transcript_id) = transcript_id {
                                        let entry = cc.attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                                        entry.clear();
                                        entry.push_str(&mut transcript_id.clone());
                                    }
                                }
                                if cc.attributes.contains_key("transcript_name") {
                                    if let Some(ref transcript_name) = transcript_name {
                                        let entry = cc.attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                                        entry.clear();
                                        entry.push_str(&mut transcript_name.clone());
                                    }
                                }
                                if let Some(id) = cc.attributes.get_mut("ID") {
                                    let mut feature_id = format!("{}.reannot", id);
                                    while ids.contains(&feature_id) {
                                        feature_id = String::from(EXON_RENAME.replace(&feature_id.as_ref(), |caps: &Captures| {
                                            format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                                    }
                                    id.clear();
                                    id.push_str(&mut feature_id);
                                    ids.insert(id.clone());
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
                            {   let entry = attributes.entry("ID".to_string()).or_insert_with(|| "".to_string());
                                entry.clear();
                                entry.push_str(&mut cassette_id.clone());
                            }
                            {   let entry = attributes.entry("exon_id".to_string()).or_insert_with(|| "".to_string());
                                entry.clear();
                                entry.push_str(&mut cassette_id.clone());
                            }
                        }
                        if let Some(ref transcript_id) = transcript_id {
                            let entry = attributes.entry("Parent".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_id.clone());
                        }
                        if let Some(ref transcript_id) = transcript_id {
                            let entry = attributes.entry("transcript_id".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_id.clone());
                        }
                        if let Some(ref transcript_name) = transcript_name {
                            let entry = attributes.entry("transcript_name".to_string()).or_insert_with(|| "".to_string());
                            entry.clear();
                            entry.push_str(&mut transcript_name.clone());
                        }
                        
                        let mut attributes = attributes.clone();
                        attributes.insert("exon_type".to_string(),"cassette".to_string());
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
                                let entry = attributes.entry("ID".to_string()).or_insert_with(|| "".to_string());
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
            if transcript.strand=="-" {
                records.sort_by(|a,b| b.start.cmp(&a.start));
            } else {
                records.sort_by(|a,b| a.start.cmp(&b.start));
            }
            let mut exontree = IntervalTree::<u64,Record>::new();
            let mut prevcdslen=0;
            let mut exon_number=1;
            for record in &mut records {
                let record = record;
                // compute the frame for CDS features
                if record.feature_type == "CDS" {
                    record.frame = (prevcdslen % 3).to_string();
                    prevcdslen += record.end-record.start+1;
                }
                // store exons in an interval tree to compute exon_number rank
                else if record.feature_type == "exon" {
                    {   let en = record.attributes.entry("exon_number".to_string()).or_insert_with(|| "".to_string());
                        en.clear();
                        en.push_str(&mut exon_number.to_string());
                    }
                    exon_number += 1;
                    exontree.insert(Interval::new((record.start - 1)..(record.end))?, record.clone());
                }
            }
            // compute exon_number rank
            for record in &mut records {
                if record.feature_type != "exon" {
                    for exon in exontree.find(record.start-1..record.end) {
                        if let Some(exon_number) = exon.data().attributes.get("exon_number") {
                            if let Some(en) = record.attributes.get_mut("exon_number") {
                                en.clear();
                                en.push_str(&mut exon_number.clone());
                            }
                        }
                    }
                }
            }
            if transcript.strand=="-" {
                records.sort_by(|a,b| a.start.cmp(&b.start));
            }
            for record in &mut records {
                writeln!(output, "{}", record.to_gff()?)?;
            }
        }
    }
    
    if let Some(ref debug_outannot_bigbed) = options.debug_outannot_bigbed {
        let newannot = IndexedAnnotation::from_gff(&outannot, 
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
    if let Some(ref debug_outannot_gtf) = options.debug_outannot_gtf {
        let newannot = IndexedAnnotation::from_gff(&outannot, 
            &options.chrmap_file,
            &options.vizchrmap_file)?;
        newannot.to_gtf(&debug_outannot_gtf)?; 
    }
    if let Some(ref debug_outannot_fasta) = options.debug_outannot_fasta {
        if let Some(ref genome_file) = options.genome_file {
            let newannot = IndexedAnnotation::from_gff(&outannot, 
                &options.chrmap_file,
                &options.vizchrmap_file)?;
            newannot.to_fasta(
                &debug_outannot_fasta, 
                &genome_file,
                &options.exon_type, 
                &options.transcript_type, 
                &options.gene_type)?; 
        }
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
        if options.debug_exon_bigbed.is_none() { options.debug_exon_bigbed = Some(format!("{}.exon.bb", options.debug_prefix)) }
        if options.debug_bigwig.is_none() { options.debug_bigwig = Some(format!("{}.intronrpkm", options.debug_prefix)) }
        if options.debug_reannot_bigbed.is_none() { options.debug_reannot_bigbed = Some(format!("{}.reannot.bb", options.debug_prefix)) }
        if options.debug_rpkm_region_bigbed.is_none() { options.debug_rpkm_region_bigbed = Some(format!("{}.rpkm.bb", options.debug_prefix)) }
        if options.debug_rpkmstats_json.is_none() { options.debug_rpkmstats_json = Some(format!("{}.rpkmstats.json", options.debug_prefix)) }
        if options.debug_trackdb.is_none() { options.debug_trackdb = Some(format!("{}.trackDb.txt", options.debug_prefix)) }
        if options.debug_outannot_bigbed.is_none() { options.debug_outannot_bigbed = Some(format!("{}.outannot.bb", options.debug_prefix)) }
        if options.debug_outannot_gtf.is_none() { options.debug_outannot_gtf = Some(format!("{}.outannot.gtf", options.debug_prefix)) }
        if options.debug_outannot_fasta.is_none() { options.debug_outannot_fasta = Some(format!("{}.outannot.fa", options.debug_prefix)) }
        if options.debug_retained_introns.is_none() { options.debug_retained_introns = Some(format!("{}.retained_introns.txt", options.debug_prefix)) }
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
        eprintln!("\nNo bam files were passed in!");
        std::process::exit(1);
    }
    if options.debug_outannot_fasta.is_some() && options.genome_file.is_none() {
        eprintln!("\nNo genome file was specified!");
        std::process::exit(1);
    }
    let transcript_type = String::from("transcript");
    let mut annot = if let Some(annotfile_gff) = options.annotfile_gff.clone() {
        eprintln!("Reading annotation file {:?}", &annotfile_gff);
        IndexedAnnotation::from_gff(
            &annotfile_gff, 
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else if let Some(annotfile_gtf) = options.annotfile_gtf.clone() {
        eprintln!("Reading annotation file {:?}", &annotfile_gtf);
        IndexedAnnotation::from_gtf(&annotfile_gtf, 
            options.gene_type.get(0).r()?, 
            options.transcript_type.get(0).unwrap_or(&transcript_type),
            &options.chrmap_file,
            &options.vizchrmap_file)?
    } else {
        Options::clap().print_help()?;
        bail!("No annotation file was given!");
    };
    eprintln!("Getting refseq lengths from bam file {:?}", &bamfiles[0]);
    let refs = match options.sizes_file.clone() {
        Some(sizes_file) => read_sizes_file(&sizes_file, &annot.chrmap)?,
        None => get_bam_refs(&bamfiles[0], &annot.chrmap)?,
    };
    annot.refs = refs;
    
    if let Some(ref debug_annot_gff) = options.debug_annot_gff {
        eprintln!("Writing annotation file to {:?}", &debug_annot_gff);
        annot.to_gff(&debug_annot_gff)?; 
    }
    if let Some(ref debug_annot_gtf) = options.debug_annot_gtf {
        eprintln!("Writing annotation file to {:?}", &debug_annot_gtf);
        annot.to_gtf(&debug_annot_gtf)?; 
    }
    if let Some(ref debug_annot_bigbed) = options.debug_annot_bigbed {
        eprintln!("Writing annotation file to {:?}", &debug_annot_bigbed);
        annot.to_bigbed(
            &debug_annot_bigbed, 
            &options.exon_type, 
            &options.cds_type, 
            &options.transcript_type,
            &options.gene_type,
            &mut trackdb)?; 
    }
    
    // get the total bam reads
    eprintln!("Running samtools idxstats to get total bam read counts");
    let total_reads = get_bam_total_reads(&bamfiles)?;
    eprintln!("Found {} total reads", total_reads);
    
    // find the constituitive exons
    eprintln!("Searching for constituitive exons in the annotation");
    let exonpairs = find_constituitive_splice_pairs(&annot, &options)?;
    
    let annot = Arc::new(annot);
    if let Some(ref debug_exon_bigbed) = options.debug_exon_bigbed {
        eprintln!("Writing constituitive pairs to a bigbed file");
        write_exon_bigbed(&exonpairs, &annot, &debug_exon_bigbed, &mut trackdb)?;
    }
        
    eprintln!("Reannotating cassette regions");
    let (reannotated_pairs, mut rpkmstats) = reannotate_regions(
        &annot,
        &exonpairs, 
        &bamfiles, 
        &bamstrand,
        total_reads,
        &options,
        &mut trackdb)?;
    if let Some(ref debug_reannot_bigbed) = options.debug_reannot_bigbed {
        eprintln!("Writing reannotation to bigbed");
        write_exon_bigbed(&reannotated_pairs, &annot, &debug_reannot_bigbed, &mut trackdb)?;
    }
    if let Some(ref debug_retained_introns) = options.debug_retained_introns {
        eprintln!("Writing retained introns to file");
        write_retained_introns(&reannotated_pairs, &annot, &debug_retained_introns)?;
    }
    
    eprintln!("Computing RPKM stats");
    if let Some(ref debug_rpkmstats_json) = options.debug_rpkmstats_json {
        eprintln!("Writing RPKM stats to {:?}", &debug_rpkmstats_json);
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if debug_rpkmstats_json == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(debug_rpkmstats_json)?) });
            serde_json::to_writer_pretty(&mut output, &rpkmstats)?;
    }
    eprintln!("Writing RPKM stats to {:?}", &options.outfile);
    write_rpkm_stats(&options.outfile, &mut rpkmstats)?;
    
    if let Some(ref outannot) = options.outannot {
        eprintln!("Writing output annotation file");
        write_enriched_annotation(&annot, &reannotated_pairs, &outannot, &options, &mut trackdb)?;
    }
    Ok(())
}

fn main() {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "full");

    if let Err(ref e) = run() {
        eprintln!("error: {}", e);
        eprintln!("backtrace: {:?}", e.backtrace());
        ::std::process::exit(1);
    }
}
