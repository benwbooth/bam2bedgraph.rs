#![recursion_limit="128"]
#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::str;
use std::vec::Vec;
use std::collections::HashMap;
use std::sync::Arc;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::cmp::{min, max};
use std::ops::Range;

#[macro_use] 
extern crate failure;

extern crate bam2bedgraph;
use bam2bedgraph::*;
use bam2bedgraph::error::*;
use bam2bedgraph::indexed_annotation::*;

extern crate rust_htslib;
use rust_htslib::bam::Read;
use rust_htslib::bam::IndexedReader;

extern crate structopt;
#[macro_use]
extern crate structopt_derive;
use structopt::StructOpt;

extern crate futures;
use futures::Future;
extern crate futures_cpupool;
use futures_cpupool::CpuPool;

extern crate num_cpus;

extern crate serde;
#[macro_use]
extern crate serde_derive;

extern crate csv;

#[derive(StructOpt, Debug)]
#[structopt(name = "adjusted_intron_cov", about = "Compute adjusted intron coverage and PSI values from spladder output")]
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
    #[structopt(long="input", short="i", help = "Spladder input file", name="INPUT", default_value="-")]
    input: String,
    // feature types filter
    #[structopt(long="exon_type", help = "The exon type(s) to search for", name="EXON_TYPE")]
    exon_type: Vec<String>,
    #[structopt(long="transcript_type", help = "The transcript type(s) to search for", name="TRANSCRIPT_TYPE")]
    transcript_type: Vec<String>,
    #[structopt(long="gene_type", help = "The gene type(s) to search for", name="GENE_TYPE")]
    gene_type: Vec<String>,
    // flags
    #[structopt(long="cpu_threads", short="t", help = "How many threads to use for processing", default_value="0")]
    cpu_threads: usize,
}

#[derive(Debug,Serialize,Deserialize,Clone)]
struct Row {
    contig: String,
    strand: String,
    event_id: String,
    gene_name: String,
    exon1_start: u64,
    exon1_end: u64,
    intron_start: u64,
    intron_end: u64,
    exon2_start: u64,
    exon2_end: u64,
    exon1_cov: f64,
    intron_cov: f64,
    exon2_cov: f64,
    intron_conf: u64,
    psi: String,
}

#[derive(Debug,Serialize,Deserialize,Clone)]
struct OutRow {
    contig: String,
    strand: String,
    event_id: String,
    gene_name: String,
    exon1_start: u64,
    exon1_end: u64,
    intron_start: u64,
    intron_end: u64,
    exon2_start: u64,
    exon2_end: u64,
    exon1_cov: f64,
    intron_cov: f64,
    exon2_cov: f64,
    intron_conf: u64,
    psi: String,
    intron_bases: u64,
    intron_coverage: u64,
    exon_bases: u64,
    exon_coverage: u64,
    total_bases: u64,
    total_coverage: u64,
    adjusted_total_intron_coverage: u64,
    adjusted_intron_coverage: String,
    adjusted_psi: String,
}

fn write_intron_cov(
    options: &Options,
    bamfiles: &Vec<String>,
    bamstrand: &Vec<Option<bool>>,
    annot: &Arc<IndexedAnnotation>)
    -> Result<()> 
{
    let input: Box<std::io::Read> = match options.input.as_ref() {
        "-" => Box::new(std::io::stdin()),
        _ => Box::new(std::fs::File::open(&options.input)? )};

    let output: BufWriter<Box<Write>> = BufWriter::new(
        if options.outfile == "-" { Box::new(std::io::stdout()) }
            else { Box::new(File::create(&options.outfile)?) });

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

    let num_cpus = num_cpus::get();
    let mut pair_futures = Vec::new();
    let pool = Arc::new(CpuPool::new(if options.cpu_threads==0 {num_cpus} else {options.cpu_threads}));
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(input);
    for result in rdr.deserialize().skip(1) {
        let row: Row = result?;
        let row = row.clone();
        let annot = annot.clone();
        let bamfiles = bamfiles.clone();
        let bamstrand = bamstrand.clone();
        let tidmaps = tidmaps.clone();
        let pair_future = pool.spawn_fn(move || -> Result<OutRow> {
            //get all the bam reads in parallel
            let mut coverage = vec![0u64; (row.intron_end-row.intron_start+1) as usize];
            for (i, bamfile) in bamfiles.iter().enumerate() {
                let read1strand = bamstrand[i];
                let tidmap = &tidmaps[bamfile];
                if let Some(tid) = tidmap.get(&row.contig) {
                    let mut bam = IndexedReader::from_path(bamfile)?;
                    bam.fetch(*tid, (row.intron_start - 1) as u32, row.intron_end as u32)?;
                    for read in bam.records() {
                        let read = read?;
                        // make sure the read's strand matches
                        if let Some(is_read1strand) = read1strand {
                            if !(((is_read1strand == read.is_first_in_template()) == !read.is_reverse()) == (row.strand == "+"))
                            {
                                continue;
                            }
                        }
                        let exons = cigar2exons(&read.cigar(), read.pos() as u64)?;
                        for exon in exons {
                            if exon.start < row.intron_end && row.intron_start - 1 < exon.end {
                                for i in max(exon.start, row.intron_start-1)..min(exon.end, row.intron_end) {
                                    coverage[(i-(row.intron_start-1)) as usize] += 1
                                }
                            }
                        }
                    }
                }
            }
            // get the annotated exons
            let mut annotated_exons = Vec::new();
            let chr = annot.chrmap.get(&row.contig).unwrap_or(&row.contig);
            if let Some(tree) = annot.tree.get(chr) {
                annotated_exons = tree.find(row.intron_start - 1..row.intron_end).
                    map(|t| &annot.rows[*t.data()]).
                    filter(|t| t.feature_type == "exon" && t.strand == row.strand).
                    collect::<Vec<_>>();
            }
            annotated_exons.sort_by_key(|e| e.start);

            // merge the annotated exons
            let mut merged_exons = Vec::<Range<u64>>::new();
            for exon in annotated_exons {
                if let Some(ref mut last) = merged_exons.last_mut() {
                    if (exon.start-1) <= last.end && (last.start-1) <= exon.end {
                        last.end = exon.end;
                        continue;
                    }
                }
                merged_exons.push((exon.start-1)..exon.end);
            }

            let mut intron_bases = 0u64;
            let mut intron_coverage = 0u64;
            let mut exon_bases = 0u64;
            let mut exon_coverage = 0u64;
            for exon in &merged_exons {
                if (row.intron_start-1) < exon.start {
                    intron_bases += exon.start-(row.intron_start-1);
                    for j in (row.intron_start-1)..exon.start {
                        intron_coverage += coverage[(j-(row.intron_start-1)) as usize]
                    }
                }
                exon_bases += max(row.intron_end, exon.end)-min(row.intron_start-1, exon.start);
                for j in max(row.intron_start-1, exon.start)..min(row.intron_end, exon.end) {
                    exon_coverage += coverage[(j-(row.intron_start-1)) as usize];
                }
            }
            if let Some(last) = merged_exons.last() {
                if last.end < row.intron_end {
                    intron_bases += row.intron_end-last.end;
                    for j in last.end..row.intron_end {
                        intron_coverage += coverage[(j-(row.intron_start-1)) as usize];
                    }
                }
            }
            else {
                intron_bases += row.intron_end-(row.intron_start-1);
                for j in (row.intron_start-1)..row.intron_end {
                    intron_coverage += coverage[(j-(row.intron_start-1)) as usize];
                }
            }
            let total_coverage = intron_coverage + exon_coverage;
            let total_bases = intron_bases + exon_bases;
            let adjusted_total_intron_coverage = intron_coverage + ((intron_coverage as f64 / intron_bases as f64) * exon_bases as f64).floor() as u64;
            let adjusted_intron_coverage = adjusted_total_intron_coverage as f64 / total_bases as f64;
            let adjusted_psi = adjusted_intron_coverage / (adjusted_intron_coverage + row.intron_conf as f64);
            Ok(OutRow{
                contig: row.contig.clone(),
                strand: row.strand.clone(),
                event_id: row.event_id.clone(),
                gene_name: row.gene_name.clone(),
                exon1_start: row.exon1_start,
                exon1_end: row.exon1_end,
                intron_start: row.intron_start,
                intron_end: row.intron_end,
                exon2_start: row.exon2_start,
                exon2_end: row.exon2_end,
                exon1_cov: row.exon1_cov,
                intron_cov: row.intron_cov,
                exon2_cov: row.exon2_cov,
                intron_conf: row.intron_conf,
                psi: row.psi.clone(),
                intron_bases,
                intron_coverage,
                exon_bases,
                exon_coverage,
                total_bases,
                total_coverage,
                adjusted_total_intron_coverage,
                adjusted_intron_coverage: format!("{:.*}", 2, adjusted_intron_coverage),
                adjusted_psi: format!("{:.*}", 2, adjusted_psi),
            })
        });
        pair_futures.push(pair_future);
    }
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .quote_style(csv::QuoteStyle::Necessary)
        .from_writer(output);
    for future in pair_futures {
        match future.wait() {
            Ok(outrow) => {
                wtr.serialize(outrow)?;
            }
            Err(ref e) => {
                eprintln!("Got Err in write_exon_cov: {:?}", e);
            }
        }
    }
    wtr.flush()?;
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
    // set debug options if --debug flag is set

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
        eprintln!("\n\nNo bam files were passed in!");
        std::process::exit(1)
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
        eprintln!("\n\nNo annotation file was given!");
        std::process::exit(1);
    };
    eprintln!("Getting refseq lengths from bam file {:?}", &bamfiles[0]);
    let refs = match options.sizes_file.clone() {
        Some(sizes_file) => read_sizes_file(&sizes_file, &annot.chrmap)?,
        None => get_bam_refs(&bamfiles[0], &annot.chrmap)?,
    };
    annot.refs = refs;
    
    // get the total bam reads
    eprintln!("Running samtools idxstats to get total bam read counts");
    let total_reads = get_bam_total_reads(&bamfiles)?;
    eprintln!("Found {} total reads", total_reads);
    write_intron_cov(&options, &bamfiles, &bamstrand,  &Arc::new(annot))?;
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
