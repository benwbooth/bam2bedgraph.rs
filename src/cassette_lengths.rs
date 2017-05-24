#![recursion_limit="128"]
#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::io::stderr;
use std::io::Write;
use std::collections::HashMap;

extern crate bam2bedgraph;
use bam2bedgraph::errors::*;
use bam2bedgraph::indexed_annotation::IndexedAnnotation;

extern crate bio;

extern crate structopt;
#[macro_use]
extern crate structopt_derive;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "cassette_lengths", about = "Histogram of cassette exon lengths within the intron of tr of another gene")]
struct Options {
    // input files
    #[structopt(long="overlapping_genes", help = "Only show cassettes which overlap other genes?")]
    overlapping_genes: bool,
    // input files
    #[structopt(long="gff", help = "A genome annotation file in gff3 format", name="ANNOT_GFF_FILE")]
    annotfile_gff: Option<String>,
}

fn run() -> Result<()> {
    let options = Options::from_args();
    let annot = IndexedAnnotation::from_gff(
        &options.annotfile_gff.r()?, 
        &None,
        &None)?;
    println!("row\tgene_name\tother_gene_name\tcassette_location\tcassette_length");
    for (row, record) in annot.rows.iter().enumerate() {
        if let Some(exon_type) = record.attributes.get("exon_type") {
            if exon_type == "cassette" && record.feature_type == "exon" {
                if let Some(record_gene_name) = record.attributes.get("gene_name") {
                    let mut other_genes = HashMap::<String,u64>::new();
                    for overlap_row in annot.tree[&record.seqname].find(record.start-1..record.end) {
                        let overlaps = &annot.rows[*overlap_row.data()];
                        if overlaps.strand == record.strand && 
                            (overlaps.feature_type == "gene" || overlaps.feature_type == "exon")
                        {
                            if let Some(overlaps_gene_name) = overlaps.attributes.get("gene_name") {
                                if overlaps_gene_name != record_gene_name && 
                                    &format!("{}.reannot", overlaps_gene_name) != record_gene_name 
                                {
                                    if overlaps.feature_type == "exon" {
                                        *other_genes.entry(overlaps_gene_name.clone()).or_insert(0) += 1;
                                    }
                                    else {
                                        other_genes.entry(overlaps_gene_name.clone()).or_insert(0);
                                    }
                                }
                            }
                        }
                    }
                    let mut overlapping_genes = Vec::<String>::new();
                    for (gene_name,overlapping_exons) in other_genes {
                        if overlapping_exons == 0 {
                            overlapping_genes.push(gene_name.clone());
                        }
                    }
                    if !options.overlapping_genes || !overlapping_genes.is_empty() {
                        println!("{}\t{}\t{}\t{}:{}..{}:{}\t{}", 
                            row,
                            record_gene_name,
                            overlapping_genes.join(","),
                            record.seqname,
                            record.start,
                            record.end,
                            record.strand,
                            record.end-record.start+1);
                    }
                }
            }
        }
    }
    
    Ok(())
}

fn main() {
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
