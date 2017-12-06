#![recursion_limit = "1024"]
use std::vec::Vec;
use std::ops::Range;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead};

extern crate linked_hash_map;
use linked_hash_map::LinkedHashMap;

#[macro_use] 
extern crate lazy_static;

extern crate regex;

extern crate rust_htslib;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read;
use rust_htslib::bam::IndexedReader;

extern crate bio;

extern crate csv;

#[macro_use]
extern crate failure;

extern crate structopt;

#[macro_use]
extern crate duct;

extern crate serde;
extern crate serde_derive;
extern crate serde_json;

#[macro_use]
extern crate url;

extern crate itertools;

extern crate unindent;

pub mod indexed_annotation;

pub mod error {
    pub type Result<T> = ::std::result::Result<T, ::failure::Error>;

    #[derive(Debug, Fail)]
    enum NoneError {
        #[fail(display = "Option value is None")]
        NoneError {}
    }

    pub trait ToResult<T> {
        fn r(self) -> Result<T>;
    }
    impl<T> ToResult<T> for Option<T> {
        fn r(self) -> Result<T> {
            match self {
                Some(v) => Ok(v),
                None => Err(NoneError::NoneError{}.into()),
            }
        }
    }
}
use ::error::*;

pub mod power_set {
    pub struct PowerSet<'a, T: 'a> {
        source: &'a [T],
        position: usize
    }
    
    impl<'a, T> PowerSet<'a, T> where T: Clone {
        pub fn new(source: &'a [T]) -> PowerSet<'a, T> {
            PowerSet { source: source, position: 0 }
        }
    }

    impl<'a, T> Iterator for PowerSet<'a, T> where T: Clone {
        type Item = Vec<T>;

        fn next(&mut self) -> Option<Self::Item> {
            if 2usize.pow(self.source.len() as u32) <= self.position {
                None
            } else {
                let res = self.source.iter().enumerate().filter(|&(i, _)| (self.position >> i) % 2 == 1)
                                                        .map(|(_, element)| element.clone()).collect();
                self.position = self.position + 1;
                Some(res)
            }
        }
    }
}

pub fn cigar2exons(cigar: &[Cigar], pos: u64) -> Result<Vec<Range<u64>>> {
    let mut exons = Vec::<Range<u64>>::new();
    let mut pos = pos;
    for op in cigar {
        match op {
            &Cigar::Match(length) => {
                pos += length as u64;
                if length > 0 {
                    exons.push(Range{start: pos - length as u64, end: pos});
                }
            }
            &Cigar::RefSkip(length) |
            &Cigar::Del(length) |
            &Cigar::Equal(length) |
            &Cigar::Diff(length) => {
                pos += length as u64;
            }
            &Cigar::Back(length) => {
                pos -= length as u64;
            }
            &Cigar::Ins(_) |
            &Cigar::SoftClip(_) |
            &Cigar::HardClip(_) |
            &Cigar::Pad(_) => (),
        };
    }
    Ok(exons)
}

pub fn read_sizes_file(sizes_file: &str, chrmap: &HashMap<String,String>) -> Result<LinkedHashMap<String,u64>> {
    let mut refs = HashMap::<String,u64>::new();
    let f = File::open(&sizes_file)?;
    let mut file = BufReader::new(&f);
    let mut buf = String::new();
    while file.read_line(&mut buf)? > 0 {
        {   let line = buf.trim_right_matches('\n').trim_right_matches('\r');
            let cols: Vec<&str> = line.split('\t').collect();
            if let Some(chr) = cols.get(0) {
                let chr = String::from(*chr);
                let chr = chrmap.get(&chr).unwrap_or(&chr);
                if let Some(size) = cols.get(1) {
                    if let Ok(size) = size.parse::<u64>() {
                        refs.insert(chr.clone(), size);
                    }
                    else {
                        bail!("Could not parse size \"{}\" for chr \"{}\" from line \"{}\" of file \"{}\"", size, chr, line, sizes_file);
                    }
                }
            }
        }
        buf.clear();
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

pub fn get_bam_refs(bamfile: &str, chrmap: &HashMap<String,String>) -> Result<LinkedHashMap<String,u64>> {
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

pub fn get_bam_total_reads(bamfiles: &[String]) -> Result<u64> {
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
