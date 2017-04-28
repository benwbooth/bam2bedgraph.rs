#![recursion_limit = "1024"]
use std::vec::Vec;
use std::ops::Range;

#[macro_use] 
extern crate lazy_static;

extern crate regex;

extern crate rust_htslib;
use rust_htslib::bam::record::Cigar;

extern crate bio;

extern crate csv;

#[macro_use]
extern crate error_chain;

extern crate structopt;

#[macro_use]
extern crate duct;

extern crate serde;
extern crate serde_derive;
extern crate serde_json;

extern crate linked_hash_map;

#[macro_use]
extern crate url;

extern crate itertools;

extern crate unindent;

pub mod indexed_annotation;

pub mod errors {
    error_chain!{
        foreign_links {
            Io(::std::io::Error) #[cfg(unix)];
            Utf8(::std::str::Utf8Error);
            FromUtf8(::std::string::FromUtf8Error);
            Regex(::regex::Error);
            ReaderPath(::rust_htslib::bam::ReaderPathError);
            Read(::rust_htslib::bam::ReadError);
            IndexedReaderPath(::rust_htslib::bam::IndexedReaderPathError);
            ParseInt(::std::num::ParseIntError);
            Seek(::rust_htslib::bam::SeekError);
            Interval(::bio::utils::IntervalError);
            Csv(::csv::Error);
            Json(::serde_json::Error);
            StructOpt(::structopt::clap::Error);
        }
        errors {
            NoneError
        }
    }
    pub trait ToResult<T> {
        fn r(self) -> Result<T>;
    }
    impl<T> ToResult<T> for Option<T> {
        fn r(self) -> Result<T> {
            match self {
                Some(v) => Ok(v),
                None => Err(ErrorKind::NoneError.into()),
            }
        }
    }
}

use errors::*;

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
