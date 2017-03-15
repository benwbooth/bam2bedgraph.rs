#![recursion_limit = "1024"]
use std::str;
use std::vec::Vec;

extern crate regex;

extern crate rust_htslib;
use rust_htslib::bam::record::Cigar;

extern crate bio;

extern crate csv;

#[macro_use]
extern crate error_chain;

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

pub fn cigar2exons(exons: &mut Vec<(u64, u64)>, cigar: &[Cigar], pos: u64) -> Result<()> {
    let mut pos = pos;
    for op in cigar {
        match op {
            &Cigar::Match(length) => {
                pos += length as u64;
                exons.push((pos - length as u64, pos));
                Ok(())
            }
            &Cigar::RefSkip(length) |
            &Cigar::Del(length) => {
                pos += length as u64;
                Ok(())
            }
            &Cigar::Ins(_) |
            &Cigar::SoftClip(_) |
            &Cigar::HardClip(_) |
            &Cigar::Pad(_) => Ok(()),
            c => Err(format!("Bad CIGAR string: {:?}", c)),
        }?;
    }
    Ok(())
}
