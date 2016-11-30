#![recursion_limit = "1024"]
use std::ascii::AsciiExt;
use std::collections::HashMap;
use std::io::Write;
use std::io::BufWriter;
use std::io::stderr;
use std::fs::File;
use std::str;
use std::path::{PathBuf, Path};
use std::process::Command;
use std::vec::Vec;
use std::marker::Send;

extern crate argparse;
use argparse::{ArgumentParser, StoreTrue, StoreFalse, Store};

extern crate regex;
use regex::Regex;

extern crate rust_htslib;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Reader;

#[macro_use]
extern crate error_chain;

mod intervaltree;

error_chain!{
    foreign_links {
        ::std::io::Error, Io;
        ::std::str::Utf8Error, Utf8;
        ::regex::Error, Regex;
    }
}

fn cigar2exons(exons: &mut Vec<(i32, i32)>, cigar: &[Cigar], pos: i32) -> Result<()> {
    let mut pos = pos;
    for op in cigar {
        match op {
            &Cigar::Match(length) => {
                pos += length as i32;
                exons.push((pos - length as i32, pos));
                Ok(())
            }
            &Cigar::RefSkip(length) |
            &Cigar::Del(length) => {
                pos += length as i32;
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

fn open_file(options: &Options,
             read_number: i32,
             strand: &str,
             split_strand: &str,
             fhs: &mut HashMap<String, Option<File>>)
             -> Result<String> {
    let mut prefix = PathBuf::new();
    prefix.set_file_name(&options.bamfile);
    prefix.set_extension("");
    let track_name = vec![if !options.trackname.is_empty() {
                              options.trackname.clone()
                          } else {
                              let p = prefix.as_path().to_str().ok_or("Failed to parse path!")?;
                              p.to_string()
                          },
                          if options.split_read && read_number > 0 {
                              format!(".r{}", read_number)
                          } else {
                              "".to_string()
                          },
                          if split_strand != "uu" && !strand.is_empty() {
                              format!(".{}", strand)
                          } else {
                              "".to_string()
                          }]
        .join("");

    let filename = vec![if !options.out.is_empty() {
                            options.out.clone()
                        } else {
                            let p = prefix.as_path().to_str().ok_or("Failed to parse path!")?;
                            p.to_string()
                        },
                        if options.split_read && read_number > 0 {
                            format!(".r{}", read_number)
                        } else {
                            "".to_string()
                        },
                        if split_strand != "uu" && !strand.is_empty() {
                            format!(".{}", strand)
                        } else {
                            "".to_string()
                        },
                        ".bedgraph".to_string()]
        .join("");

    // initialize the file if needed
    if !fhs.contains_key(&filename) {
        let mut f = File::create(&filename)?;
        if options.trackline {
            writeln!(f,
                     "track type=bedGraph name=\"{}\" description=\"{}\" visibility=full",
                     track_name,
                     track_name)?;
        }
        fhs.insert(filename.clone(), Some(f));
    }
    Ok(filename)
}

fn write_chr(options: &Options,
             chr: &(u32, String),
             histogram: &HashMap<(i32, String), Vec<i32>>,
             fhs: &mut HashMap<String, Option<File>>,
             split_strand: &str)
             -> Result<()> {
    for (key, histo) in histogram {
        let read_number = key.0;
        let strand = &key.1;
        let filename = open_file(options, read_number, strand, split_strand, fhs)?;
        let mut f = fhs.get_mut(&filename).ok_or("Error!")?;
        let file = f.as_mut().ok_or("Error!")?;
        let mut writer = BufWriter::new(file);

        // scan the histogram to produce the bedgraph data
        let mut start: usize = 0;
        let mut end: usize = 0;
        let ref_length: usize = chr.0 as usize;
        while start < ref_length {
            while (if end < histo.len() { histo[end] } else { 0 }) ==
                  (if start < histo.len() { histo[start] } else { 0 }) &&
                  end < ref_length {
                end += 1
            }
            if options.zero || (if start < histo.len() { histo[start] } else { 0 }) > 0 {
                writeln!(writer,
                         "{}\t{}\t{}\t{}",
                         chr.1,
                         start,
                         end,
                         if strand == "-" {
                             -histo[start]
                         } else {
                             histo[start]
                         })?;
            }
            start = end;
        }
    }
    Ok(())
}

fn analyze_bam(options: &Options,
               split_strand: &str,
               autostrand_pass: bool,
               intervals: &Option<HashMap<String, intervaltree::IntervalTree<u8>>>)
               -> Result<()> {
    if !Path::new(&options.bamfile).exists() {
        return Err(format!("Bam file {} could not be found!", &options.bamfile).into());
    }
    let bam = (match Reader::from_path(&options.bamfile) {
        Ok(reader) => Ok(reader),
        Err(_) => Err("Error: BGZFError"),
    })?;
    let header = bam.header();

    let mut refs: Vec<(u32, String)> = Vec::new();
    refs.resize(header.target_count() as usize, (0, "".to_string()));
    let target_names = header.target_names();
    for target_name in target_names {
        let tid = header.tid(target_name).ok_or("Error!")?;
        let target_len = header.target_len(tid).ok_or("Error!")?;
        let target_name = std::str::from_utf8(target_name)?;
        refs[tid as usize] = (target_len, target_name.to_string());
    }

    if options.fixchr {
        for r in &mut refs {
            let regex = Regex::new(r"^(chr|Zv9_)")?;
            if regex.is_match(&r.1) {
                let refname = r.1.to_string();
                r.1.clear();
                r.1.push_str(&format!("chr{}", refname));
            }
        }
    }
    if autostrand_pass {
        writeln!(stderr(),
                 "Running strand detection phase on {}",
                 options.bamfile)?;
    } else {
        writeln!(stderr(), "Building histograms for {}", options.bamfile)?;
    }

    // build a lookup map for the refseqs
    let mut refmap: HashMap<String, usize> = HashMap::new();
    for (i, _) in refs.iter().enumerate() {
        refmap.insert(refs[i].1.to_string(), i);
    }

    let mut lastchr: i32 = -1;
    let mut fhs: HashMap<String, Option<File>> = HashMap::new();
    let mut histogram: HashMap<(i32, String), Vec<i32>> = HashMap::new();

    let mut autostrand_totals: HashMap<char, i64> = HashMap::new();
    autostrand_totals.insert('s', 0);
    autostrand_totals.insert('r', 0);
    let mut autostrand_totals2: HashMap<char, i64> = HashMap::new();
    autostrand_totals2.insert('s', 0);
    autostrand_totals2.insert('r', 0);

    let mut read = rust_htslib::bam::record::Record::new();
    while bam.read(&mut read).is_ok() {
        // if we've hit a new chr, write out the bedgraph data and clear the histogram
        if lastchr == -1 || read.tid() != lastchr {
            if !autostrand_pass && !histogram.is_empty() && lastchr != -1 {
                write_chr(options,
                          &refs[lastchr as usize],
                          &histogram,
                          &mut fhs,
                          split_strand)?;
            }
            histogram.clear();
            lastchr = read.tid();
        }

        // skip this read if it's no good
        let paired = read.is_paired();
        let proper = read.is_proper_pair();
        let primary = !read.is_secondary();
        if (options.paired_only && !paired) || (options.primary_only && !primary) ||
           (options.proper_only && !proper) {
            continue;
        }

        // skip if it's not unique and we want unique alignments
        if options.uniq {
            let hits = read.aux("NH".to_string().as_bytes());
            if hits == None || hits.ok_or("Error!")?.integer() != 1 {
                continue;
            }
        }

        let mut exons: Vec<(i32, i32)> = Vec::new();
        let mut get_exons: Vec<(i32, i32)> = Vec::new();
        cigar2exons(&mut get_exons, &read.cigar(), read.pos())?;
        if !options.split_exons && !get_exons.is_empty() {
            let first = get_exons.get(0).ok_or("Error!")?;
            let last = get_exons.get(get_exons.len() - 1).ok_or("Error!")?;
            exons = vec![(first.0, last.1)];
        }
        let read_number = if read.is_secondary() { 2 } else { 1 };

        // attempt to determine the strandedness of the transcript
        // read numbers match, is not reverse, is not flipped
        // let xs = read.aux("XS".as_bytes());
        let strand = 
        //if xs.is_some() {
         //   str::from_utf8(xs.ok_or("Error!")?.string())?
        //} else 
        if read_number == 1 {
                if split_strand.chars().nth(0).ok_or("Error!")? == 'r' {
                    if read.is_reverse() { "+" } else { "-" }
                } else if split_strand.chars().nth(0).ok_or("Error!")? == 's' {
                    if read.is_reverse() { "-" } else { "+" }
                } else {
                    ""
                }
            } else if read_number == 2 {
                if split_strand.chars().nth(1).ok_or("Error!")? == 's' {
                    if read.is_reverse() { "-" } else { "+" }
                } else if split_strand.chars().nth(1).ok_or("Error!")? == 'r' {
                    if read.is_reverse() { "+" } else { "-" }
                } else {
                    ""
                }
            } else {
                ""
            };

        let read_num = if options.split_read { read_number } else { 0 };

        let ref_length = refs[read.tid() as usize].0;

        // add the read to the histogram
        for exon in exons {
            // try to determine the strandedness of the data
            if autostrand_pass {
                if intervals.is_some() {
                    let intervals = intervals.as_ref().ok_or("Error!")?;
                    if intervals.contains_key(&refs[lastchr as usize].1) {
                        let mut overlapping_annot: Vec<intervaltree::Interval<u8>> = Vec::new();
                        intervals[&refs[lastchr as usize].1]
                            .find_overlapping(exon.0 + 1, exon.1, &mut overlapping_annot);

                        for interval in overlapping_annot {
                            let overlap_length = std::cmp::min(exon.1, interval.stop) -
                                                 std::cmp::max(exon.0, interval.start - 1);

                            let strandtype = if read.is_reverse() == (interval.value == b'-') {
                                's'
                            } else {
                                'r'
                            };
                            if read_number == 1 {
                                let at = autostrand_totals.get_mut(&strandtype).ok_or("Error!")?;
                                *at += overlap_length as i64
                            } else if read_number == 2 {
                                let at2 = autostrand_totals2.get_mut(&strandtype).ok_or("Error!")?;
                                *at2 += overlap_length as i64
                            }
                        }
                    }
                }
            } else {
                let tuple = (read_num, strand.to_string());
                if !histogram.contains_key(&tuple) {
                    histogram.insert(tuple.clone(), Vec::new());
                }
                // keep track of chromosome sizes
                if ref_length < exon.1 as u32 {
                    refs[read.tid() as usize].0 = exon.1 as u32;
                }
                if histogram[&tuple].len() < ref_length as usize {
                    let h = histogram.get_mut(&tuple).ok_or("Error!")?;
                    h.resize(ref_length as usize, 0);
                }

                for pos in exon.0..exon.1 {
                    let h = histogram.get_mut(&tuple).ok_or("Error!")?;
                    (*h)[pos as usize] += 1;
                }
            }
        }
    }

    if !autostrand_pass && !histogram.is_empty() && lastchr != -1 {
        write_chr(options,
                  &refs[lastchr as usize],
                  &histogram,
                  &mut fhs,
                  split_strand)?;
    }

    // make sure empty files were created
    if histogram.is_empty() && !autostrand_pass {
        for read_number in if options.split_read {
            vec![1, 2]
        } else {
            vec![0]
        } {
            for s in if split_strand != "uu" {
                vec!["+", "-"]
            } else {
                vec![""]
            } {
                open_file(options, read_number, s, split_strand, &mut fhs)?;
            }
        }
    }

    // close the filehandles
    for (_, fh) in &mut fhs {
        *fh = None;
    }

    if autostrand_pass {
        // get the read 1 and read2 totals
        let mut total1: i64 = 0;
        let mut total2: i64 = 0;
        for i in &autostrand_totals {
            total1 += *i.1;
        }
        for i in &autostrand_totals2 {
            total2 += *i.1;
        }
        // figure out the best and second-best strand types for reads 1 and 2
        let mut best1: (char, i64) = ('\0', 0);
        let mut second_best1: (char, i64) = ('\0', 0);
        let mut best2: (char, i64) = ('\0', 0);
        let mut second_best2: (char, i64) = ('\0', 0);
        for i in autostrand_totals {
            if i.1 > 0 {
                writeln!(stderr(),
                         "Total evidence for read 1 strand type {}: {}",
                         i.0,
                         i.1)?;
            }
            if best1.1 < i.1 {
                second_best1 = best1;
                best1 = i;
            } else if second_best1.1 < i.1 {
                second_best1 = i;
            }
        }
        for i in autostrand_totals2 {
            if i.1 > 0 {
                writeln!(stderr(),
                         "Total evidence for read 2 strand type {}: {}",
                         i.0,
                         i.1)?;
            }
            if best2.1 < i.1 {
                second_best2 = best2;
                best2 = i;
            } else if second_best2.1 < i.1 {
                second_best2 = i;
            }
        }
        let threshold: f64 = 0.0; //threshold isn't working, set to zero
        let strand1 = if total1 > 0 {
            if threshold < (best1.1 - second_best1.1) as f64 / (total1) as f64 {
                best1.0
            } else {
                'u'
            }
        } else {
            'u'
        };
        let strand2 = if total2 > 0 {
            if threshold < (best2.1 - second_best2.1) as f64 / (total2) as f64 {
                best2.0
            } else {
                'u'
            }
        } else {
            'u'
        };

        let best_strand = format!("{}{}", strand1, strand2);
        writeln!(stderr(),
                 "autostrand_pass found best strand type: {}",
                 best_strand)?;

        // re-run analyzeBam with the strand type indicated
        analyze_bam(options, &best_strand, false, intervals)?;
    }

    if !autostrand_pass && options.bigwig {
        // Convert bedgraph file to bigwig file
        for fh in &fhs {
            // write the genome file for bigwigs
            let genome_filename = format!("{}.genome", fh.0);
            {
                let mut genome_fh = File::create(&genome_filename)?;
                for r in &refs {
                    writeln!(genome_fh, "{}\t{}", r.1, r.0)?;
                }
            }

            // run bedGraphToBigWig
            let regex = Regex::new(r"\.bedgraph$")?;
            let bigwig_file = regex.replace(fh.0, ".bw");
            let sorted_bedgraph = regex.replace(fh.0, ".sorted.bedgraph");
            for command in vec![Command::new("sort")
                                    .args(&["-k1,1", "-k2,2n", "-o", &sorted_bedgraph, fh.0])
                                    .env("LC_COLLATE", "C"),
                                Command::new("bedGraphToBigWig")
                                    .args(&[&sorted_bedgraph, &genome_filename, &bigwig_file])] {
                let mut child = command.spawn()?;
                let exit_code = child.wait()?;
                if !exit_code.success() {
                    if exit_code.code().is_some() {
                        let code = exit_code.code().ok_or("Error!")?;
                        return Err(format!("Nonzero exit code {} returned from \
                                                       command: {:?}",
                                           code,
                                           command)
                            .into());
                    } else {
                        return Err(format!("Command was interrupted: {:?}", command).into());
                    }
                }
            }
            // remove the bedgraph file
            std::fs::remove_file(fh.0)?;
            // remove the sorted bedgraph file
            std::fs::remove_file(sorted_bedgraph)?;
            // remove the genome file
            std::fs::remove_file(&genome_filename)?;
        }
    };
    Ok(())
}

struct Options {
    split_exons: bool,
    split_read: bool,
    zero: bool,
    paired_only: bool,
    proper_only: bool,
    primary_only: bool,
    trackline: bool,
    bigwig: bool,
    uniq: bool,
    fixchr: bool,
    bamfile: String,
    trackname: String,
    out: String,
    autostrand: String,
    split_strand: String,
}
impl Options {
    fn default() -> Options {
        Options {
            split_exons: false,
            split_read: false,
            zero: false,
            paired_only: false,
            proper_only: false,
            primary_only: false,
            trackline: false,
            bigwig: false,
            uniq: false,
            fixchr: false,
            bamfile: "".to_string(),
            trackname: "".to_string(),
            out: "".to_string(),
            autostrand: "".to_string(),
            split_strand: "uu".to_string(),
        }
    }
}

fn run() -> Result<()> {
    // enable stack traces
    std::env::set_var("RUST_BACKTRACE", "1");

    let mut options = Options { ..Options::default() };
    {
        let mut ap = ArgumentParser::new();
        ap.set_description("Convert a bam file into a bedgraph/bigwig file.");
        ap.refer(&mut options.bamfile)
            .add_argument("BAMFILE", Store, "Input BAM filename")
            .required();
        ap.refer(&mut options.split_exons)
            .add_option(&["--split"],
                        StoreTrue,
                        "Use CIGAR string to split alignment into separate exons (default)")
            .add_option(&["--nosplit"], StoreFalse, "");
        ap.refer(&mut options.autostrand)
            .add_option(&["--autostrand"],
                        Store,
                        "Attempt to determine the strandedness of the input data using an \
                         annotation file. Must be a .bam file.")
            .metavar("ANNOT_BAMFILE");
        ap.refer(&mut options.split_strand)
            .add_option(&["--strand"],
                        Store,
                        "Split output bedgraph by strand: Possible values: u s r uu us ur su ss \
                         sr ru rs rr, first char is read1, second is read2, u=unstranded, \
                         s=stranded, r=reverse")
            .metavar("[TYPE]");
        ap.refer(&mut options.split_read)
            .add_option(&["--read"],
                        StoreTrue,
                        "Split output bedgraph by read number")
            .add_option(&["--noread"], StoreFalse, "(default)");
        ap.refer(&mut options.zero)
            .add_option(&["--zero"], StoreTrue, "Pad output bedgraph with zeroes")
            .add_option(&["--nozero"], StoreFalse, "(default)");
        ap.refer(&mut options.fixchr)
            .add_option(&["--fixchr"],
                        StoreTrue,
                        "Transform chromosome names to be UCSC-compatible")
            .add_option(&["--nofixchr"], StoreFalse, "(default)");
        ap.refer(&mut options.paired_only)
            .add_option(&["--paired"],
                        StoreTrue,
                        "Only output paired read alignments")
            .add_option(&["--nopaired"], StoreFalse, "(default)");
        ap.refer(&mut options.proper_only)
            .add_option(&["--proper"],
                        StoreTrue,
                        "Only output proper-paired read alignments")
            .add_option(&["--noproper"], StoreFalse, "(default)");
        ap.refer(&mut options.primary_only)
            .add_option(&["--primary"], StoreTrue, "Only output primary alignments")
            .add_option(&["--noprimary"], StoreFalse, "(default)");
        ap.refer(&mut options.bigwig)
            .add_option(&["--bigwig"],
                        StoreTrue,
                        "Output bigwig files (requires bedGraphToBigWig in $PATH)")
            .add_option(&["--nobigwig"], StoreFalse, "(default)");
        ap.refer(&mut options.uniq)
            .add_option(&["--uniq"],
                        StoreTrue,
                        "Keep only unique alignments (NH:i:1)")
            .add_option(&["--nouniq"], StoreFalse, "(default)");
        ap.refer(&mut options.out)
            .add_option(&["--out"], Store, "Output file prefix")
            .metavar("FILE");
        ap.refer(&mut options.trackline)
            .add_option(&["--trackline"],
                        StoreTrue,
                        "Output a UCSC track line (default)")
            .add_option(&["--notrackline"], StoreFalse, "");
        ap.refer(&mut options.trackname)
            .add_option(&["--trackname"], Store, "Name of track for the track line")
            .metavar("TRACKNAME");

        if ap.parse(std::env::args().collect(),
                   &mut std::io::sink(),
                   &mut std::io::sink())
            .is_err() {
            let name = if std::env::args().count() > 0 {
                std::env::args().nth(0).ok_or("Error!")?
            } else {
                "unknown".to_string()
            };
            ap.print_help(&name, &mut stderr())?;
            std::process::exit(1);
        }
    }
    options.split_strand = options.split_strand.to_ascii_lowercase();
    if options.split_strand.len() == 1 {
        options.split_strand = options.split_strand + "u";
    }
    let regex = Regex::new(r"^[usr][usr]$")?;
    if !regex.is_match(&options.split_strand) {
        return Err(format!("Invalid value for split_strand: \"{}\": values must be \
                                       one of: u s r uu us ur su ss sr ru rs rr",
                           options.split_strand)
            .into());
    }

    // read in the annotation file
    let mut intervals: Option<HashMap<String, intervaltree::IntervalTree<u8>>> = None;
    if !options.autostrand.is_empty() {
        if !Path::new(&options.autostrand).exists() {
            return Err(format!("Autostrand Bam file {} could not be found!",
                               &options.autostrand)
                .into());
        }
        let bam = match rust_htslib::bam::Reader::from_path(&options.autostrand) {
            Ok(r) => Ok(r),
            Err(_) => Err("BGZFError!"),
        }?;
        let header = bam.header();

        let mut refs: Vec<(u32, String)> = Vec::new();
        refs.resize(header.target_count() as usize, (0, "".to_string()));
        let target_names = header.target_names();
        for target_name in target_names {
            let tid = header.tid(target_name).ok_or("Error!")?;
            let target_len = header.target_len(tid).ok_or("Error!")?;
            let target_name = std::str::from_utf8(target_name)?;
            refs[tid as usize] = (target_len, target_name.to_string());
        }
        if options.fixchr {
            for r in &mut refs {
                let regex = Regex::new(r"^(chr|Zv9_)")?;
                if regex.is_match(&r.1) {
                    let refname = r.1.to_string();
                    r.1.clear();
                    r.1.push_str(&format!("chr{}", refname));
                }
            }
        }

        let mut interval_lists: HashMap<String, Vec<intervaltree::Interval<u8>>> = HashMap::new();
        let mut read = rust_htslib::bam::record::Record::new();
        while bam.read(&mut read).is_ok() {
            let chr = refs[read.tid() as usize].1.clone();
            if !interval_lists.contains_key(&chr) {
                interval_lists.insert(chr.clone(), Vec::new());
            }

            let mut exons: Vec<(i32, i32)> = Vec::new();
            cigar2exons(&mut exons, &read.cigar(), read.pos())?;

            if !exons.is_empty() {
                let interval_list = interval_lists.get_mut(&chr).ok_or("Error!")?;
                interval_list.push(intervaltree::Interval::new(read.pos() + 1,
                                                               exons[exons.len() - 1].1,
                                                               if read.is_reverse() {
                                                                   b'-'
                                                               } else {
                                                                   b'+'
                                                               }));
            }
        }
        for (chr, list) in &interval_lists {
            if intervals.is_none() {
                intervals = Some(HashMap::new());
            }
            let interval = intervals.as_mut().ok_or("Error!")?;
            interval.insert(chr.clone(), intervaltree::IntervalTree::new_from(list));
        }
    }

    // // analyze the bam file and produce histograms
    if !options.autostrand.is_empty() {
        // make both stranded and unstranded files
        analyze_bam(&options,
                    &options.split_strand,
                    !options.autostrand.is_empty(),
                    &intervals)?;
        analyze_bam(&options, "uu", false, &intervals)?;
    } else {
        analyze_bam(&options,
                    &options.split_strand,
                    !options.autostrand.is_empty(),
                    &intervals)?;
    }
    Ok(())
}

fn main() {
    let result = run();
    if result.is_err() {
        writeln!(stderr(),
                 "Backtrace: {:?}",
                 result.err().unwrap().backtrace())
            .unwrap();
        std::process::exit(1);
    }
    std::process::exit(0);
}
