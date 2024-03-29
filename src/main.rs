#![cfg_attr(feature = "cargo-clippy", allow(cyclomatic_complexity, trivial_regex))]
use std::collections::BTreeMap;
use std::io::Write;
use std::io::BufWriter;
use std::fs::File;
use std::str;
use std::path::{PathBuf, Path};

use regex::Regex;

use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;

use structopt::StructOpt;

use std::vec::Vec;

use rust_htslib::bam::Read;
use rust_htslib::bam::Reader;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;

use duct::cmd;
use bio::io::gff::GffType;

use lazy_static::lazy_static;
use lazy_regex::{regex as re};

use shell_words;
use rust_htslib::bam::record::Cigar;

use eyre::eyre;
pub type Report = eyre::Report<color_eyre::Context>;
pub type Result<T, E = Report> = core::result::Result<T, E>;
trait ToResult<T> {
    fn r(self) -> Result<T>;
}
impl<T> ToResult<T> for Option<T> {
    fn r(self) -> Result<T> {
        self.ok_or_else(|| eyre!("NoneError"))
    }
}

fn aligned_blocks(record: &bam::Record) -> Vec<[i64; 2]> {
    let mut result = Vec::new();
    let mut pos = record.pos();
    for entry in record.cigar().iter() {
        match entry {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                result.push([pos, pos + *len as i64]);
                pos += *len as i64;
            }
            Cigar::Del(len) => pos += *len as i64,
            Cigar::RefSkip(len) => pos += *len as i64,
            _ => (),
        }
    }
    result
}

#[derive(StructOpt, Debug)]
#[structopt(name = "bam2bedgraph", about = "Convert bam files to bedgraph/bigWig format")]
struct Options {
    #[structopt(long = "nosplit", help = "Do not use CIGAR string to split alignment into separate exons")]
    nosplit_exons: bool,
    #[structopt(long = "read", help = "Split output bedgraph by read number")]
    split_read: bool,
    #[structopt(long = "zero", help = "Pad output bedgraph with zeroes")]
    zero: bool,
    #[structopt(long = "paired", help = "Only output paired read alignments")]
    paired_only: bool,
    #[structopt(long = "proper", help = "Only output proper-paired read alignments")]
    proper_only: bool,
    #[structopt(long = "primary", help = "Only output primary read alignments")]
    primary_only: bool,
    #[structopt(long = "notrackline", help = "Do not output a UCSC track line")]
    notrackline: bool,
    #[structopt(long = "bigwig", help = "Output bigwig files (requires bedGraphToBigWig in $PATH)")]
    bigwig: bool,
    #[structopt(long = "uniq", help = "Keep only unique alignments (NH:i:1)")]
    uniq: bool,
    #[structopt(help = "Convert a bam file into a bedgraph/bigwig file.", name="BAMFILE")]
    bamfile: String,
    #[structopt(long = "trackname", help = "Name of track for the track line", name="TRACKNAME", default_value="")]
    trackname: String,
    #[structopt(long = "out", help = "Output file prefix", name="PREFIX", default_value="")]
    out: String,
    #[structopt(long = "autostrand", help = 
        "Attempt to determine the strandedness of the input data using an \
        annotation file. Must be a .gff/.gtf file.", name="AUTOSTRAND_FILE", default_value="")]
    autostrand: String,
    #[structopt(long = "split_strand", help =
        "Split output bedgraph by strand: Possible values: u s r uu us ur su ss \
        sr ru rs rr, first char is read1, second is read2, u=unstranded, \
        s=stranded, r=reverse", default_value = "uu", name="DESCRIPTION")]
    split_strand: String,
}

fn open_file(options: &Options,
             read_number: i32,
             strand: &str,
             split_strand: &str,
             fhs: &mut BTreeMap<String, Option<BufWriter<File>>>)
             -> Result<String> {
    let mut prefix = PathBuf::new();
    prefix.set_file_name(&options.bamfile);
    prefix.set_extension("");
    let track_name = vec![if !options.trackname.is_empty() {
                              options.trackname.clone()
                          } else {
                              let p = prefix.as_path().to_str().r()?;
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
                            let p = prefix.as_path().to_str().r()?;
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
        let mut f = BufWriter::new(File::create(&filename)?);
        if !options.notrackline && !options.bigwig {
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
             chr: &(u64, String),
             histogram: &BTreeMap<(i32, String), Vec<i32>>,
             fhs: &mut BTreeMap<String, Option<BufWriter<File>>>,
             split_strand: &str)
             -> Result<()> {
    for (key, histo) in histogram {
        let read_number = key.0;
        let strand = &key.1;
        let filename = open_file(options, read_number, strand, split_strand, fhs)?;
        let f = fhs.get_mut(&filename).r()?;
        let file = f.as_mut().r()?;
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
               intervals: &Option<BTreeMap<String, IntervalTree<i64, u8>>>)
               -> Result<()> {
    if !Path::new(&options.bamfile).exists() {
        Err(eyre!("Bam file {} could not be found!", &options.bamfile))?;
    }
    let mut bam = Reader::from_path(&options.bamfile)?;
    let header = bam.header().clone();
    let mut refs = vec![(0, "".to_string()); header.target_count() as usize];
    let target_names = header.target_names();
    for target_name in target_names {
        let tid = header.tid(target_name).r()?;
        let target_len = header.target_len(tid).r()?;
        let target_name = std::str::from_utf8(target_name)?;
        refs[tid as usize] = (target_len, target_name.to_string());
    }

    if autostrand_pass {
        eprintln!("Running strand detection phase on {}", options.bamfile);
    } else {
        eprintln!("Building {} histograms for {}",
                  if split_strand=="uu" {"unstranded"} else {"stranded"},
                  options.bamfile);
    }

    // build a lookup map for the refseqs
    let mut refmap: BTreeMap<String, usize> = BTreeMap::new();
    for (i, _) in refs.iter().enumerate() {
        refmap.insert(refs[i].1.to_string(), i);
    }

    let mut lastchr: i32 = -1;
    let mut fhs: BTreeMap<String, Option<BufWriter<File>>> = BTreeMap::new();
    let mut histogram: BTreeMap<(i32, String), Vec<i32>> = BTreeMap::new();

    let mut autostrand_totals: BTreeMap<char, i64> = BTreeMap::new();
    autostrand_totals.insert('s', 0);
    autostrand_totals.insert('r', 0);
    let mut autostrand_totals2: BTreeMap<char, i64> = BTreeMap::new();
    autostrand_totals2.insert('s', 0);
    autostrand_totals2.insert('r', 0);

    let mut read = bam::Record::new();
    while bam.read(&mut read)? {
        // skip unaligned reads
        if read.tid() < 0 { continue }

        // if we've hit a new chr, write out the bedgraph data and clear the histogram
        if lastchr == -1 || read.tid() != lastchr {
            if !autostrand_pass && !histogram.is_empty() && lastchr != -1 {
                write_chr(options,
                          &refs[lastchr as usize],
                          &histogram,
                          &mut fhs,
                          split_strand)?;
            }
            println!("Reading chr {}", &refs[read.tid() as usize].1);
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
            if hits == None || hits.r()?.integer() != 1 {
                continue;
            }
        }

        let mut exons: Vec<[i64; 2]> = Vec::new();
        let mut get_exons = aligned_blocks(&read);
        if options.nosplit_exons && !get_exons.is_empty() {
            let first = get_exons.get(0).r()?;
            let last = get_exons.get(get_exons.len() - 1).r()?;
            // if the exon does not have positive width, skip it
            if last[1] - first[0] <= 0 {
                continue;
            }
            exons = vec![[first[0], last[1]]];
        } else {
            exons.append(&mut get_exons);
        }
        let read_number = if read.is_last_in_template() { 2 } else { 1 };

        // attempt to determine the strandedness of the transcript
        // read numbers match, is not reverse, is not flipped
        // let xs = read.aux("XS".as_bytes());
        let strand = 
        //if xs.is_some() {
         //   str::from_utf8(xs.r()?.string())?
        //} else 
        if read_number == 1 {
                if split_strand.chars().nth(0).r()? == 'r' {
                    if read.is_reverse() { "+" } else { "-" }
                } else if split_strand.chars().nth(0).r()? == 's' {
                    if read.is_reverse() { "-" } else { "+" }
                } else {
                    ""
                }
            } else if read_number == 2 {
                if split_strand.chars().nth(1).r()? == 's' {
                    if read.is_reverse() { "-" } else { "+" }
                } else if split_strand.chars().nth(1).r()? == 'r' {
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
                    let intervals = intervals.as_ref().r()?;
                    if intervals.contains_key(&refs[lastchr as usize].1) {
                        for r in intervals[&refs[lastchr as usize].1].find(exon[0]..exon[1]) {
                            let overlap_length = std::cmp::min(exon[1], r.interval().end) -
                                                 std::cmp::max(exon[0], r.interval().start);

                            let strandtype = if read.is_reverse() == (*r.data() == b'-') {
                                's'
                            } else {
                                'r'
                            };
                            if read_number == 1 {
                                let at = autostrand_totals.get_mut(&strandtype).r()?;
                                *at += overlap_length as i64
                            } else if read_number == 2 {
                                let at2 = autostrand_totals2.get_mut(&strandtype).r()?;
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
                if ref_length < exon[1] as u64 {
                    refs[read.tid() as usize].0 = exon[1] as u64;
                }
                if histogram[&tuple].len() < ref_length as usize {
                    let h = histogram.get_mut(&tuple).r()?;
                    h.resize(ref_length as usize, 0);
                }

                let h = histogram.get_mut(&tuple).r()?;
                for pos in std::cmp::max(0i64, exon[0])..std::cmp::min(ref_length as i64, exon[1])
                {
                    h[pos as usize] += 1;
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
    for fh in &mut fhs {
        eprintln!("Wrote file: {}", fh.0);
        *fh.1 = None;
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
                eprintln!(
                         "Total evidence for read 1 strand type {}: {}",
                         i.0,
                         i.1);
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
                eprintln!("Total evidence for read 2 strand type {}: {}",
                         i.0,
                         i.1);
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
        eprintln!(
                 "autostrand_pass found best strand type: {}",
                 best_strand);

        // re-run analyzeBam with the strand type indicated
        analyze_bam(options, &best_strand, false, intervals)?;
    }

    if !autostrand_pass && options.bigwig {
        // Convert bedgraph file to bigwig file
        for fh in &fhs {
            // write the genome file for bigwigs
            let genome_filename = format!("{}.genome", fh.0);
            {
                let mut genome_fh = BufWriter::new(File::create(&genome_filename)?);
                for r in &refs {
                    writeln!(genome_fh, "{}\t{}", r.1, r.0)?;
                }
            }

            // run bedGraphToBigWig
            let bigwig_file = String::from(re!(r"\.bedgraph$").replace(fh.0, ".bw"));
            let sorted_bedgraph = String::from(re!(r"\.bedgraph$").replace(fh.0, ".sorted.bedgraph"));
            let sort_cmd = ["sort","-k1,1","-k2,2n","-o", &sorted_bedgraph, fh.0];
            eprintln!("Running command: {}", shell_words::join(&sort_cmd));
            cmd(sort_cmd[0], &sort_cmd[1..]).env("LC_COLLATE","C").run()?;
            let convert_cmd = ["bedGraphToBigWig",&sorted_bedgraph,&genome_filename,&bigwig_file];
            eprintln!("Running command: {}", shell_words::join(&convert_cmd));
            cmd(convert_cmd[0], &convert_cmd[1..]).run()?;
            eprintln!("Wrote file: {}", bigwig_file);
            // remove the bedgraph file
            eprintln!("Removing file: {}", fh.0);
            std::fs::remove_file(fh.0)?;
            // remove the sorted bedgraph file
            eprintln!("Removing file: {}", sorted_bedgraph);
            std::fs::remove_file(&sorted_bedgraph)?;
            // remove the genome file
            eprintln!("Removing file: {}", genome_filename);
            std::fs::remove_file(&genome_filename)?;
        }
    };
    Ok(())
}

fn run() -> Result<()> {
    let mut options = Options::from_args();
    options.split_strand = options.split_strand.to_lowercase();
    if options.split_strand.len() == 1 {
        options.split_strand = options.split_strand + "u";
    }
    if !re!(r"^[usr][usr]$").is_match(&options.split_strand) {
        Err(eyre!("Invalid value for split_strand: \"{}\": values must be \
                                       one of: u s r uu us ur su ss sr ru rs rr",
                           options.split_strand))?;
    }

    // read in the annotation file
    let mut intervals: Option<BTreeMap<String, IntervalTree<i64, u8>>> = None;
    if !options.autostrand.is_empty() {
        if !Path::new(&options.autostrand).exists() {
            Err(eyre!("Autostrand Bam file {} could not be found!", &options.autostrand))?;
        }
        let mut interval_lists: BTreeMap<String, Vec<(Interval<i64>, u8)>> = BTreeMap::new();
        let format = if re!(r"(?i)\.gtf").is_match(&options.autostrand) { GffType::GTF2 } else { GffType::GFF3 };
        let mut reader = bio::io::gff::Reader::from_file(&options.autostrand, format)?;

        eprintln!("Loading autostrand file: {}", options.autostrand);
        for record in reader.records() {
            let record = record?;
            let chr = record.seqname();
            if !interval_lists.contains_key(chr) {
                interval_lists.insert(chr.to_string(), Vec::new());
            }

            let interval_list = interval_lists.get_mut(chr).r()?;
            interval_list.push((Interval::new((*record.start()-1) as i64..*record.end() as i64)?,
                                if record.strand().is_some() &&
                                   record.strand().r()?.strand_symbol() == "-"
                                { b'-' } else { b'+' }));
        }
        for (chr, list) in &interval_lists {
            if intervals.is_none() {
                intervals = Some(BTreeMap::new());
            }
            let interval = intervals.as_mut().r()?;
            let mut tree = IntervalTree::<i64, u8>::new();
            for l in list {
                tree.insert(l.0.clone(), l.1);
            }
            interval.insert(chr.clone(), tree);
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

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "full");
    std::env::set_var("RUST_LIB_BACKTRACE", "1");
    Ok(run()?)
}
