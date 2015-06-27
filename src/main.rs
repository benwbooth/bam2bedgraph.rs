use std::ascii::AsciiExt;
use std::collections::HashMap;
use std::io::Write;
use std::io::stderr;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::cmp::max;
use std::str;
use std::path::PathBuf;

extern crate argparse;
use argparse::{ArgumentParser, StoreTrue, StoreFalse, Store};

extern crate regex;
use regex::Regex;

extern crate interval;
use interval::interval_set::IntervalSet;
use interval::ops::Range;
use std::ops::Add;

extern crate rust_htslib;
use rust_htslib::bam::Reader;
use rust_htslib::bam::Read;

// void cigar2exons(vector<pair<size_t,size_t>> &exons, const vector<CigarOp> &cigar, size_t pos) {
//   for (auto &op : cigar) {
//     if (op.Type == 'M') {
//       pos += op.Length;
//       exons.push_back({pos - op.Length, pos});
//     }
//     else if (op.Type == 'N' || op.Type == 'D') {
//       pos += op.Length;
//     }
//     else if (op.Type == 'I' || op.Type == 'S' || op.Type == 'H' || op.Type == 'P') {}
//     else {
//       throw string("Bad CIGAR string: ")+op.Type;
//     }
//   }
// }

fn open_file(
    options: &Options,
    read_number: i32,
    strand: &str,
    split_strand: &str,
    fhs: &mut HashMap<String,File>)
    -> String
{
    let mut prefix = PathBuf::new();
    prefix.set_file_name(&options.bamfile);
    prefix.set_extension("");
    let track_name = vec!(
        if !options.trackname.is_empty() {options.trackname.clone()} else {prefix.as_path().to_str().unwrap().to_string()},
        if options.split_read && read_number>0 {format!(".r{}", read_number)} else {"".to_string()},
        if options.split_strand != "uu" && !strand.is_empty() {format!(".{}", strand)} else {"".to_string()},
        ).connect("");

    let filename = vec!(
        if !options.out.is_empty() {options.out.clone()} else {prefix.as_path().to_str().unwrap().to_string()},
        if options.split_read && read_number>0 {format!(".r{}", read_number)} else {"".to_string()},
        if options.split_strand != "uu" && !strand.is_empty() {format!(".{}", strand)} else {"".to_string()},
        ".bedgraph".to_string(),
        ).connect("");

    // initialize the file if needed
    if !fhs.contains_key(&filename) {
        let mut f = File::create(&filename).unwrap();
        if options.trackline {
            writeln!(f, "track type=bedGraph name=\"{}\" description=\"{}\" visibility-full", track_name, track_name).unwrap();
        }
        fhs.insert(filename.clone(), f);
    }
    filename
}

fn write_chr(
    options: &Options,
    chr: &(u32, String),
    histogram: &HashMap<(i32,String),Vec<i32>>,
    fhs: &mut HashMap<String,File>,
    split_strand: &str)
{
   for (key, histo) in histogram {
       let read_number = key.0;
       let strand = &key.1;
       let filename = open_file(options, read_number, strand, split_strand, fhs);
// 
//     string filename = open_file(read_number, strand, split_strand, fhs);
//     auto &fh = *fhs[filename];
//     // scan the histogram to produce the bedgraph data
//     size_t start = 0;
//     size_t end = 0;
//     size_t ref_length = numeric_cast<size_t>(chr.RefLength);
//     while (start < ref_length) {
//       while ((end < histo.size()? histo[end] : 0) == (start < histo.size()? histo[start] : 0) && end < ref_length) ++end;
//       if (zero || (start < histo.size()? histo[start] : 0)) {
//         fh << chr.RefName << "\t" <<
//               start << "\t" <<
//               end << "\t" <<
//               (strand == "-"? -histo[start] : histo[start]) << endl;
//       }
//       start = end;
//     }
   }
}

// if autostrand is set, read the input annotation file into an interval tree
fn read_annotation_file(autostrand: &str, intervals: &mut HashMap<(String, char), IntervalSet<i32>>) {
    if ! autostrand.is_empty() {
        let mut chr_col=0;
        let mut start_col=0;
        let mut stop_col=0;
        let mut strand_col=0;
        let mut base=0;
        for cap in Regex::new(r"\.([^.]+)$").unwrap().captures_iter(autostrand) {
            if Regex::new(r"(?i)^bed$").unwrap().is_match(cap.at(1).unwrap_or("")) {
                writeln!(&mut stderr(), "Reading bed annotation file {}", autostrand).unwrap();
                chr_col = 0;
                start_col = 1;
                stop_col = 2;
                strand_col = 5;
                base = 0;
            }
            else if Regex::new(r"(?i)^(gff|gtf|gff3)$").unwrap().is_match(cap.at(1).unwrap_or("")) {
                writeln!(&mut stderr(), "Reading gff annotation file {}", autostrand).unwrap();
                chr_col = 0;
                start_col = 3;
                stop_col = 4;
                strand_col = 6;
                base = 1;
            }
            else {
                panic!("Could not recognize file extension of annotation file {}: {}",
                       autostrand, cap.at(1).unwrap_or(""));
            }
        }
        let annot_fh = BufReader::new(File::open(autostrand).unwrap());
        let comment_line = Regex::new(r"^[ ]*#").unwrap();
        for line in annot_fh.lines() {
            let line = line.unwrap();
            if comment_line.is_match(&line) { continue }
            let cols: Vec<&str> = Regex::new(r"\t").unwrap().split(&line).collect();
            let max_col = [chr_col, start_col, stop_col, strand_col].iter().fold(0, |a, &b| max(a, b));
            if max_col < cols.len() {
                let chr = cols[chr_col].to_string();
                let start = cols[start_col].trim().parse::<i32>().unwrap();
                let stop = cols[stop_col].trim().parse::<i32>().unwrap();
                let strand_str = cols[strand_col];
                let strand =
                    if strand_str == "-1" || strand_str == "-" {'-'} else
                    if strand_str == "1" || strand_str == "+" {'+'} else {'\0'};
                if strand != '\0' {
                    let key = (chr, strand);
                    let interval = IntervalSet::new(start-(base-1), stop);
                    if !intervals.contains_key(&key) {
                        intervals.insert(key, interval);
                    }
                    else {
                        let ref intervals = intervals[&key];
                        intervals.add(interval);
                    }
                }
            }
        }
    }
}

fn analyze_bam(
    options: &Options,
    split_strand: &str,
    autostrand_pass: bool,
    intervals: &HashMap<(String, char), IntervalSet<i32>>)
{
    let bam = rust_htslib::bam::Reader::new(&options.bamfile);
    let ref header = bam.header;
    let mut refs: Vec<(u32, String)> = Vec::new();
    for target_name in header.target_names() {
        refs.push((header.tid(&target_name).unwrap(), str::from_utf8(target_name).unwrap().to_string()));
    }
    refs.sort_by(|a, b| a.0.cmp(&b.0));

    if options.fixchr {
        for r in &mut refs {
            if Regex::new(r"^(chr|Zv9_)").unwrap().is_match(&r.1) {
                let refname = r.1.to_string();
                r.1.clear();
                r.1.push_str(&format!("chr{}", refname));
            }
        }
    }
    if autostrand_pass {
        writeln!(stderr(), "Running strand detection phase on {}", options.bamfile).unwrap();
    }
    else {
        writeln!(stderr(), "Building histograms for {}", options.bamfile).unwrap();
    }

    // build a lookup map for the refseqs
    let mut refmap: HashMap<String,usize> = HashMap::new();
    for i in 0..refs.len() {
        refmap.insert(refs[i].1.to_string(), i);
    }

    let mut lastchr: i32 = -1;
    let mut fhs: HashMap<String, File> = HashMap::new();
    let mut histogram: HashMap<(i32, String), Vec<i32>> = HashMap::new();

    let mut autostrandTotals: HashMap<char, i64> = HashMap::new();
    let mut autostrandTotals2: HashMap<char, i64> = HashMap::new();

    let mut read = rust_htslib::bam::record::Record::new();
    while bam.read(&mut read).is_ok() {
        // if we've hit a new chr, write out the bedgraph data and clear the histogram
        if lastchr == -1 || read.tid() != lastchr {
            if !autostrand_pass && !histogram.is_empty() && lastchr != -1 {
                write_chr(options, &refs[lastchr as usize], &histogram, &mut fhs, split_strand);
            }
            histogram.clear();
            lastchr = read.tid();
        }

    }
//   while (bam.GetNextAlignment(read)) {
// 
//     // skip this read if it's no good
//     bool paired = read.IsPaired();
//     bool proper = read.IsProperPair();
//     bool primary  = read.IsPrimaryAlignment();
//     if ((paired_only && !paired) || (primary_only && !primary) || (proper_only && !proper)) continue;
// 
//     // skip if it's not unique and we want unique alignments
//     if (uniq) {
//       int32_t hits=0;
//       if (!read.GetTag("NH",hits) || hits != 1) continue;
//     }
// 
//     vector<pair<size_t,size_t>> exons;
//     if (split_exons) cigar2exons(exons, read.CigarData, read.Position);
//     else exons.push_back({read.Position, read.GetEndPosition()});
//     int32_t read_number = read.IsSecondMate()? 2 : 1;
// 
//     // attempt to determine the strandedness of the transcript
//     uint8_t xs = 0;
//     // read numbers match, is not reverse, is not flipped
//     string strand =
//         read.GetTag("XS",xs) && xs? string{static_cast<char>(xs)} :
//           read_number == 1? split_strand[0] == 'r'? read.IsReverseStrand()? "+" : "-" :
//                             split_strand[0] == 's'? read.IsReverseStrand()? "-" : "+" :
//                             "" :
//           read_number == 2? split_strand[1] == 's'? read.IsReverseStrand()? "-" : "+" :
//                             split_strand[1] == 'r'? read.IsReverseStrand()? "+" : "-" :
//                             "" :
//         "";
// 
//     int32_t read_num = split_read? read_number : 0;
// 
//     size_t ref_length = numeric_cast<size_t>(refs[read.RefID].RefLength);
//     // add the read to the histogram
//     for (auto &exon : exons) {
//       // try to determine the strandedness of the data
//       if (autostrandPass) {
//         if (intervals.count(refs[lastchr].RefName)) {
//           vector<Interval<int8_t>> overlappingAnnot;
//           intervals[refs[lastchr].RefName].findOverlapping(exon.first+1, exon.second, overlappingAnnot);
//           for (const auto& interval : overlappingAnnot) {
//             size_t overlap_length = std::min(exon.second, numeric_cast<size_t>(interval.stop)) -
//                                     std::max(exon.first, numeric_cast<size_t>(interval.start-1));
// 
//             char strandtype = read.IsReverseStrand() == (interval.value == '-')? 's' : 'r';
//             if (read_number == 1) autostrandTotals[strandtype] += overlap_length;
//             else if (read_number == 2) autostrandTotals2[strandtype] += overlap_length;
//           }
//         }
//       }
//       else {
//         auto tuple = std::make_tuple(read_num, strand);
//         if (!histogram.count(tuple)) histogram.insert({tuple, vector<int32_t>{}});
//         // keep track of chromosome sizes
//         if (ref_length < exon.second) refs[read.RefID].RefLength = numeric_cast<int32_t>(exon.second);
//         if (histogram[tuple].size() < ref_length) histogram[tuple].resize(ref_length);
// 
//         for (size_t pos=exon.first; pos < exon.second; ++pos) {
//           histogram[tuple][pos]++;
//         }
//       }
//     }
//   }
//   bam.Close();
// 
//   if (!autostrandPass && !histogram.empty() && lastchr != -1) {
//     write_chr(refs[lastchr], histogram, fhs, split_strand);
//   }
// 
//   // make sure empty files were created
//   if (histogram.empty() && !autostrandPass) {
//     for (auto &read_number : split_read? vector<int32_t>{1,2} : vector<int32_t>{0}) {
//       for (auto &s : split_strand != "uu"? vector<string>{"+","-"} : vector<string>{""}) {
//         open_file(read_number, s, split_strand, fhs);
//       }
//     }
//   }
// 
//   // close the filehandles
//   for (auto &fh : fhs) {
//     fh.second->close();
//   }
//   if (autostrandPass) {
//     // get the read 1 and read2 totals
//     int64_t total1 = 0;
//     int64_t total2 = 0;
//     for (auto &i : autostrandTotals) {
//       total1 += i.second;
//     }
//     for (auto &i : autostrandTotals2) {
//       total2 += i.second;
//     }
//     // figure out the best and second-best strand types for reads 1 and 2
//     pair<int32_t,int64_t> best1;
//     pair<int32_t,int64_t> second_best1;
//     pair<int32_t,int64_t> best2;
//     pair<int32_t,int64_t> second_best2;
//     for (auto &i : autostrandTotals) {
//       cerr << "Total evidence for read 1 strand type " << i.first << ": " << i.second << endl;
//       if (best1.second < i.second) {
//         second_best1 = best1;
//         best1 = i;
//       }
//       else if (second_best1.second < i.second) {
//         second_best1 = i; 
//       }
//     }
//     for (auto &i : autostrandTotals2) {
//       cerr << "Total evidence for read 2 strand type " << i.first << ": " << i.second << endl;
//       if (best2.second < i.second) {
//         second_best2 = best2;
//         best2 = i;
//       }
//       else if (second_best2.second < i.second) {
//         second_best2 = i;
//       }
//     }
//     const double threshold = 0.0; //threshold isn't working, set to zero
//     char strand1 = (total1 > 0.0)?
//       threshold < numeric_cast<double>(best1.second - second_best1.second) / numeric_cast<double>(total1)? best1.first : 'u'
//                   : 'u';
//     char strand2 = (total2 > 0.0)?
//       threshold < numeric_cast<double>(best2.second - second_best2.second) / numeric_cast<double>(total2)? best2.first : 'u'
//                   : 'u';
//     string best_strand {strand1, strand2};
//     cerr << "autostrandPass found best strand type: " << best_strand << endl;
// 
//     // re-run analyzeBam with the strand type indicated
//     analyzeBam(best_strand, false, intervals);
//   }
//   if (!autostrandPass && bigwig) {
//     // Convert bedgraph file to bigwig file
//     for (auto &fh : fhs) {
//       // write the genome file for bigwigs
//       string genome_filename = fh.first+".genome";
//       ofstream genome_fh {};
//       genome_fh.exceptions(std::ios::badbit | std::ios::failbit);
//       genome_fh.open(genome_filename);
//       for (auto& ref : refs) {
//         genome_fh << ref.RefName << "\t" << ref.RefLength << endl;
//       }
//       genome_fh.close();
// 
//       string cmd = (format("bedGraphToBigWig '%s' '%s' '%s'") %
//           boost::regex_replace(fh.first, regex(R"(')"), "'\\''") %
//           boost::regex_replace(genome_filename, regex(R"(')"), "'\\''") %
//           boost::regex_replace(boost::regex_replace(fh.first, regex(R"(\.bedgraph$)"),"")+".bw", regex(R"(')"), "'\\''")).str();
//       int result = system(cmd.c_str());
//       if (WEXITSTATUS(result) != 0) {
//         throw format("Command \"%s\" returned bad exit status: %d") % cmd % WEXITSTATUS(result);
//       }
// 
//       // remove the bedgraph file
//       boost::filesystem::remove(fh.first);
//       // remove the genome file
//       boost::filesystem::remove(genome_filename);
//     }
//   }
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
    split_strand: String
}

fn main() {
    let mut options: Options = Options {
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
        split_strand: "uu".to_string()
    };
    {
        let mut ap = ArgumentParser::new();
        ap.set_description("Convert a bam file into a bedgraph/bigwig file.");
        ap.refer(&mut options.bamfile).
            add_argument("BAMFILE", Store, "Input BAM filename").required();
        ap.refer(&mut options.split_exons).
            add_option(&["--split"], StoreTrue, "Use CIGAR string to split alignment into separate exons (default)").
            add_option(&["--nosplit"], StoreFalse, "");
        ap.refer(&mut options.autostrand).
            add_option(&["--autostrand"], Store, "Attempt to determine the strandedness of the input data using an annotation file. Can take GFF2/3, GTF, BED formatted files").
            metavar("ANNOT_FILE");
        ap.refer(&mut options.split_strand).
            add_option(&["--strand"], Store, "Split output bedgraph by strand: Possible values: u s r uu us ur su ss sr ru rs rr, first char is read1, second is read2, u=unstranded, s=stranded, r=reverse").
            metavar("[TYPE]");
        ap.refer(&mut options.split_read).
            add_option(&["--read"], StoreTrue, "Split output bedgraph by read number").
            add_option(&["--noread"], StoreFalse, "(default)");
        ap.refer(&mut options.zero).
            add_option(&["--zero"], StoreTrue, "Pad output bedgraph with zeroes").
            add_option(&["--nozero"], StoreFalse, "(default)");
        ap.refer(&mut options.fixchr).
            add_option(&["--fixchr"], StoreTrue, "Transform chromosome names to be UCSC-compatible").
            add_option(&["--nofixchr"], StoreFalse, "(default)");
        ap.refer(&mut options.paired_only).
            add_option(&["--paired"], StoreTrue, "Only output paired read alignments").
            add_option(&["--nopaired"], StoreFalse, "(default)");
        ap.refer(&mut options.proper_only).
            add_option(&["--proper"], StoreTrue, "Only output proper-paired read alignments").
            add_option(&["--noproper"], StoreFalse, "(default)");
        ap.refer(&mut options.primary_only).
            add_option(&["--primary"], StoreTrue, "Only output primary alignments").
            add_option(&["--noprimary"], StoreFalse, "(default)");
        ap.refer(&mut options.bigwig).
            add_option(&["--bigwig"], StoreTrue, "Output bigwig files (requires bedGraphToBigWig in $PATH)").
            add_option(&["--nobigwig"], StoreFalse, "(default)");
        ap.refer(&mut options.uniq).
            add_option(&["--uniq"], StoreTrue, "Keep only unique alignments (NH:i:1)").
            add_option(&["--nouniq"], StoreFalse, "(default)");
        ap.refer(&mut options.out).
            add_option(&["--out"], Store, "Output file prefix").metavar("FILE");
        ap.refer(&mut options.trackline).
            add_option(&["--trackline"], StoreTrue, "Output a UCSC track line (default)").
            add_option(&["--notrackline"], StoreFalse, "");
        ap.refer(&mut options.trackname).
            add_option(&["--trackname"], Store, "Name of track for the track line").metavar("TRACKNAME");
        ap.parse_args_or_exit();
    }
    options.split_strand = options.split_strand.to_ascii_lowercase();
    if options.split_strand.len() == 1 {
        options.split_strand = options.split_strand + "u";
    }
    if ! Regex::new(r"^[usr][usr]$").unwrap().is_match(&options.split_strand) {
        panic!("Invalid value for split_strand: \"{}\": values must be one of: u s r uu us ur su ss sr ru rs rr", options.split_strand);
    }

    // read in the annotation file
    let mut intervals: HashMap<(String, char), IntervalSet<i32>> = HashMap::new();
    if ! options.autostrand.is_empty() {
        read_annotation_file(&options.autostrand, &mut intervals);
    }

    // // analyze the bam file and produce histograms
    if ! options.autostrand.is_empty() {
        // make both stranded and unstranded files
        analyze_bam(&options, &options.split_strand, !options.autostrand.is_empty(), &intervals);
        analyze_bam(&options, "uu", false, &intervals);
    }
    else {
         analyze_bam(&options, &options.split_strand, !options.autostrand.is_empty(), &intervals);
    }
}

