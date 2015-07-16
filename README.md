bam2bedgraph
==============

Convert .bam alignment files to bedgraph or bigwig format. 

Requires:

- nightly version of Rust (use multirust to install it)

To build:

Run ```cargo build```

Help documentation:

```
Usage:
    ./target/release/bam2bedgraph [OPTIONS] BAMFILE

Convert a bam file into a bedgraph/bigwig file.

positional arguments:
  BAMFILE               Input BAM filename

optional arguments:
  -h,--help             show this help message and exit
  --split               Use CIGAR string to split alignment into separate exons
                        (default)
  --nosplit
  --autostrand ANNOT_BAMFILE
                        Attempt to determine the strandedness of the input data
                        using an annotation file. Must be an indexed .bam file.
  --strand [TYPE]       Split output bedgraph by strand: Possible values: u s r
                        uu us ur su ss sr ru rs rr, first char is read1, second
                        is read2, u=unstranded, s=stranded, r=reverse
  --read                Split output bedgraph by read number
  --noread              (default)
  --zero                Pad output bedgraph with zeroes
  --nozero              (default)
  --fixchr              Transform chromosome names to be UCSC-compatible
  --nofixchr            (default)
  --paired              Only output paired read alignments
  --nopaired            (default)
  --proper              Only output proper-paired read alignments
  --noproper            (default)
  --primary             Only output primary alignments
  --noprimary           (default)
  --bigwig              Output bigwig files (requires bedGraphToBigWig in
                        $PATH)
  --nobigwig            (default)
  --uniq                Keep only unique alignments (NH:i:1)
  --nouniq              (default)
  --out FILE            Output file prefix
  --trackline           Output a UCSC track line (default)
  --notrackline
  --trackname TRACKNAME Name of track for the track line
```
