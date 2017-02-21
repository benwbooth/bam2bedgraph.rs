bam2bedgraph
==============

Convert .bam alignment files to bedgraph or bigwig format. 

Requires:

- latest Rust stable (or nightly)

To build:

Run ```cargo build```

Help documentation:

```
bam2bedgraph 1.0.0
Ben Booth <benwbooth@gmail.com>
Convert bam files to bedgraph/bigWig format

USAGE:
    bam2bedgraph [FLAGS] [OPTIONS] <BAMFILE>

FLAGS:
        --bigwig     Output bigwig files (requires bedGraphToBigWig in $PATH)
        --fixchr     Transform chromosome names to be UCSC-compatible
    -h, --help       Prints help information
        --paired     Only output paired read alignments
        --primary    Only output primary read alignments
        --proper     Only output proper-paired read alignments
        --read       Split output bedgraph by read number
        --uniq       Keep only unique alignments (NH:i:1)
    -V, --version    Prints version information
        --zero       Pad output bedgraph with zeroes

OPTIONS:
        --autostrand <AUTOSTRAND_FILE>    Attempt to determine the strandedness of the input data using an annotation file. Must be a .bam file. [default: ]
        --split_strand <DESCRIPTION>      Split output bedgraph by strand: Possible values: u s r uu us ur su ss sr ru rs rr, first char is read1, second is
                                          read2, u=unstranded, s=stranded, r=reverse [default: uu]
        --out <PREFIX>                    Output file prefix [default: ]
        --trackname <TRACKNAME>           Name of track for the track line [default: ]
        --split <split_exons>             Use CIGAR string to split alignment into separate exons (default) [default: true]
        --trackline <trackline>           Output a UCSC track line (default) [default: true]

ARGS:
    <BAMFILE>    Convert a bam file into a bedgraph/bigwig file.
```
