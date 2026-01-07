# xamdst -- a BAM/CRAM Depth Stat Tool

This project is a hard fork of [bamdst by shiquan](https://github.com/shiquan/bamdst). It has been refactored to use HTSlib instead of the legacy samlib/bam.h, bringing support for multi-threading and CRAM format.

xamdst is a lightweight tool to stat the depth coverage of target regions of BAM/CRAM/SAM file(s).

Input files should be properly sorted, and the probe file (bed file) and the output dir

must be specified in the first place.

## Dependencies

- **htslib** - Required for BAM/CRAM/SAM file reading and writing
- **zlib** - Required for compression
- **pthread** - Required for multi-threading support

On Debian/Ubuntu:
```bash
apt-get install libhts-dev zlib1g-dev
```

On macOS with Homebrew:
```bash
brew install htslib
```

## Install
```
git clone https://github.com/pzweuj/xamdst
cd xamdst
make
```

or Docker
```
docker pull ghcr.io/pzweuj/mapping:2026Jan
```

## USAGE

Normal:

	xamdst -p <probe.bed> -o ./ in1.bam

With multi-threading (8 threads):

	xamdst -p <probe.bed> -o ./ --threads 8 in1.bam

Pipeline mode:

	samtools view in1.bam -u | xamdst -p x.bed -o ./ -

## PARAMETERS

-o / --outdir [dir]

set the output dir [mandatory]

-p / --bed [file]

the probe or captured target region file, these regions will be merged first [mandatory]

## OPTIONAL PARAMETERS

-f / --flank [num]

if you want calculate the coverage of flank region, set this value, default is 200

--maxdepth [num]

for some projects, the depths of sepcial region are very high, if you don't want show

these unnormal depths in cumulation distrbution file, set the cutoff value to filter them.

default is 0 (no filter).

--cutoffdepth [num]

for some projects, people care about the coverage of specified depth, like 10000x etc.

bamdst just calculate the coverage of 0x, 4x, 10x, 30x, 100x, so you can set this value

to show the specified coverage in the coverage.report file. Default is 0.

--isize [num]

for bad mapped paired reads, the inferred insert size is very huge. So set a cutoff

value for reasonal visual purpose. Default is 2000.

--uncover [num]

set this cutoff value for calculate the bad covered region. Default is <5.

--depthratio [ratios]

specify depth ratios for coverage calculation relative to average depth.

Use comma-separated values like "0.1,0.2,0.5". Default is "0.2,0.5".

This calculates coverage at positions with depth > (ratio × average_depth).

--threads [num]

set the number of threads for BAM/CRAM I/O operations. Uses htslib's native
multi-threading for faster decompression and compression. Default is 0 (single-threaded mode).

## OUTPUT FILES

Eight files will be created in the output direction. There are:

-**coverage.report**

-**coverage.report.json** (NEW: JSON format for programmatic access)

-**cumu.plot**

-**insert.plot**

-**chromosome.report**

-**region.tsv.gz**

-**depth.tsv.gz**

-**uncover.bed**

## DETAILS of each file

**coverage.report**

This file contains all the coverage information of target and

flank region, and reads stat information of the input file(s). 

Here is the full details of each entry.

     [Total] Raw Reads (All reads) // All reads in the bam file(s).
     [Total] QC Fail reads // Reads number failed QC, this flag is marked by other software,like bwa. See flag in the bam structure.
     [Total] Raw Data(Mb) // Total reads data in the bam file(s).
 	[Total] Paired Reads // Paired reads numbers.
	[Total] Mapped Reads // Mapped reads numbers.
	[Total] Fraction of Mapped Reads // Ratio of mapped reads against raw reads.
	[Total] Mapped Data(Mb) // Mapped data in the bam file(s).
	[Total] Fraction of Mapped Data(Mb) // Ratio of mapped data against raw data.
	[Total] Properly paired // Paired reads with properly insert size. See bam format protocol for details.
	[Total] Fraction of Properly paired // Ratio of properly paired reads against mapped reads
	[Total] Read and mate paired // Read (read1) and mate read (read2) paired.
	[Total] Fraction of Read and mate paired // Ratio of read and mate paired against mapped reads
	[Total] Singletons // Read mapped but mate read unmapped, and vice versa.
	[Total] Read and mate map to diff chr // Read and mate read mapped to different chromosome, usually because mapping error and structure variants.
	[Total] Read1 // First reads in mate paired sequencing
	[Total] Read2 // Mate reads
	[Total] Read1(rmdup) // First reads after remove duplications.
	[Total] Read2(rmdup) // Mate reads after remove duplications.
	[Total] forward strand reads // Number of forward strand reads.
	[Total] backward strand reads // Number of backward strand reads.
	[Total] PCR duplicate reads // PCR duplications.
	[Total] Fraction of PCR duplicate reads // Ratio of PCR duplications.
	[Total] Map quality cutoff value // Cutoff map quality score, this value can be set by -q. default is 20, because some variants caller like GATK only consider high quality reads.
	[Total] MapQuality above cutoff reads // Number of reads with higher or equal quality score than cutoff value.
	[Total] Fraction of MapQ reads in all reads // Ratio of reads with higher or equal Q score against raw reads.
	[Total] Fraction of MapQ reads in mapped reads // Ratio of reads with higher or equal Q score against mapped reads.
	[Target] Target Reads // Number of reads covered target region (specified by bed file).
	[Target] Fraction of Target Reads in all reads // Ratio of target reads against raw reads.
	[Target] Fraction of Target Reads in mapped reads // Ratio of target reads against mapped reads.
	[Target] Target Data(Mb) // Total bases covered target region. If a read covered target region partly, only the covered bases will be counted.
	[Target] Target Data Rmdup(Mb) // Total bases covered target region after remove PCR duplications. 
	[Target] Fraction of Target Data in all data // Ratio of target bases against raw bases.
	[Target] Fraction of Target Data in mapped data // Ratio of target bases against mapped bases.
	[Target] Len of region // The length of target regions.
	[Target] Average depth // Average depth of target regions. Calculated by "target bases / length of regions".
	[Target] Average depth(rmdup) // Average depth of target regions after remove PCR duplications.
	[Target] Coverage (>0x) // Ratio of bases with depth greater than 0x in target regions, which also means the ratio of covered regions in target regions.
	[Target] Coverage (>=4x) // Ratio of bases with depth greater than or equal to 4x in target regions.
	[Target] Coverage (>=10x) // Ratio of bases with depth greater than or equal to 10x in target regions.
	[Target] Coverage (>=30x) // Ratio of bases with depth greater than or equal to 30x in target regions.
	[Target] Coverage (>=100x) // Ratio of bases with depth greater than or equal to 100x in target regions.
	[Target] Coverage (>=Nx) // This is addtional line for user self-defined cutoff value, see --cutoffdepth
	[Target] Target Region Count // Number of target regions. In normal practise,it is the total number of exomes.
	[Target] Region covered > 0x // The number of these regions with average depth greater than 0x.
	[Target] Fraction Region covered > 0x // Ratio of these regions with average depth greater than 0x.
	[Target] Fraction Region covered >= 4x // Ratio of these regions with average depth greater than or equal to 4x.
	[Target] Fraction Region covered >= 10x // Ratio of these regions with average depth greater than or equal to 10x.
	[Target] Fraction Region covered >= 30x // Ratio of these regions with average depth greater than or equal to 30x.
	[Target] Fraction Region covered >= 100x // Ratio of these regions with average depth greater than or equal to 100x.
	[flank] flank size // The flank size will be count. 200 bp in default. Oligos could also capture the nearby regions of target regions.
	[flank] Len of region (not include target region) // The length of flank regions (target regions will not be count).
	[flank] Average depth // Average depth of flank regions.
	[flank] flank Reads // The total number of reads covered the flank regions. Note: some reads covered the edge of target regions, will be count in flank regions also. 
	[flank] Fraction of flank Reads in all reads // Ratio of reads covered in flank regions against raw reads.
	[flank] Fraction of flank Reads in mapped reads // Ration of reads covered in flank regions against mapped reads.
	[flank] flank Data(Mb) // Total bases in the flank regions.
	[flank] Fraction of flank Data in all data // Ratio of total bases in the flank regions against raw data.
	[flank] Fraction of flank Data in mapped data // Ratio of total bases in the flank regions against mapped data.
	[flank] Coverage (>0x) // Ratio of flank bases with depth greater than 0x.
	[flank] Coverage (>=4x) // Ratio of flank bases with depth greater than or equal to 4x.
	[flank] Coverage (>=10x) // Ratio of flank bases with depth greater than or equal to 10x.
	[flank] Coverage (>=30x) // Ratio of flank bases with depth greater than or equal to 30x.


**cumu.plot**

Depth distrbution for plot.

**insert.plot**

Inferred insert size distribution for plot.

**chromosome.report**

Depth and coverage information of each chromosome.

**region.tsv.gz**

For each region in probe file (in.bed), average depth, median

depth and coverage of these regions will be listed in the file.

**depth.tsv.gz**

For each position in the probe file(in.bed), three kinds depth

value will be calculated, including raw depth, rmdup depth and

the coverage depth. The raw depth is retrieved from input bam 

file(s) without any restriction. And the rmdup depth is 

calculated after remove duplicated reads and secondary alignment

reads and low map quality reads(mapQ < 20), this value is similar

with the output depth of `samtools depth`. The coverage depth is

the raw depth with consider of deletion region, so this value 

should be equal to or greated than the raw depth. We usw raw depth

to stat the coverage information in the coverage.report file. If

you want use rmdup depth to calculate the coverage, please use 

the parameter "--use_rmdup".

**uncover.bed**

This bed file contains the bad covered or uncovered region in the

input bam file(s) against the probe file. Set the cutoff value of

uncover by parameter "--uncover"

