# bamdst -- a BAM Depth Stat Tool

Bamdst is a lightweight tool to stat the depth coverage of  target regions of bam file(s).

Bam file(s) should be properly sorted, and the probe file (bed file) and the output dir

must be specified in the first place.

## USAGE

Normal:

	bamdst -p <probe.bed> -o ./ in1.bam

Pipeline mode:

	samtools view in1.bam -u | bamdst -p x.bed -o ./ -

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

--use_rmdup (an invalid parament since v1.0.0 )

Use rmdup depth instead of cover depth to calculate the coverage of target regions and

so on.

## OUTPUT FILES

Seven files will be created in the output direction. There are:

-**coverage.report**

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

# Known bugs

For a large region, like whole genome region, bamdst may go crash with a segmental fault. I have noticed issues like this, and this bug can be tolerated by split a large region into several small pieces. However, this bug may not be fixed until next major update.