bamdst -- a BAM Depth Stat Tool
================================
Bamdst is a lightweight tool to stat the depth coverage of  target regions of bam file(s).
Bam file(s) should be properly sorted, and the probe file (bed file) and the output dir
must be specified in the first place.

USAGE
------

Normal:

	bamdst -p <probe.bed> -o ./ in1.bam

Pipeline mode:

	samtools view in1.bam -u | bamdst -p x.bed -o ./ -

PARAMETERS
-----------

-o / --outdir [dir]

set the output dir [mandatory]

-p / --bed [file]

the probe or captured target region file, these regions will be merged first [mandatory]

OPTIONAL PARAMETERS
-------------------

-f / --flank [num]

if you want calculate the coverage of flank region, set this value, default is 200

--maxdepth [num]

for some projects, the depths of sepcial region are very high, if you don't want show
these unnormal depths in cumulation distrbution file, set the cutoff value to filter them.
default is 0 (no filter).

--cutoffdepth [num]

for some projects, people care about the coverage of specified depth, like 10000x etc.
bamdst just calculate the coverage of 0x, 4x, 10x, 30x, 100x, so you can set this value
to show the specified coverage in the coverage.report file. default is 0.

--isize [num]

for bad mapped paired reads, the inferred insert size is very huge. So set a cutoff
value for reasonal visual purpose. default is 2000.

OUTPUT FILES
------------
Seven files will be created in the output direction. There are:

* ** coverage.report **
* ** cumu.plot **
* ** insert.plot **
* ** chromosome.report **
* ** region.tsv.gz **
* ** depth.tsv.gz **
* ** uncover.bed **
