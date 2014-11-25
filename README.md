bamdst -- a BAM Depth Stat Tool
================================
Bamdst is a lightweight tool to stat the depth coverage of 
target regions of bam file(s).
Bam file(s) should be properly sorted, and the probe file (bed
file) and the output dir must be specified in the first place.

USAGE
------

Normal:

	- bamdst -p <probe.bed> -o . in1.bam

Pipeline mode:

	- samtools view in1.bam -u | bamdst -p x.bed -o . -

PARAMETERS
-----------

-o / --outdir

set the output dir [mandatory]

-p / --bed

the probe or captured target region file, these regions will be merged
first [mandatory]



