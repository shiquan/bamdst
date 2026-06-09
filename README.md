# bamdst — BAM Depth Statistics for Targeted Sequencing

bamdst is a lightweight tool to compute depth and coverage statistics for target
regions (e.g. exome capture panels) from sorted BAM files.

## Background

bamdst was originally developed as a personal learning project when I first started
working with BAM files and the samtools codebase. At the time, reading through
samtools' BAM parsing routines was one of the best ways to understand how BAM files
are structured at the binary level, including the BGZF block format, the BAM record
layout, the CIGAR encoding, and the coordinate-sorted traversal pattern.

Over the years, samtools and htslib have evolved into powerful, production-grade
libraries. However, they have also grown substantially in complexity, adding CRAM
support, multi-threading, remote file access, and elaborate indexing schemes. For a
beginner trying to understand the fundamentals of BAM file I/O, reading the modern
htslib source is no longer the approachable experience it once was.

This project takes a deliberate, opinionated stance.

- **No htslib dependency.** bamdst bundles a minimal BAM parsing library (`samlib/`)
  and a BGZF implementation (`bgzf.c`) that are small enough to read and understand
  in an afternoon. There are no external build requirements beyond a C compiler and
  zlib.
- **No multi-threading.** The processing loop is single-threaded and sequential,
  mirroring the natural sorted-BAM traversal pattern. This makes the control flow
  easy to trace and debug.
- **Self-contained.** The entire codebase is ~6500 lines of C, including all
  dependencies. You can read it from top to bottom.

## Install

```bash
git clone https://github.com/shiquan/bamdst
cd bamdst
make
```

bamdst is also hosted on conda via a third-party package `xdgene::bamdst` (thanks @xdgene), but may not be the latest version.

## Usage

```bash
# Normal mode
bamdst -p probe.bed -o output_dir/ in1.bam

# Multiple BAM files
bamdst -p probe.bed -o output_dir/ in1.bam in2.bam

# Pipeline mode (stdin)
samtools view in1.bam -u | bamdst -p probe.bed -o output_dir/ -

# CSV output
bamdst -F csv -p probe.bed -o output_dir/ in1.bam

# JSON output
bamdst -F json -p probe.bed -o output_dir/ in1.bam

# With custom depth ratio thresholds
bamdst --depthratio 0.1,0.3,0.5 -p probe.bed -o output_dir/ in1.bam
```

## Parameters

### Mandatory

| Option | Description |
|--------|-------------|
| `-o`, `--outdir DIR` | Output directory. Must already exist and be writable. bamdst will not create it. |
| `-p`, `--bed FILE` | Probe/target regions in BED format (0-based). Regions are merged before calculation. |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-f`, `--flank` | 200 | Number of flanking bases to extend each target region. |
| `-q` | 20 | MAPQ cutoff. Reads with MAPQ below this value are excluded from rmdup depth and `n_qual` counts. |
| `--maxdepth` | 0 (no limit) | Maximum depth value to include in the cumulative depth distribution plot. |
| `--cutoffdepth` | - | Comma-separated custom depth thresholds (e.g. `50,200,500`). Adds coverage lines like `Coverage (>=50x)`. Max 10 values. |
| `--depthratio` | - | Comma-separated ratios of average depth (e.g. `0.2,0.5`). Adds coverage lines like `Coverage (>0.20*Avg)`. Values must be in (0, 1]. |
| `--isize` | 2000 | Maximum inferred insert size for the insert size distribution plot. Larger values are excluded for visual clarity. |
| `--uncover` | 5 | Depth threshold for uncover regions. Positions with coverage depth below this value are reported in `uncover.bed`. |
| `--bamout FILE` | - | Export reads overlapping target regions to a BAM file. |
| `-F`, `--format` | txt | Output report format, one of `txt` (default), `csv`, or `json`. |
| `-1` | - | Treat the BED file as 1-based instead of 0-based. |
| `-h`, `--help` | - | Print help message. |

## Output Files

| File | Format | Description |
|------|--------|-------------|
| `coverage.report` | TXT (or `.csv` / `report.json`) | Full coverage and read statistics report (see below). |
| `chromosomes.report` | TXT (or `.csv`) | Per-chromosome depth and coverage summary table. |
| `depth_distribution.plot` | TXT (or `.csv`) | Depth distribution data for plotting (depth, count, fraction, cumulative). |
| `insertsize.plot` | TXT (or `.csv`) | Inferred insert size distribution for plotting. |
| `depth.tsv.gz` | BGZF-compressed TSV | Per-base depth values for every target position (chromosome, position, raw depth, rmdup depth, coverage depth). |
| `region.tsv.gz` | BGZF-compressed TSV | Per-region summary (chromosome, start, end, mean depth, median, coverage %, coverage-fixed %). |
| `uncover.bed` | BED | Regions with coverage depth below `--uncover` threshold. |

When using `-F csv`, `.report` and `.plot` files are replaced with `.csv` equivalents.
When using `-F json`, a single `report.json` is produced containing all statistics.

## coverage.report Field Reference

### `[Total]` Global Read Statistics

| Field | Description |
|-------|-------------|
| Raw Reads (All reads) | Total number of reads across all input BAM files. |
| QC Fail reads | Reads flagged with `BAM_FQCFAIL` (0x200). Marked by aligners like BWA for low-quality reads. |
| Raw Data(Mb) | Total bases from all reads (read length × read count) in megabases. |
| Paired Reads | Reads with the paired-end flag (`BAM_FPAIRED`, 0x1). |
| Mapped Reads | Reads that are mapped (`!BAM_FUNMAP`, i.e., not unmapped). |
| Fraction of Mapped Reads | `n_mapped / n_reads × 100%`. Overall mapping rate. |
| Mapped Data(Mb) | Total bases from mapped reads. |
| Fraction of Mapped Data(Mb) | `mapped_bases / total_bases × 100%`. |
| Properly paired | Reads where both mates are mapped with proper insert size and orientation (`BAM_FPROPER_PAIR`). |
| Fraction of Properly paired | `n_proper_pair / n_reads × 100%`. |
| Read and mate paired | Reads where both the read and its mate are mapped (regardless of orientation). |
| Fraction of Read and mate paired | `n_pair_map / n_reads × 100%`. |
| Singletons | Reads that are mapped but whose mate is unmapped (or vice versa). |
| Read and mate map to diff chr | Read pairs mapped to different chromosomes, which may indicate structural variants or mapping errors. |
| Read1 | Count of first-in-pair reads (`BAM_FREAD1`). |
| Read2 | Count of second-in-pair reads (`BAM_FREAD2`). |
| Read1(rmdup) | Read1 count after excluding duplicates, secondary alignments, and low-MAPQ reads. |
| Read2(rmdup) | Read2 count after same filtering as Read1(rmdup). |
| forward strand reads | Reads mapped to the forward strand. |
| backward strand reads | Reads mapped to the reverse strand. |
| PCR duplicate reads | Reads flagged as PCR/optical duplicates (`BAM_FDUP`, 0x400). |
| Fraction of PCR duplicate reads | `n_dup / n_mapped × 100%`. Library complexity indicator. |
| Map quality cutoff value | The MAPQ threshold set by `-q` (default 20). |
| MapQuality above cutoff reads | Number of mapped reads with MAPQ ≥ the cutoff. |
| Fraction of MapQ reads in all reads | `n_qual / n_reads × 100%`. |
| Fraction of MapQ reads in mapped reads | `n_qual / n_mapped × 100%`. |

### `[Insert size]` Inferred Insert Size

| Field | Description |
|-------|-------------|
| Average | Mean inferred insert size from properly paired reads (up to `--isize` limit). `N/A` if no paired reads present. |
| Median | Median inferred insert size. |

### `[Target]` Target Region Statistics

These statistics use the **coverage depth** (`covdep`), which counts all reads plus
positions covered by deletions (CIGAR `D` operations).

| Field | Description |
|-------|-------------|
| Target Reads | Number of reads whose aligned bases overlap with target regions. A read spanning multiple target regions is counted once. |
| Fraction of Target Reads in all reads | `n_tgt / n_reads × 100%`. Enrichment efficiency indicator. |
| Fraction of Target Reads in mapped reads | `n_tgt / n_mapped × 100%`. On-target rate. |
| Target Data(Mb) | Total aligned bases falling within target regions (megabases). |
| Target Data Rmdup(Mb) | Target bases after duplicate removal (rmdup depth × bases). |
| Fraction of Target Data in all data | `target_bases / total_bases × 100%`. |
| Fraction of Target Data in mapped data | `target_bases / mapped_bases × 100%`. |
| Len of region | Total length of merged target regions (base pairs). |
| Average depth | `target_bases / region_length`. Mean sequencing depth over target regions. |
| Average depth(rmdup) | Mean depth after duplicate removal, an estimate of unique molecule coverage. |
| Coverage (>0.2*(Average depth)x) | Fraction of positions with depth > 20% of the mean depth. Indicates uniformity. Low values suggest uneven coverage. |
| Coverage (>0.5*(Average depth)x) | Fraction of positions with depth > 50% of the mean depth. Another uniformity metric. |
| Coverage (>0x) | Fraction of target positions with at least 1 read covering them. Also called "breadth of coverage". |
| Coverage (>=4x) | Positions with depth ≥ 4, a commonly used threshold for variant detection sensitivity. |
| Coverage (>=10x) | Positions with depth ≥ 10. |
| Coverage (>=30x) | Positions with depth ≥ 30, a typical clinical WES threshold. |
| Coverage (>=100x) | Positions with depth ≥ 100. |

### `[Target] Coverage(rmdup)` Duplicate-Removed Coverage

These use the **rmdup depth** distribution, which excludes PCR duplicates,
secondary/supplementary alignments, and reads with MAPQ below `-q`. This represents
unique molecule coverage.

| Field | Description |
|-------|-------------|
| Coverage(rmdup) (>0x) | Fraction of positions covered by at least one unique molecule. |
| Coverage(rmdup) (>=4x) | Unique molecule coverage at ≥4x. |
| Coverage(rmdup) (>=10x) | Unique molecule coverage at ≥10x. |
| Coverage(rmdup) (>=30x) | Unique molecule coverage at ≥30x. |
| Coverage(rmdup) (>=100x) | Unique molecule coverage at ≥100x. |

### `[Target] Coverage (>R*Avg)` Ratio-Based Coverage

When `--depthratio` is set, additional lines appear for each ratio R showing the
fraction of positions with depth ≥ R × average depth. Requires `--depthratio`.
For example, `--depthratio 0.2,0.5` produces these additional fields.

| Field | Description |
|-------|-------------|
| Coverage (>0.20*Avg) | Fraction of positions with depth > 20% of mean depth. |
| Coverage (>0.50*Avg) | Fraction of positions with depth > 50% of mean depth. |

### `[Target]` Region-Level Statistics

These metrics are computed per **target region** (each merged BED interval), not per
base. They describe how many regions achieve certain mean depth thresholds.

| Field | Description |
|-------|-------------|
| Target Region Count | Total number of merged target regions (e.g., number of exons). |
| Region covered > 0x | Number of regions with mean depth > 0. |
| Fraction Region covered > 0x | `regions_covered / total_regions × 100%`. |
| Fraction Region covered >= 4x | Regions with mean depth ≥ 4. |
| Fraction Region covered >= 10x | Regions with mean depth ≥ 10. |
| Fraction Region covered >= 30x | Regions with mean depth ≥ 30. |
| Fraction Region covered >= 100x | Regions with mean depth ≥ 100. |

### `[flank]` Flank Region Statistics

Flank regions are the areas adjacent to target regions (extended by `--flank` bases
on each side). These statistics exclude the target regions themselves.

| Field | Description |
|-------|-------------|
| flank size | The flank extension size in bp (set by `-f`, default 200). |
| Len of region (not include target region) | Total length of flank-only regions. |
| Average depth | Mean depth in flank regions. |
| flank Reads | Number of reads overlapping flank regions. A read may be counted in both target and flank if it spans the boundary. |
| Fraction of flank Reads in all reads | `n_flk / n_reads × 100%`. |
| Fraction of flank Reads in mapped reads | `n_flk / n_mapped × 100%`. |
| flank Data(Mb) | Total aligned bases in flank regions. |
| Fraction of flank Data in all data | `flank_bases / total_bases × 100%`. |
| Fraction of flank Data in mapped data | `flank_bases / mapped_bases × 100%`. |
| Coverage (>0x) | Fraction of flank bases with depth > 0. |
| Coverage (>=4x) | Flank bases with depth ≥ 4. |
| Coverage (>=10x) | Flank bases with depth ≥ 10. |
| Coverage (>=30x) | Flank bases with depth ≥ 30. |
| Coverage (>=100x) | Flank bases with depth ≥ 100. |

## Other Output Files

### chromosomes.report

TSV table (or CSV with `-F csv`) with one row per chromosome.

| Column | Description |
|--------|-------------|
| Chromosome | Chromosome/contig name. |
| DATA(%) | Percentage of total target data contributed by this chromosome. |
| Avg depth | Mean target depth on this chromosome. |
| Median | Median target depth on this chromosome. |
| Coverage% | Breadth of coverage (>0x) on this chromosome. |
| Cov 4x/10x/30x/100x % | Coverage at each threshold level. |

### depth_distribution.plot

Five-column data file for plotting cumulative depth distributions.

| Column | Description |
|--------|-------------|
| depth | Depth value (0, 1, 2, ...). |
| count | Number of target positions with exactly this depth. |
| fraction | `count / total_positions`. |
| cumulative_count | Remaining positions with depth > current value. |
| cumulative_fraction | `cumulative_count / total_positions`. |

### insertsize.plot

Same five-column format as depth_distribution.plot, but for inferred insert sizes
(0 to `--isize` limit).

### depth.tsv.gz

BGZF-compressed, tabix-indexable file with per-base depth values.

| Column | Description |
|--------|-------------|
| Chromosome | Chromosome/contig name. |
| Position | 1-based genomic coordinate. |
| Raw depth | Unfiltered depth (all reads). |
| Rmdup depth | Depth after duplicate removal and MAPQ filter. |
| Cover depth | Coverage depth (raw plus CDEL coverage). Used for coverage calculations. |

### region.tsv.gz

BGZF-compressed, tabix-indexable file with per-region summary.

| Column | Description |
|--------|-------------|
| Chromosome | Chromosome/contig name. |
| Start | 0-based region start. |
| Stop | 1-based region end. |
| Avg depth | Mean coverage depth for this region. |
| Median | Median coverage depth for this region. |
| Coverage | Breadth of coverage (>0x) for this region. |
| Coverage(FIX) | Alternative coverage using covdep array (accounts for deletions). |

### uncover.bed

BED file listing regions where coverage depth falls below the `--uncover` threshold
(default 5). These are candidate regions for probe design improvement.

## JSON Output Structure

When using `-F json`, a single `report.json` is produced.

```json
{
  "program": "bamdst",
  "version": "1.2.0",
  "input_files": ["in1.bam"],
  "reads": { "total": ..., "mapped": ..., ... },
  "insert_size": { "average": ..., "median": ... },
  "target": {
    "reads": ...,
    "avg_depth": ...,
    "coverage": { "gt_0x": ..., "ge_4x": ..., ... },
    "coverage_rmdup": { "gt_0x": ..., "ge_4x": ..., ... },
    "coverage_ratio": { "gt_0.20_avg": ..., ... }
  },
  "flank": { "size": ..., "avg_depth": ..., "coverage": { ... } },
  "chromosomes": [ { "name": "chr1", "avg_depth": ..., ... } ],
  "depth_distribution": [ { "depth": 0, "count": ..., "fraction": ... } ],
  "insert_size_distribution": [ { "size": 0, "count": ..., "fraction": ... } ]
}
```

## Notes

- The output directory specified by `-o` must already exist and be writable. bamdst
  will not create it. Create it beforehand with `mkdir -p /path/to/outdir`.
- BAM files **must be sorted by coordinate**. bamdst checks this at runtime and aborts
  if reads appear out of order.
- The probe BED file is 0-based by default. Use `-1` if your BED is 1-based.
- When processing multiple BAM files, they must share the same reference genome
  and chromosome ordering (the header from the first file is used).
- `depth.tsv.gz` and `region.tsv.gz` are compressed with bgzip and can be indexed
  with `tabix` for fast random access.
