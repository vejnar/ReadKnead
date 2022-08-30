# <img src="https://raw.githubusercontent.com/vejnar/ReadKnead/main/img/logo.svg" alt="ReadKnead" width="45%" />

[![MPLv2](https://img.shields.io/aur/license/ReadKnead?color=1793d1&style=for-the-badge)](https://mozilla.org/MPL/2.0/)

<p align="center" style="text-align:center;"><i>Knead your sequencing reads before baking</i></p>

ReadKnead **clips**, **trims**, **demultiplexes**, **filters** (e.g. by length), **selects** (e.g. randomly) and **renames** reads from FASTQ files.
* Choice of algorithm for trimming: fast and accurate [bit-masked k-difference matching](https://git.sr.ht/~vejnar/bktrim), Needleman–Wunsch, search or match.
* Demultiplexing using internal barcodes (user-defined positions in the reads)

For testing read preparation pipelines and quality control, ReadKnead:
* Plots read-length barplot
* Plots quality scores boxplot
* Provide detailed stats for each operation (clip, trim etc)

## Why?

* **Integrated pipeline for read preparation** When individual tools are used for read preparation, demultiplexing can't depend on successful trimming for example. With **ReadKnead**, users define read preparation pipeline with multiple dependent steps.
* **Speed** ReadKnead is written in Go and parallelized.

## Download

See [refs](https://git.sr.ht/~vejnar/ReadKnead/refs) page for tarball and executable.

## Running tests

```bash
go test -v ./cmd/readknead/...
```

## Examples

Input files for following examples can be found in the `cmd/readknead/testdata` directory. Please note that FASTQ files found in that directory are uncompressed, while following examples demonstrate how to use compressed FASTQ files using multiple compressor (Gzip, Zstandard and LZ4). We suggest to run the tests (see above) to see the output of these examples.

**Highly recommended to use the `-verbose_level 20` argument to test pipelines.**

### Single-end clipping and trimming

First, define a pipeline in `clip_trim.json` file for single-end reads that will:
1. Clip 10 nucleotides from the 5' end of reads
2. Trim the 3' end of reads using the [bit-masked k-difference matching](https://git.sr.ht/~vejnar/bktrim) keeping all reads including untrimmed reads (`no_trim`).
3. Remove reads shorter than 20 nucleotides

```json
[{"name": "clip",
  "end": 5,
  "length": 10},
 {"name": "trim",
  "end": 3,
  "algo": "bktrim",
  "epsilon": 0.15,
  "epsilon_indel": 0.1,
  "keep": ["no_trim",
           "trim_exact",
           "trim_align"],
  "sequence":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"},
 {"name": "length",
  "min_length": 20}]
```

Then run the pipeline defined above in the JSON file on input file `sample1_R1.fastq.gz`:
```bash
readknead -fq_fnames_r1 "sample1_R1.fastq.gz" \
          -fq_command_in "zcat" \
          -fq_path_out "output" \
          -fq_fname_out_r1 "sample1_R1.fastq.lz4" \
          -fq_command_out "lz4,-" \
          -ops_r1_path "clip_trim.json" \
          -report_path "output/preparing_report.json"
```

### Paired-end clipping, trimming, filtering and renaming

First, define a pipeline in `paired_end_trim.json` file for paired-end reads that will:
1. Rename reads to shorten their names using the pattern `sample2.##` where `##` will be replaced by the read number. In case read names contain a barcode (prefixed with `#`), it will be kept in the renamed reads.
2. Trim the reads on both ends using the [bit-masked k-difference matching](https://git.sr.ht/~vejnar/bktrim) keeping all reads including untrimmed reads (`no_trim`).
3. Remove reads shorter than 20 nucleotides

```json
[{"name": "rename",
  "base36": true,
  "keep_barcode": true,
  "new_name": "sample2."},
 {"name": "trim",
  "algo": "bktrim_paired",
  "epsilon": 0.15,
  "epsilon_indel": 0.1,
  "keep": ["no_trim",
           "trim_exact",
           "trim_align"],
  "sequence":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
  "sequence_paired":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"},
 {"name": "length",
  "min_length": 20}]
```

To generate quality-control stats, quality score encoding (Illumina 1.8+ Phred+33) argument is added. Four workers are used. Run the pipeline defined above in the JSON file on input files `sample2_R1.fastq.gz` and `sample2_R2.fastq.gz`:
```bash
readknead -fq_fnames_r1 "sample2_R1.fastq.zst" \
          -fq_fnames_r2 "sample2_R2.fastq.zst" \
          -fq_command_in "zstdcat" \
          -fq_path_out "output" \
          -fq_fname_out_r1 "sample2_R1.fastq.zst" \
          -fq_fname_out_r2 "sample2_R2.fastq.zst" \
          -fq_command_out "zstd,-,-fo" \
          -ops_r1_path "paired_end_trim.json" \
          -report_path "output/preparing_report.json" \
          -label "WT replicate 2" \
          -ascii_min "33" \
          -stats_in_path "output/stats_in" \
          -stats_out_path "output/stats_out" \
          -num_worker 4
```

The input FASTQ files are compressed using [Zstandard](https://github.com/facebook/zstd).

### Demultiplexing

First, define a pipeline in `demultiplex.json` file for paired-end reads that will:
1. Trim/Filter the reads on 5' end with a match algorithm with one of the sequences listed in `sequences`. Only perfect match and >70% match reads are kept. Trimmed sequences are not added to the read names (`"add_trimmed": false`) and matched sequences are not actually trimmed from the reads (`"apply_trim_seq": false`): this step is only a filter.
2. Clip 10 nucleotides from the 5' end of reads and add the clipped sequence to the read names
3. Reads are separated (i.e. demultiplexed) matching one of the barcodes listed in `barcodes` (at position 10 from the 5' end with maximum of 1 mismatch). One FASTQ file will be created per barcode (replacing the sequence of the barcode in the output FASTQ file in `[DPX]` in the `-fq_fname_out_r1` and `-fq_fname_out_r2` arguments).
4. Clip 19 nucleotides from the 5' end of reads and don't add the clipped sequence to the read names
5. Remove reads shorter than 20 nucleotides

```json
[{"name": "trim",
  "end": 5,
  "algo": "match",
  "position": 6,
  "min_score": 0.7,
  "add_trimmed": false,
  "keep": ["trim_exact",
           "trim_align"],
  "sequences": ["CATTGCTTATGG",
                "GTACGGGACTTA"],
  "apply_trim_seq": false},
 {"name": "clip",
  "end": 5,
  "length": 10,
  "add_clipped": true},
 {"name": "demultiplex",
  "end": 5,
  "max_mismatch": 1,
  "barcodes": ["GAGTA",
               "CTGAG"]},
 {"name": "clip",
  "end": 5,
  "length": 19,
  "add_clipped": false},
 {"name": "length",
  "min_length": 20}]
```

Then run the pipeline defined above in the JSON file on input files `sample2_R1.fastq.zst` and `sample2_R2.fastq.zst`:
```bash
readknead -fq_fnames_r1 "sample2_R1.fastq.zst" \
          -fq_fnames_r2 "sample2_R2.fastq.zst" \
          -fq_command_in "zstdcat" \
          -fq_path_out "output" \
          -fq_fname_out_r1 "sample2_[DPX]_R1.fastq.zst" \
          -fq_fname_out_r2 "sample2_[DPX]_R2.fastq.zst" \
          -fq_command_out "zstdcat" \
          -ops_r1_path "demultiplex.json" \
          -report_path "output/preparing_report.json" \
          -label "WT replicate 2" \
          -ascii_min "33" \
          -stats_in_path "output/stats_in" \
          -stats_out_path "output/stats_out" \
          -num_worker 4
```

## Command-line arguments

* Input
    * `-fq_fnames_r1` Path to read 1 FASTQ files (comma separated)
    * `-fq_fnames_r2` Path to read 2 FASTQ files (comma separated)
    * `-fq_command_in` Command line to execute for opening each input file (comma separated)
    * `-buf_size` Buffer IO size (default 41943040)
* Output
    * `-fq_path_out`  Path to output FASTQ files
    * `-fq_fname_out_r1` Output read 1 FASTQ file
    * `-fq_fname_out_r2` Output read 2 FASTQ file
    * `-fq_command_out` Command line to execute for opening each output file (comma separated)
* Pipeline
    * `-ops_r1` Operation(s) for read1
    * `-ops_r1_path` Path to operation(s) for read1
    * `-ops_r2` Operation(s) for read2
    * `-ops_r2_path` Path to operation(s) for read2
* Pipeline report
    * `-report_path` Write report to path (stdout with -)
* Plotting 
    * Activate plotting for input and/or output:
        * `-stats_in_path` Path to statistics of input FASTQ files
        * `-stats_out_path` Path to statistics of output FASTQ files
    * Configure plotting and read characteristics:
        * `-label` Label
        * `-ascii_min` ASCII coding of the minimum quality score (Ex: 33 for Phred+33) (default: 33)
        * `-max_quality` Highest nucleotide quality (default: 43)
        * `-max_read_length` Maximum read length. No precision required (default: 1000)
* Other options
    * `-num_worker` Number of worker(s) (default 1)
    * `-verbose` Verbose
    * `-verbose_level` Verbose level. This option is useful for testing pipelines.
    * `-version` Print version and quit

## Operations

| Operation   | Parameter        | Type      | Default                 |                                                                               |
|-------------|------------------|-----------|-------------------------|-------------------------------------------------------------------------------|
| clip        | length           | integer   |                         | Number of nucleotide to clip                                                  |
|             | end              | integer   |                         | End of read to clip: 5 or 3                                                   |
|             | add_clipped      | boolean   | false                   | Copy clipped nucleotide to read name (#-prefixed)                             |
|             | add_separator    | boolean   | true                    | Add prefix (#) before clipped sequence in read name                           |
| demultiplex | barcodes         | []strings |                         | List of barcode sequences                                                     |
|             | end              | integer   |                         | End of read to clip: 5 or 3                                                   |
|             | barcode_idx      | integer   |                         | Index (first: 0) of #-prefixed sequence (barcode or UMI) in read name         |
|             | max_mismatch     | integer   | 0                       | Maximum number of mismatch between read and barcode                           |
|             | length_ligand    | integer   | 0                       | Clip if barcode found                                                         |
| length      | min_length       | integer   | -1                      | Minimum read length                                                           |
|             | max_length       | integer   | -1                      | Maximum read length                                                           |
| random      | probability      | float     | 1.                      | Probability to keep read (between 0 and 1)                                    |
| rename      | new_name         | string    |                         | New read name                                                                 |
|             | base36           | boolean   | false                   | Convert read number to shorter base36                                         |
|             | keep_barcode     | boolean   | false                   | Keep #-prefixed sequences                                                     |
|             | merge_barcode    | boolean   | false                   | Merge #-prefixed sequences to one                                             |
|             | all_reads        | boolean   | true                    | Rename all reads                                                              |
| trim        | sequence         | string    |                         | Sequence to trim (for pair-end reads: downstream sequence)                    |
|             | sequences        | []strings |                         | Use for multiple-sequence trimming                                            |
|             | sequence_paired  | string    |                         | Upstream sequence to trim for paired-end reads                                |
|             | sequences_paired | []strings |                         | Use for multiple-sequence paired-end reads trimming                           |
|             | add_trimmed      | boolean   | false                   | Copy trimmed nucleotide to read name (#-prefixed)                             |
|             | add_trimmed_ref  | boolean   | false                   | Copy reference trimming sequence to read name (#-prefixed)                    |
|             | add_separator    | boolean   | true                    | Add prefix (#) before trimmed sequence in read name                           |
|             | algo             | string    | bktrim or bktrim_paired | Algorithms: *align*, *bktrim*, *bktrim_paired*, *search* or *match*           |
|             | end              | integer   |                         | End of read to trim : 5 or 3 (only for *align*, *bktrim* and *search* `algo`) |
|             | min_sequence     | integer   | 0                       | Length of perfect match (starting at trimming position) in trimming alignment |
|             | min_score        | float     | 0.8                     | Minimum alignment score (only for *align*, *search* and *match* `algo`)       |
|             | position         | integer   | 0                       | Position in reads to match trimming sequence (only for *match* `algo`)        |
|             | epsilon          | float     | 0.1                     | Maximum mismatch ratio (only for *bktrim* and *bktrim_paired* `algo`)         |
|             | epsilon_indel    | float     | 0.03                    | Maximum indel ratio (only for *bktrim* and *bktrim_paired* `algo`)            |
|             | min_overlap      | integer   | 3                       | Minimum overlap length (only for *bktrim* and *bktrim_paired* `algo`)         |

## License

*ReadKnead* is distributed under the Mozilla Public License Version 2.0 (see /LICENSE).

Copyright (C) 2017-2022 Charles E. Vejnar
