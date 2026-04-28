# Influenza Genomes Pipeline

Bash-based pipeline for processing Illumina paired-end reads and generating segment-level consensus genomes for Influenza A and B viruses.

The workflow is organized in two main stages. First, raw reads are cleaned with `fastp`. Second, cleaned reads are mapped against multiple influenza reference genomes. The best reference is selected automatically using mapping support on segment 4 (HA) and segment 6 (NA), followed by consensus and variant calling with `iVar`.

## Included scripts

- `scripts/preProcess_list_v2.sh`: batch read preprocessing with `fastp`; suitable for SGE-based clusters.
- `scripts/ConsensoInfluenza_v3.sh`: minimal multi-reference mapping, best-reference selection using S4/S6 (HA/NA), and consensus generation with `iVar`.

## Requirements

- Linux environment; tested on SGE clusters.
- Required tools: `fastp`, `bowtie2`, `samtools`, `bcftools`, `bamtools`, and `ivar`.
- SGE variables used by the scripts: `$SGE_TASK_ID` and `$NSLOTS`.

## Input structure

- Paired FASTQ or FASTQ.GZ files in the same directory, using an `R1_` / `R2_` naming convention, for example:

```text
SampleX_R1_001.fastq.gz
SampleX_R2_001.fastq.gz
```

- A `FileNames.txt` file listing paired reads in consecutive lines. For each sample, R1 must be listed in line N and R2 in line N+1. This pairing is handled through `$SGE_TASK_ID`.

An example file is provided in `examples/FileNames.txt`.

## References

Place all reference genomes to be evaluated in `refs/`, and index each one with Bowtie2:

```bash
bowtie2-build A_H1N1.fasta A_H1N1
bowtie2-build A_H3N2.fasta A_H3N2
bowtie2-build B_Victoria.fasta B_Victoria
bowtie2-build B_Yamagata.fasta B_Yamagata
```

For each reference, include:

- `{ref}.fasta` and its FASTA index (`{ref}.fasta.fai`).
- `{ref}.bed` with segment coordinates. Segment names must be consistent, especially for `S4_HA` and `S6_NA`.
- `{ref}.gff3` files per segment, if required by the workflow.

Best-reference selection is based on the number of alignments supporting S4/S6 in each mapping. Make sure that segment names in FASTA, BED, and GFF3 files match the names expected by the script. By default, the expected segment names are `S4_HA` and `S6_NA`.

## Stage 1: read preprocessing

Example command; adjust task ranges, queues, and paths according to your cluster:

```bash
qsub -R y -l h_rt=23:59:59 -pe thread 4 -t 1-100 \
  scripts/preProcess_list_v2.sh \
  /path/to/FASTQ \
  FileNames.txt \
  /path/to/Results
```

Output: cleaned paired files named `{Sample}_Q.fastq.gz` in the `Results` directory, along with logs containing pre- and post-filtering read counts.

## Stage 2: consensus and variant calling

Run the consensus script using the directory containing cleaned FASTQ files and the same `FileNames.txt` file:

```bash
qsub -R y -l h_rt=23:59:59 -pe thread 4 -t 1-100 \
  scripts/ConsensoInfluenza_v3.sh \
  /path/to/Results \
  FileNames.txt \
  /absolute/path/to/refs
```

Per-sample workflow:

1. Map reads against each reference using `bowtie2 --very-sensitive-local`.
2. Count alignments supporting `S4_HA` and `S6_NA`.
3. Select the reference with the highest mapping support.
4. Generate sorted and indexed BAM files, VCF files, segment-level consensus FASTA files, and variant tables.
5. Organize outputs into `FASTA/`, `SAM/`, `BAM/`, and `VIRUS/` directories.
6. Generate a concatenated FASTA file named `{sample}_All_ivar2.fasta`.

## Troubleshooting

- If no reference is selected, check that segment names match exactly across FASTA, BED, and GFF3 files.
- For very short reads (<60 nt), consider adjusting Bowtie2 `--score-min` parameters.
- By default, `ivar consensus` uses `-q 20`, `-m 5`, and `-t 0`, masking low-coverage positions with `N`. Adjust these thresholds according to your analysis criteria.
- Verify that `$NSLOTS` is defined by your SGE queue. Otherwise, manually set the number of threads in the scripts.

## Repository structure

```text
influenza-genomes-pipeline/
├─ scripts/
│  ├─ preProcess_list_v2.sh
│  └─ ConsensoInfluenza_v3.sh
├─ refs/
├─ examples/
│  └─ FileNames.txt
├─ docs/
│  └─ PROTOCOL.md
└─ README.md
```

## Citation

If you use this pipeline, please cite the associated publication or repository version when available.

## Acknowledgements

This work was supported by Universidad Nacional Autónoma de México (UNAM) through PAPIIT-DGAPA-IN230523 awarded to Blanca Taboada, and by Secretaría de Educación, Ciencia, Tecnología e Innovación de la Ciudad de México (SECTEI) through project SECTEI/138/2024 awarded to Selene Zárate.
