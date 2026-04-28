#!/bin/bash

# Before running this script, all influenza reference genomes must be indexed with Bowtie2.
# Reference FASTA, FAI, BED, and GFF3 files must be located in the same reference directory.
# A BED and GFF3 file must exist for each segment of each reference genome, when required.
# The output directory should contain or allow creation of the SAM, BAM, FASTA, and VIRUS folders.
#
# Script arguments:
# 1) Path to the directory containing the input FASTQ/FASTQ.GZ files.
# 2) Name of the text file containing the list of paired-end read files.
#    R1 must be listed in line N and R2 in line N+1.
# 3) Path to the directory containing reference FASTA files and Bowtie2 indexes.
#
# Example:
# qsub -R y -l h_rt=23:59:59 -t 1-10:2 ConsensoInfluenza_v3.sh \
#   /scratch/btaboada/Influenza/Results \
#   FileNames.txt \
#   /scratch/btaboada/DBs/Influenza

source "$HOME/.bashrc"
module load programs/bowtie2-2.5.0
module load programs/samtools-1.10
module load programs/ivar-1.3.1

#$ -N InfluenzaConse
#$ -j y
#$ -o /scratch/btaboada/Programs
#$ -pe thread 4

if [ $# -ge 3 ]; then
  cd "$1" || exit 1

  FILE1=$(cat "$2" | head -n "$SGE_TASK_ID" | tail -1)
  echo "$FILE1"

  num=$((SGE_TASK_ID + 1))
  FILE2=$(cat "$2" | head -n "$num" | tail -1)
  echo "$FILE2"

  # Get file extension.
  nSe=${FILE1%.*}
  tipo=${FILE1#$nSe.}
  nSe=${FILE1%%.*}

  # Build output file name.
  TEMP=''
  R1='R1_'
  FLT=${nSe/$R1/$TEMP}

  echo "1. INITIAL READS:"

  # Check file extension and count reads.
  if [ "$tipo" == 'fastq' ]; then
    echo "Initial number of reads in $FILE1 $FILE2:"
    numLin=$(cat "$FILE1" "$FILE2" | wc -l)
    expr "$numLin" / 4
  elif [ "$tipo" == 'gz' ]; then
    echo "Initial number of reads in $FILE1 $FILE2:"
    numLin=$(zcat "$FILE1" "$FILE2" | wc -l)
    expr "$numLin" / 4
  else
    echo "The provided file is not fastq or fastq.gz"
  fi

  cd "$3" || exit 1
  Refs=$(ls *.fasta)
  cd "$1" || exit 1

  for i in $Refs; do
    nameRef=${i%.*}
    nSe=${FLT}_${nameRef}
    bwRef=${3}${nameRef}

    echo "2. MAPPING READS TO REFERENCE WITH BOWTIE..."

    bowtie2 \
      --very-sensitive-local \
      --soft-clipped-unmapped-tlen \
      -t \
      --sam-no-qname-trunc \
      -X 3000 \
      -x "$bwRef" \
      -1 "$FILE1" \
      -2 "$FILE2" \
      -S "${nSe}.sam" \
      --threads "$NSLOTS" \
      && echo "bowtie2 mapping completed successfully for $nameRef" \
      || { echo "bowtie2 mapping failed for $nameRef"; exit 1; }

    awk 'BEGIN {FS = OFS = "\t"}!/^@/{if ($2 ~ /77/){print $1,$10,$11}}' "${nSe}.sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' > "${nSe}_R1_noM.fastq"
    awk 'BEGIN {FS = OFS = "\t"}!/^@/{if ($2 ~ /141/){print $1,$10,$11}}' "${nSe}.sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' > "${nSe}_R2_noM.fastq"
    awk 'BEGIN {FS = OFS = "\t"}!/^@/{if (($2 !~ /141/ && $2 !~ /77/) && $1 ~ /1:N:0/ ){print $1,$10,$11}}' "${nSe}.sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' > "${nSe}_R1_M.fastq"
    awk 'BEGIN {FS = OFS = "\t"}!/^@/{if (($2 !~ /141/ && $2 !~ /77/) && $1 ~ /2:N:0/ ){print $1,$10,$11}}' "${nSe}.sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' > "${nSe}_R2_M.fastq"

    echo ""
  done

  Sams=$(ls "${FLT}"*.sam)
  MapsMax=0

  for j in $Sams; do
    echo "$j"
    Maps=$(grep -c -E 'S6_NA|S4_HA' "$j")
    echo "$Maps"

    if [ "$Maps" -gt "$MapsMax" ]; then
      MapsMax=$Maps
      nameSeq=${j%.sam}
      echo "$nameSeq"
      BowRef=${nameSeq#*Dp_}
      echo "$BowRef"
    fi
  done

  nameBase=${nameSeq%_Dp*}
  echo "Maximum mapping support: $nameSeq $MapsMax $BowRef"

  find . -type f -name "*$nameBase*" -not -name "*Dp.fastq.gz" -not -name "*$BowRef*" -exec rm {} \;

  echo ""
  echo "3. CONVERTING SAM TO BAM"

  samtools view \
    -bT "$3/${BowRef}.fasta" \
    "${nameSeq}.sam" \
    -o "${nameSeq}.bam" \
    --threads 2 \
    && echo "samtools view completed successfully" \
    || { echo "samtools view failed"; exit 1; }

  echo ""
  echo "4. SORTING BAM"

  samtools sort "${nameSeq}.bam" -o "${nameSeq}.srt.bam" --threads 2 \
    && echo "samtools sort completed successfully" \
    || { echo "samtools sort failed"; exit 1; }

  rm "${nameSeq}.bam"

  samtools index "${nameSeq}.srt.bam" \
    && echo "samtools index completed successfully" \
    || { echo "samtools index failed"; exit 1; }

  echo ""
  echo "5. CALLING VARIANTS AND CLEANING CONSENSUS FILES"

  segmen=$(cut -f1 "$3/${BowRef}.bed")

  for j in $segmen; do
    echo "$j"

    samtools mpileup -aa -A -B -d 0 -Q 10 \
      -r "$j" \
      -l "$3/${BowRef}.bed" \
      -f "$3/${BowRef}.fasta" \
      "${nameSeq}.srt.bam" \
      | ivar consensus -q 20 -t 0 -m 5 -n N -p "${nameSeq}_${j}.ivar" \
      && echo "Consensus generated successfully for $j" \
      || { echo "Consensus generation failed for $j"; exit 1; }

    samtools mpileup -A -B -d 0 -Q 10 \
      -r "$j" \
      -l "$3/${j}.bed" \
      -f "$3/${BowRef}.fasta" \
      "${nameSeq}.srt.bam" \
      | ivar variants \
        -p "${nameSeq}_${j}.ivar" \
        -q 20 \
        -t 0.0 \
        -r "$3/${BowRef}.fasta" \
        -g "$3/${j}.gff3" \
      && echo "Variants called successfully for $j" \
      || { echo "Variant calling failed for $j"; exit 1; }

    perl /scratch/btaboada/Programs/cleanConsensusSeg_v1.pl \
      "${nameSeq}_${j}.ivar.fa" \
      "${nameSeq}_${j}.ivar.tsv" \
      "${nameSeq}_${j}.ivar2.fasta"

    sed -i 's/Consensus_//g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/\.ivar_threshold_0_quality_20//g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/Q_Dp_//g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/A_H3N2_A_H3N2/A_H3N2/g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/A_H1N1_A_H1N1/A_H1N1/g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/B_Yama_B_Yama/B_Yama/g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/B_Vict_B_Vict/B_Vict/g' "${nameSeq}_${j}.ivar2.fasta"
    sed -i 's/_S[0-9]\+//' "${nameSeq}_${j}.ivar2.fasta"
  done

  rm $(ls "${nameSeq}"*.fa)

  zipF=$(ls "${nameSeq}"*.fastq "${nameSeq}"*.sam "${nameSeq}"*qual.txt "${nameSeq}"*.ivar2.fasta "${nameSeq}"*.ivar.tsv 2>/dev/null)
  for z in $zipF; do
    gzip -q "$z"
  done

  mkdir -p FASTA SAM BAM VIRUS

  zcat $(ls "${nameSeq}"*ivar2.fasta.gz) > "${nameSeq}_All_ivar2.fasta"

  mv $(ls "${nameSeq}"*.sam.gz) SAM
  mv $(ls "${nameSeq}"*.bam*) BAM
  mv $(ls "${nameSeq}"*.ivar2.fasta.gz) FASTA
  mv $(ls "${nameSeq}"*.ivar.tsv.gz) FASTA
  mv $(ls "${nameSeq}"*.ivar.vcf.gz 2>/dev/null) FASTA 2>/dev/null || true
  mv $(ls "${nameSeq}"*.ivar.qual.txt.gz 2>/dev/null) FASTA 2>/dev/null || true
  mv $(ls "${nameSeq}"*_M.fastq.gz) VIRUS

  cd FASTA || exit 1
  rename 'A_H3N2_A_H3N2' 'A_H3N2' *.gz 2>/dev/null || true
  rename 'A_H1N1_A_H1N1' 'A_H1N1' *.gz 2>/dev/null || true
  rename '_B_Vict_B_Vict' '_B_Vict' *.gz 2>/dev/null || true
  rename 'B_Yama_B_Yama' 'B_Yama_B' *.gz 2>/dev/null || true
else
  echo "ERROR: three input parameters are required."
  exit 1
fi
