#!/bin/bash

# Script arguments:
# 1) Path to the directory containing the input FASTQ/FASTQ.GZ files.
# 2) Name of the text file containing the list of paired-end read files.
#    R1 must be listed in line N and R2 in line N+1.
# 3) Path to the output directory.
#
# Example:
# qsub -R y -l h_rt=23:59:59 -t 1-886:2 preProcess_list_v2.sh \
#   /scratch/btaboada/Data/Ojeda \
#   FileNames.txt \
#   /scratch/btaboada/Data/Ojeda/Results

#$ -N preProcess_list
#$ -j y
#$ -o /scratch/btaboada/Programs
#$ -pe thread 4 -l h_vmem=4G

source "$HOME/.bashrc"
module load programs/fastp-0.20.0

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
  nSe1=${FILE1%%.*}
  nSe2=${FILE2%%.*}

  echo "$nSe1 $nSe2"

  # Build output file name.
  TEMP=''
  R1='R1_'
  FLT=${nSe1/$R1/$TEMP}

  echo ""
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

  echo ""
  echo "2. STARTING FASTP PROCESS"

  fastp \
    -i "$FILE1" \
    -I "$FILE2" \
    -o "$3/${nSe1}_Q.fastq.gz" \
    -O "$3/${nSe2}_Q.fastq.gz" \
    -V \
    -w "$NSLOTS" \
    -q 20 \
    --cut_tail \
    --cut_mean_quality 10 \
    -x --poly_x_min_len=40 \
    -g --poly_g_min_len=10 \
    -l 40 \
    -n 15 \
    -y -Y 40 \
    && echo "fastp completed successfully" \
    || { echo "fastp failed"; exit 1; }

  echo "END"
  echo ""
  echo "3. FINAL CLEAN READS:"

  if [ "$tipo" == 'fastq' ]; then
    echo "Final number of reads in ${nSe1}_Q.fastq.gz ${nSe2}_Q.fastq.gz:"
    numLin=$(cat "$3/${nSe1}_Q.fastq.gz" "$3/${nSe2}_Q.fastq.gz" | wc -l)
    expr "$numLin" / 4
  elif [ "$tipo" == 'gz' ]; then
    echo "Final number of reads in ${nSe1}_Q.fastq.gz ${nSe2}_Q.fastq.gz:"
    numLin=$(zcat "$3/${nSe1}_Q.fastq.gz" "$3/${nSe2}_Q.fastq.gz" | wc -l)
    expr "$numLin" / 4
  else
    echo "The provided file is not fastq or fastq.gz"
  fi
else
  echo "ERROR: three input parameters are required."
  exit 1
fi
