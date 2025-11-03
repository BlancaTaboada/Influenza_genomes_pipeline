#!/bin/bash

#Parameters:
## 1) Path of the file (txt) with the names of the libraries to analysis. 
	##It must be in the same place of the libraries. Example: /scratch/btaboada/Data/Flavi
# 2) Name of the file with the names of the libraries to analyze. Example: ISA_Todos.txt
# 3) Path of the output results. /scratch/btaboada/Data/Flavi
## Example: qsub -R y -l h_rt=23:59:59 -t 1-886:2 preProcess_list_v2.sh /scratch/btaboada/Data/Ojeda Filenames.txt /scratch/btaboada/Data/Ojeda

# Process name
#$ -N preProcess_list

# Salida estandar y de errror juntas en la salida estandar y las coloca en la carpeta de Programas
#$ -j y
#$ -o /scratch/btaboada/Programs
#$ -pe thread 4 
-l h_vmem=4G

source $HOME/.bashrc
module load programs/fastp-0.20.0

if [ $# -ge 3 ]; then
	cd $1
	FILE1=$(cat $2 | head -n $SGE_TASK_ID| tail -1)
	echo $FILE1
	num=$(($SGE_TASK_ID + 1))
	FILE2=$(cat $2 | head -n $num| tail -1)
	echo $FILE2

	# Obtener extensión del archivo
	nSe=${FILE1%.*}
	tipo=${FILE1#$nSe.}
	nSe1=${FILE1%%.*}	
	nSe2=${FILE2%%.*}	

	echo $nSE1 $nSE2

	# Nombre de archivos salida
	TEMP=''
	R1='R1_'
	FLT=${nSe1/$R1/$TEMP}

	echo " "
	echo "1. INITIAL READS: "
		# Revisar extensión y contar lecturas
		if [ $tipo == 'fastq' ]; then
			echo "Numero de lecturas iniciales en $FILE1 $FILE2 :"
			numLin=$(cat $FILE1 $FILE2 | wc -l)
			expr $numLin / 4
		elif [ $tipo == 'gz' ]; then
			echo "Numero de lecturas iniciales en $FILE1 $FILE2 :"
			numLin=$(zcat $FILE1 $FILE2 | wc -l)
			expr $numLin / 4
		else
			echo "El archivo dado NO es fastq o fastq.gz"
		fi

	echo ""
	echo "2. STARTING FASTP PROCESS"
		fastp -i $FILE1 -I $FILE2 -o $3/$nSe1"_Q.fastq.gz" -O $3/$nSe2"_Q.fastq.gz" -V -w $NSLOTS -q 20 --cut_tail --cut_mean_quality 10 -x --poly_x_min_len=40 -g --poly_g_min_len=10 -l 40 -n 15 -y -Y 40
		echo "fastp -i $FILE1 -I $FILE2 -o $3/$nSe1"_Q.fastq.gz" -O $3/$nSe2"_Q.fastq.gz" -V -w $NSLOTS -q 20 --cut_tail --cut_mean_quality 10 -x --poly_x_min_len=40 -g --poly_g_min_len=10 -l 40 -n 15 -y -Y 50 ... O.K" || 
		( echo "fastp -i $FILE1 -I $FILE2 -o $3/$nSe1"_Q.fastq.gz" -O $3/$nSe2"_Q.fastq.gz" -V -w $NSLOTS -q 20 --cut_tail --cut_mean_quality 10-x --poly_x_min_len=40 -g --poly_g_min_len=10 -l 40 -n 15 -y -Y 50 ... ERROR" ; exit 1 ) 
		echo "END"

	echo " "
	echo "3. FIAL CLEAN READS: "
		if [ $tipo == 'fastq' ]; then
			echo "Numero de lecturas finales en $nSe1"_Q.fastq.gz" $nSe2"_Q.fastq.gz" :"
			numLin=$(cat $3/$nSe1"_Q.fastq.gz" $3/$nSe2"_Q.fastq.gz" | wc -l)
			expr $numLin / 4
		elif [ $tipo == 'gz' ]; then
			echo "Numero de lecturas finales en $nSe1"_Q.fastq.gz" $nSe2"_Q.fastq.gz" :"
			numLin=$(zcat $3/$nSe1"_Q.fastq.gz" $3/$nSe2"_Q.fastq.gz" | wc -l)
			expr $numLin / 4
		else
			echo "El archivo dado NO es fastq o fastq.gz"
		fi

fi
