#!/bin/bash
#Antes de correr el script debe estar los 4 genomas de referencia de influenza indexados para bowtie 
# y estar en la misma carpeta todos asi como sus correspondientes fasta, fai, bed y gff3. 
# Debe existir un bed y un gff3 por cada segemnto de cada uno de los genomas.
# Deben existir una carpeta SAM, BAM, FASTA en la carpeta de resultados.

# Argumentos del script
##Parameters:
## 1) Path del archivo txt con los nombres de las muestras analizar. /path/
## 2) Nombre del archivo txt con los nombres de las muestras analizar. ej. FileName.txt
## 3) Path donde estan los genomas fasta de referencia con la extension fasta. ej. Geneoma.fasta
##		Deben tener los fasta el mismo nombre que los indexados por bowtie. Va a comparar 
##		con todos los genoma que esten en la carpeta. En este caso 4.
## qsub -R y -l h_rt=23:59:59 -t 1-10:2 ConsensoInfluenza_v2.sh /scratch/btaboada/2022_CursoMeta/Result FilesName.txt /scratch/btaboada/DBs/Influenza
## FileName.txt: Nombre de archivos analizar.

source $HOME/.bashrc
module load programs/bowtie2-2.5.0
module load programs/samtools-1.10
module load programs/ivar-1.3.1
# Salida estandar y de errror juntas en la salida estandar y las coloca en la carpeta de Programas
#$ -N InfluenzaConse
#$ -j y
#$ -o /scratch/btaboada/Programs
#$ -pe thread 4

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
	nSe=${FILE1%%.*}	
	# Nombre de archivos salida
	TEMP=''
	R1='R1_'
	FLT=${nSe/$R1/$TEMP}

	echo "1. INITIAL READS: " 
	#Revisar extensión y contar lecturas
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
	cd $3
	Refs=$(ls *.fasta)
	cd $1
	for i in $Refs
	do
		nameRef=${i%.*} #Nombre base del genoma de refencia
		nSe=$FLT"_"$nameRef #echo "Nombre base de referencia $nameRef"
		bwRef=${3}${nameRef} #echo "Nombre base de resultado $nSe"
		
		echo "2. MAPPING READS TO REFERENCE WITH BOWTIE...." #G,30 cuando son de 75
		bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen -t --sam-no-qname-trunc -X 3000 -x $bwRef -1 $FILE1 -2 $FILE2 -S $nSe".sam" --threads $NSLOTS &&
		echo "bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen -t --sam-no-qname-trunc -X 3000 -x $bwRef -1 $FILE1 -2 $FILE2 -S $nSe".sam" --threads $NSLOTS " ||
		( echo "bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen -t --sam-no-qname-trunc -X 3000 -x $bwRef -1 $FILE1 -2 $FILE2 -S $nSe".sam" --threads $NSLOTS " ; exit 1 )
		awk 'BEGIN {FS = OFS = "\t"}!/^@/{if ($2 ~ /77/){print $1,$10,$11}}' $nSe".sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' >$nSe"_R1_noM.fastq"
		awk 'BEGIN {FS = OFS = "\t"}!/^@/{if ($2 ~ /141/){print $1,$10,$11}}' $nSe".sam" |awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' >$nSe"_R2_noM.fastq"
		awk 'BEGIN {FS = OFS = "\t"}!/^@/{if (($2 !~ /141/ && $2 !~ /77/) && $1~ /1:N:0/ ){print $1,$10,$11}}' $nSe".sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' >$nSe"_R1_M.fastq"
		awk 'BEGIN {FS = OFS = "\t"}!/^@/{if (($2 !~ /141/ && $2 !~ /77/) && $1~ /2:N:0/ ){print $1,$10,$11}}' $nSe".sam" | awk '{ print "@"$1" "$2"\n"$3"\n""+\n"$4}' >$nSe"_R2_M.fastq"
	echo ""
	done
	Sams=$(ls "$FLT"*.sam"")
	MapsMax=0
	for j in $Sams 
	do 
		echo $j 
		Maps=$(grep -c -E 'S6_NA|S4_HA' $j)
		echo $Maps
		if [ ${Maps} -gt ${MapsMax} ]; then
			 MapsMax=$(($Maps))
			 nameSeq=${j%.sam}
			 echo "$nameSeq"
			 BowRef=${nameSeq#*Dp_}
			 echo "$BowRef"
		fi
	done
	nameBase=${nameSeq%_Dp*}
	echo "El mapeo Maximo es $nameSeq $MapsMax $BowRef"
	find -type f -name "*$nameBase*" -not -name "*Dp.fastq.gz" -not -name "*$BowRef*" -exec rm {} \;
	echo ""
	
	echo "3. CONVERTIR SAM -> BAM "
	samtools view -bT $3"/"$BowRef".fasta" $nameSeq".sam" -o $nameSeq".bam" --threads 2 && echo "samtools view -bT $3"/"$BowRef".fasta" $nameSeq".sam" -o $nameSeq".bam" --threads 2 ... OK" || 
		( echo "samtools view -bT $3"/"$BowRef".fasta" $nameSeq".sam" -o $nameSeq".bam" --threads 2 ... ERROR" ; exit 1 ) 
	echo""
	
	if [ $? -eq 0 ] 
	then
		echo "4. ORDENAR BAM ..."
		samtools sort $nameSeq".bam" -o $nameSeq".srt.bam" --threads 2 && echo "samtools sort $nameSeq".bam" -o $nameSeq".srt.bam" --threads 2 ... OK" ||
			( echo "samtools sort $nameSeq".bam" -o $nameSeq".srt.bam" --threads 2 ... ERROR"; exit 1 )
		rm $nameSeq".bam"
		samtools index $nameSeq".srt.bam" && echo "samtools index $nameSeq".srt.bam" ... OK " ||
			( echo "samtools index $nameSeq".srt.bam" ... ERROR " ; exit 1 )
		echo""
			
		if [ $? -eq 0 ] 
		then	
			echo "5. LLAMAR VARIANTES Y LIMPIAR" 
			segmen=$(cat $3"/"$BowRef".bed" | cut -f1)
			for j in $segmen
			do
				echo $j 
				samtools mpileup -aa -A -B -d 0 -Q 10 -r $j -l $3"/"$BowRef".bed" -f $3"/"$BowRef".fasta" $nameSeq".srt.bam" | ivar consensus -q 20 -t 0 -m 5 -n N -p $nameSeq"_"$j".ivar" &&
				echo "samtools mpileup -aa -A -d 0 -Q 10 -r $j -l $3"/"$BowRef".bed" -f $3"/"$BowRef".fasta" $nameSeq".srt.bam" | ivar consensus -q 20 -t 0 -m 5 -n N -p $nameSeq"_"$j".ivar" ... OK" ||
				(echo "samtools mpileup -aa -A -d 0 -Q 10 -r $j -l $3"/"$BowRef".bed" -f $3"/"$BowRef".fasta" $nameSeq".srt.bam" | ivar consensus -q 20 -t 0 -m 5 -n N -p $nameSeq"_"$j".ivar" ... ERROR"; exit 1)
		
				samtools mpileup -A -B -d 0 -Q 10 -r $j -l $3"/"$j".bed" -f $3"/"$BowRef".fasta" $nameSeq".srt.bam" | ivar variants -p $nameSeq"_"$j".ivar" -q 20 -t 0.0 -r $3"/"$BowRef".fasta" -g $3"/"$j".gff3" && 
				echo "samtools mpileup -A -B -d 0 -Q 10 -r $j -l $3"/"$j".bed" -f $3"/"$BowRef".fasta" $nameSeq".srt.bam" | ivar variants -p $nameSeq"_"$j".ivar" -q 20 -t 0.0 -r $3"/"$BowRef".fasta" -g $3"/"$j".gff3" ... OK" ||
				(echo "samtools mpileup -A -B -d 0 -Q 10 -r $j -l $3"/"$j".bed" -f $3"/"$BowRef".fasta" $nameSeq".srt.bam" | ivar variants -p $nameSeq"_"$j".ivar" -q 20 -t 0.0 -r $3"/"$BowRef".fasta" -g $3"/"$j".gff3" ... ERROR"; exit 1)
				perl /scratch/btaboada/Programs/cleanConsensusSeg_v1.pl $nameSeq"_"$j".ivar.fa" $nameSeq"_"$j".ivar.tsv" $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/Consensus_//g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/\.ivar_threshold_0_quality_20//g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/Q_Dp_//g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/A_H3N2_A_H3N2/A_H3N2/g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/A_H1N1_A_H1N1/A_H1N1/g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/B_Yama_B_Yama/B_Yama/g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/B_Vict_B_Vict/B_Vict/g' $nameSeq"_"$j".ivar2.fasta"
				sed -i 's/_S[0-9]\+//' $nameSeq"_"$j".ivar2.fasta"
			done 
			if [ $? -eq 0 ] 
			then
				rm $(ls $nameSeq*.fa)
				zipF=$(ls "$nameSeq"*.fastq"" "$nameSeq"*.sam"" "$nameSeq"*qual.txt"" "$nameSeq"*.ivar2.fasta"" "$nameSeq"*.ivar.tsv"")
				for z in $zipF
				do
					gzip -q $z
				done
				if [ ! -d FASTA ]; then
					mkdir FASTA
				fi
				if [ ! -d SAM ]; then
					mkdir SAM
				fi
				if [ ! -d BAM ]; then
					mkdir BAM
				fi
				if [ ! -d VIRUS ]; then
					mkdir VIRUS
				fi
				zcat $(ls $nameSeq*ivar2.fasta.gz) >$nameSeq"_All_ivar2.fasta"
				mv $(ls $nameSeq*.sam.gz) SAM
				mv $(ls $nameSeq*.bam*) BAM
				mv $(ls $nameSeq*.ivar2.fasta.gz) FASTA
				mv $(ls $nameSeq*.ivar.tsv.gz) FASTA
				mv $(ls $nameSeq*.ivar.vcf.gz) FASTA
				mv $(ls $nameSeq*.ivar.qual.txt.gz) FASTA
				mv $(ls $nameSeq*_M.fastq.gz) VIRUS
				cd FASTA
				rename 'A_H3N2_A_H3N2' 'A_H3N2_A_H3N2' *.gz
				rename 'A_H1N1_A_H1N1' 'A_H1N1_A_H1N1' *.gz
				rename '_B_Vict_B_Vict' '_B_Vict' *.gz
				rename 'B_Yama_B_Yama' 'B_Yama_B'  *.gz
			else
				echo "ERROR EN.. 7) LLAMAR VARIANTES Y LIMPIAR!!"
			fi
		else
			echo "ERROR EN.. 4) ORDENAR BAM!!"
		fi
	else
		echo "ERROR EN.. 3) CONVERTIR SAM!!"
	fi
else
	echo "ERROR. SON 3 PARAMETROS DE ENTRADA!!"
fi
