#################################################################
#################     USER SCRIPT START     #####################
#################################################################

start_time=`date +%s`
echo "start time:" $start_time

#################################################################
#################   DECLARATIONS / PATHS   ######################
#################################################################
proj="rice" #manually set job name here before running, one project can have multiole jobs
job="heat_drought" #manually set job name here before running


fastq_f="/mnt/d/rna_seq/ena_files" #fastq files (source) directory (no trailing /)
ref_g="/mnt/d/rna_seq/ref_gen/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa" #full path of your (source) reference genome file ending in .fa OR .fasta
gtf_f="/mnt/d/rna_seq/gtf/Oryza_sativa.IRGSP-1.0.50.chr.gtf" #full path of (source) gtf file ending in .gtf

#################################################################
################ Start of Pipline workflow_PE ###################
#################################################################
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----STARTING THE PIPELINE------"

#:<<'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "----------Creating Directory hierachy------"
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastqc/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastq_files/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastp/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastp_reports/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/featureCounts/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/hisat2/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/samtools/

mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/reference_genome/
mkdir -p /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/gtf/

#END

fastqc=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastqc
fastq_files=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastq_files
fastp=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastp
fastp_reports=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/fastp_reports
featureCounts=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/featureCounts
hisat2=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/hisat2
samtools=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/results/samtools
reference_genome=/mnt/f/projects/$proj/pipeline_result/$job/data/workflow_SE/reference_genome


#ADD A HASH BEFORE THE COLON ABOVE IF INPUT FILES ARE SRA (ALSO FIND THE CORRESPONDING TAG BELOW)
#REMOVE THE HASH (IF ANY) BEFORE THE COLON ABOVE IF INPUT FILES ARE FASTQ/FASTQ.GZ (ALSO FIND THE CORRESPONDING TAG BELOW)

#echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Running fastq-dump-----------"

:<<'END'

cd $fastq_files



##SRA to Fastq files

for file in $(ls $sra)
do

fasterq-dump $sra/$file -v -p -b 100 -c 1024 -m 10240 -e 50       #resource limit

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Conversion to FastQ complete. Compressing now.-----------"



pigz -v *.fastq #multhreaded gzip #resource limit

END

#echo [`date +"%Y-%m-%d %H:%M:%S"`] "-----All SRAs split and converted to fastq.gz----"

#ADD A HASH BEFORE THE END ABOVE IF INPUT FILES ARE SRA (ALSO FIND THE CORRESPONDING TAG ABOVE)
#REMOVE THE HASH (IF ANY) BEFORE THE END ABOVE IF INPUT FILES ARE FASTQ/FASTQ.GZ (ALSO FIND THE CORRESPONDING TAG ABOVE)

##Trimming by fastp
echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Trimming by fastp------"

cd $fastp

find $fastq_f -name "*.fastq.gz" | sort | while read A
do

	a=`basename ${A} | sed 's/.sra_1/_1/' | awk -F "." '{print $1}'`
	echo "$a"
	echo "Processing $a"

	fastp --thread=2 --length_required=10 --qualified_quality_phred=32 --in1=${A} --out1=$a\_trimmed.fastq.gz --json=$a.json --html=$a.html
  #--thread= number of worker threads (max 16)
  #USE your required adapter after -a, default = automatic detection

	a1=$a
	a1+=".html"
	a2=$a
	a2+=".json"
	mv -v $a.json $fastp_reports/$a2 | awk -F "_" '{print $1".json"}'`
	mv -v $a.html $fastp_reports/$a1 | awk -F "_" '{print $1".html"}'`

	echo "$a. Report generated and moved to results fastp_reports!!!"

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------All trimming completed!!!----"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing QC after trimming------"

cd $fastp

f_ls="$(ls)"

fastqc -q -t 2 $f_ls #resource limit
#-t = number of threads (250 MB memory required per thread)

ls | awk -F "." '{print $1}' | uniq | while read report
do

	mv -v $report*.html $fastqc/`echo $report | awk -F "." '{print $1}'`.html
	mv -v $report*.zip $fastqc/`echo $report | awk -F "." '{print $1}'`.zip

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc!!----------"


#END

#: << 'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Preparing for Alignment----------"

##Index building and Read alignment using hisat2## 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Locating and placing reference genome------------"

cp -v $ref_g $reference_genome/ref_genome.fa
cd $reference_genome

echo [`date +"%Y-%m-%d %H:%M:%S"`] "-----Building indices-----"

##building index
hisat2-build -p 2 $reference_genome/ref_genome.fa index #resource limit
#-p = number of threads

#END

##Doing Alignment
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----------Aligning with indices--------"

cd $fastp

find $fastp -name "*_trimmed.fastq.gz" | sort | while read A 

do

	a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`
	hisat2 --threads 2 --dta -x $reference_genome/index -U ${A} -S $hisat2/$a.sam
 	#--threads = number of simultaneous alignments
done

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done alignment and placed SAM files in results hisat2------"

#END

##-----------------------------------------##

#: << 'END'

##Converting sam files to bam files using SAMtools

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Running SAM Tools to Convert SAM to BAM---------"

cd $hisat2

for file in $(ls)
do

	a=`echo ${file} | awk -F "." '{print $1}'`

	echo "Processing ${file}"
	samtools sort -@ 2 -o $a.bam ${file} 	#resource limit
	#-@ number of threads (in addition to main thread)
	echo "${file} converted"

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "-----Converted all SAM files to BAM-----"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Moving BAM files to results samtools-----------"

mv -v *.bam $samtools/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Moved BAM files-----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Done with SAM tools !!!---------"

END

## ---------------------------------------------------------##


###FeatureCount tool

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Generating feature counts...-------------"

cd $samtools/

featureCounts -T 2 -t gene -g gene_id -a $gtf_f -o counts.txt -M *.bam	#resource limit
	#-T number of threads

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Done Generating count data-----------"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Moving Results...-----------"

mv -v *.txt $featureCounts/
mv -v *.summary $featureCounts/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Results moved to results featureCounts---------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Finished with featureCounts!!!-----"

#END

##-----------------------------------------##

### DeSeq2 in R

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Calling R for DeSeq2-------------"

:<<'END'

cp $deseq2 /mnt/f/projects/$proj/pipeline_result/$job/scripts/deseq2.R
cd /mnt/f/projects/$proj/pipeline_result/$job/scripts/

r deseq2.R

mv -v res.csv bak_res.csv
echo -n "", > res.csv; cat bak_res.csv >> res.csv #fixes the left shift of column names
rm bak_res.csv

mv -v res_PAdj_cutoff.csv bak_res_PAdj_cutoff.csv
echo -n "", > res_PAdj_cutoff.csv; cat bak_res_PAdj_cutoff.csv >> res_PAdj_cutoff.csv #fixes the left shift of column names
rm bak_res_PAdj_cutoff.csv

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------DeSeq2 Complete-------------"

END

##-----------------------------------------------##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Cleaning Up-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting SAM files-------------"
rm -r /mnt/f/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------SAM files deleted!!!-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Clean Up Completed !!!-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------END OF PIPELINE---------"

end_time=`date +%s`
echo "############################################################"
echo
echo "time taken:"
echo
echo "start time:" $start_time
echo "end time:" $end_time
echo
run_time=$((end_time-start_time))
echo "runtime (in seconds):" $run_time
echo -n "runtime (in minutes): "; awk "BEGIN {print $run_time/60}"
echo -n "runtime (in hours): "; awk "BEGIN {print $run_time/3600}"
echo
echo "############################################################"

exit

############################################################
#################### END OF SCRIPT #########################
############################################################
