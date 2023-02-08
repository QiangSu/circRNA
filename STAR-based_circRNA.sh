##1.pathway
path_raw_read=/home/zhouwanting/shujin/90_sample-Geneplus/raw_read
path_raw_cutadapt=/home/suqiang/90-samples/raw_cutadapt
path_raw_trimmomatic=/home/suqiang/90-samples/raw_trimmomatic
path_qc_result=/home/suqiang/90-samples/qc_result
path_rseqc_result=/home/suqiang/90-samples/rseqc_result
bam_dir=/home/suqiang/90-samples/bam_reads
##2.cutadapt
for i in $(cat $path_raw_read/samplelist.txt)
do
        cutadapt -j 20 -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGAC -A AAGTCGGATCGTAGCCATGTCGTTC -o ${path_raw_cutadapt}/${i}_cutadapt_R1.fastq.gz -p ${path_raw_cutadapt}/${i}_cutadapt_R2.fastq.gz $path_raw_read/${i}_raw_1.fq.gz $path_raw_read/${i}_raw_2.fq.gz
done 
echo cut adapter done

##3.trim
for i in $(cat $path_raw_read/samplelist.txt);
do
trimmomatic PE -threads 15 $path_raw_cutadapt/${i}_cutadapt_R1.fastq.gz $path_raw_cutadapt/${i}_cutadapt_R2.fastq.gz -baseout $path_raw_trimmomatic/${i}_cutadapt_trim.fastq.gz LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20; done

##4.STAR
for i in $(cat $path_raw_read/samplelist.txt)
do
        STAR --runThreadN 10 --genomeDir /home/suqiang/ref/STAR_index/ --readFilesIn $path_raw_trimmomatic/${i}_cutadapt_trim_1P.fastq.gz $path_raw_trimmomatic/${i}_cutadapt_trim_2P.fastq.gz --sjdbOverhang 149 --outFileNamePrefix ${bam_dir}/${i}- --outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --chimOutType Junctions SeparateSAMold --chimSegmentMin 10 --readFilesCommand zcat

done
echo STAR done

##5. Isoforms determination
Cufflinks 
for i in $(cat /home/suqiang/90-samples/bam_reads/samplelist-left-20230111.txt)
do

cufflinks -p 10 --library-type fr-firststrand -G /home/suqiang/ref/gencode.v42.annotation.gtf -o ./cufflinks/${i} ./${i}-Aligned.sortedByCoord.out.bam 
done


##6. CircRNA identification
perl /home/suqiang/circRNA_finder-master/postProcessStarAlignment.pl --starDir /home/suqiang/90-samples/bam_reads/ --minLen 100 --outDir /home/suqiang/90-samples/bam_reads/circRNA_finder/

##7.Summing up the item with same ID(removing duplicates)
 
for i in $(cat samplelist.txt)
do
awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' ${i}.csv >> ./unique/${i}-uniqu.csv
done


