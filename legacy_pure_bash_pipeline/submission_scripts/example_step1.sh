
#!/bin/bash
#$ -S /bin/bash

date ##to measure the duration
hostname



#!/bin/bash -l
#$ -S /bin/bash
#$ -o ./results/sample1/cluster/out 
#$ -e ./results/sample1/cluster/error
#$ -l h_rt=48:00:00
#$ -l tmem=34.9,h_vmem=34.9G
#$ -R y      
#$ -N step1
#$ -wd  ./results/sample1
#$ -V




/metagenomica/exec/trim_galore_0.3.7/trim_galore -q 20 --length 25 --paired  ./input/sample1_R1.fastq.gz  ./input/sample1_R2.fastq.gz  -o ./results/sample1/trim                                                                                                                                   

mv ./results/sample1/trim/*val_1.fq.gz ./results/sample1/trim/sample1_1_val_1.fq.gz                                                                                           
                                  
mv ./results/sample1/trim/*val_2.fq.gz ./results/sample1/trim/sample1_2_val_2.fq.gz                                                                                                                  

                                                                                                                            




/share/apps/genomics/bowtie2-2.4.1/bowtie2 -p 6 -x .//DBs/human_db_new/bowtie//human_GRCh38.p9_addition -1 ./results/sample1/trim/sample1_1_val_1.fq.gz -2 ./results/sample1/trim/sample1_2_val_2.fq.gz -S ./results/sample1/quickAlign_human/sample1.sam


	                                                                                                                                                                                                   
                      
 /metagenomica/exec/samtools-0.1.19/samtools view -bS  -o ./results/sample1/quickAlign_human/sample1.bam ./results/sample1/quickAlign_human/sample1.sam                                                                                                 
                                                                                                                                                                                                           
                                                     
 /metagenomica/exec/samtools-0.1.19/samtools sort  ./results/sample1/quickAlign_human/sample1.bam ./results/sample1/quickAlign_human/sample1_sorted                                                                                                     
                                                                                                                                                                                                           
                                                              
 /metagenomica/exec/samtools-0.1.19/samtools index  ./results/sample1/quickAlign_human/sample1_sorted.bam                                                                                                                                          
                                                                
 /metagenomica/exec/samtools-0.1.19/samtools view -u -f 12 -F 256 ./results/sample1/quickAlign_human/sample1.bam > ./results/sample1/quickAlign_human/sample1_unmapped.bam

 /share/apps/genomics/bedtools-2.25.0/bin//bamToFastq -i  ./results/sample1/quickAlign_human/sample1_unmapped.bam -fq ./results/sample1/quickAlign_human/sample1_1.fq -fq2 ./results/sample1/quickAlign_human/sample1_2.fq                                                                                                                                                                                                           
                       
                                                                                                                                                                           
  

export PERL5LIB=${PERL5LIB}:/metagenomica/exec/bioperl-live            
                                                                                                                                                                                                              
perl /metagenomica/exec/prinseq-lite-0.20.3/prinseq-lite.pl -fastq ./results/sample1/quickAlign_human/sample1_1.fq -fastq2 ./results/sample1/quickAlign_human/sample1_2.fq  -min_qual_mean 20 -out_good ./results/sample1/QC/Samplesample1  -trim_qual_right 10 -trim_qual_left 10  -lc_method dust -lc_threshold 7 -no_qual_header -qual_noscale -graph_data ./results/sample1/QC/Samplesample1_quality.plot  -out_bad null                                                                                                                                             
                                                                                                                                                                                                             



/metagenomica/exec/seqtk-master/seqtk seq -a ./results/sample1/QC/Samplesample1_1.fastq > ./results/sample1/QC/Samplesample1_1.fasta 
/metagenomica/exec/seqtk-master/seqtk seq -a ./results/sample1/QC/Samplesample1_2.fastq > ./results/sample1/QC/Samplesample1_2.fasta 

cat ./results/sample1/QC/Samplesample1_1.fasta ./results/sample1/QC/Samplesample1_2.fasta > ./results/sample1/QC/sample1_paired.fasta 

rm ./results/sample1/QC/Samplesample1_1.fasta  ./results/sample1/QC/Samplesample1_2.fasta 


  awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print $0}' ./results/sample1/QC/sample1_paired.fasta > ./results/sample1/QC/NM_paired_sample1.txt


  
                                                                                                                                                                                               
/share/apps/genomics/blast-2.9.0/bin/blastn -db .//DBs/human_db_new/human_GRCh38.p9_additions -query  ./results/sample1/QC/sample1_paired.fasta  -outfmt 6  -num_alignments 1   -evalue 1  -culling_limit 1 -num_threads 6  > ./results/sample1/blastn/sample1_paired.ncbiBLASTn


 
                                                                                                                                                                                                                    
                                                                                                                                                                                                      
####### Make fasta files for non overlapping reads (paired end)                                                                                                                                
                                                                                                                                                                                               
	file1=./results/sample1/blastn/sample1_paired.ncbiBLASTn
	file1tmp=./results/sample1/blastn/sample1_paired_ncbi_tmp.txt

	file2=./results/sample1/QC/NM_paired_sample1.txt             
	file2tmp=./results/sample1/QC/NM_paired_sample1_tmp.txt                                                                                                              
                                                                                                                                                                                               
        filtered=./results/sample1/blastn/sample1_paired_filtered.txt                                                                                                                                

                   
        awk  -v x="/"  'BEGIN {FS=OFS="\t"} {split($1, a, x);  print a[1], $1, $3}' ${file1} | uniq > ${file1tmp}
        awk  -v x="/"  'BEGIN {FS=OFS="\t"} {split($1, a, x);  print a[1], $1, $2}' ${file2} | uniq > ${file2tmp}                                                                                                                                                                                                                                                                                     
        echo $file1tmp                                                                                                                                                                    
        echo $file2tmp                                                                                                                                                                                                                                                                                                                                                                                 
        awk  'BEGIN {FS=OFS="\t"} NR==FNR{a[$1]=$1;next} a[$1]!=$1{print $2,$3}' ${file1tmp} ${file2tmp} > ${filtered}                                                     
        wc -l ${filtered} > ./results/sample1/blastn/NM_paired.number
        awk -F"\t" '{print ">"$1"\n"$2}' ${filtered} > ./results/sample1/blastn/NM_paired_sample1_filtered.fasta

                                                                                                                                                                                               
                                                                                                                                                                                               

                                                                                                                                                                                               
/share/apps/genomics/blast-2.9.0/bin/blastn -db .//DBs/UniVec/UniVec -query  ./results/sample1/blastn/NM_paired_sample1_filtered.fasta  -outfmt 6  -num_alignments 1   -evalue 1  -culling_limit 1 -num_threads 6  > ./results/sample1/blastn_univec/sample1_paired.ncbiBLASTn


 
                                                                                                                                                                                                                    
                                                                                                                                                                                                      
####### Make fasta files for non overlapping reads (paired end)                                                                                                                                
                                                                                                                                                                                               
	file1=./results/sample1/blastn_univec/sample1_paired.ncbiBLASTn
	file1tmp=./results/sample1/blastn_univec/sample1_paired_ncbi_tmp.txt

	file2=./results/sample1/blastn/sample1_paired_filtered.txt             
	file2tmp=./results/sample1/blastn/NM_paired_sample1_tmp.txt                                                                                                              
                                                                                                                                                                                               
        filtered=./results/sample1/blastn_univec/sample1_paired_filtered.txt                                                                                                                                

                   
        awk  -v x="/"  'BEGIN {FS=OFS="\t"} {split($1, a, x);  print a[1], $1, $3}' ${file1} | uniq > ${file1tmp}
        awk  -v x="/"  'BEGIN {FS=OFS="\t"} {split($1, a, x);  print a[1], $1, $2}' ${file2} | uniq > ${file2tmp}                                                                                                                                                                                                                                                                                     
        echo $file1tmp                                                                                                                                                                    
        echo $file2tmp                                                                                                                                                                                                                                                                                                                                                                                 
        awk  'BEGIN {FS=OFS="\t"} NR==FNR{a[$1]=$1;next} a[$1]!=$1{print $2,$3}' ${file1tmp} ${file2tmp} > ${filtered}                                                     
        wc -l ${filtered} > ./results/sample1/blastn_univec/NM_paired.number
        awk -F"\t" '{print ">"$1"\n"$2}' ${filtered} > ./results/sample1/blastn_univec/NM_paired_sample1_filtered.fasta

                                                                                                                                                                                               
                                                                                                                                                                                               




# ##############megaBLAST for nucleotides

#/share/apps/genomics/blast-2.9.0/bin/blastn -db .//DBs/refseq_nucleotide_28March2020/custom_nucleotide_refseq_nounplaced_July2019_March2020.fasta -query ./results/sample1/blastn_univec/NM_paired_sample1_filtered.fasta -outfmt 6 -max_target_seqs 10 -max_hsps 1 -num_threads 6  >  ./results/sample1/nucleotide/sample1_megaBLAST.tab





filtered=./results/sample1/blastn/sample1_paired_filtered.txt                                                                                                                                


 awk  'BEGIN {FS=OFS="\t"} {print $1}'  ${filtered}  > ./results/sample1/blastn/ids.txt                                                    
       

 /metagenomica/exec/seqtk-master/seqtk subseq ./results/sample1/QC/Samplesample1_1.fastq ./results/sample1/blastn/ids.txt > ./results/sample1/blastn/sample1_filtered_1.fastq
 /metagenomica/exec/seqtk-master/seqtk subseq ./results/sample1/QC/Samplesample1_2.fastq ./results/sample1/blastn/ids.txt > ./results/sample1/blastn/sample1_filtered_2.fastq






/share/apps/genomics/bowtie2-2.4.1/bowtie2 -p 6 -x .//DBs/silva_rRNA/new_db/bowtie//SILVA_128_SSU_LSU_UniVec -1 ./results/sample1/blastn/sample1_filtered_1.fastq -2 ./results/sample1/blastn/sample1_filtered_2.fastq -S ./results/sample1/quickAlign_rRNA/sample1.sam


	                                                                                                                                                                                                   
                      
/metagenomica/exec/samtools-0.1.19/samtools view -bS  -o ./results/sample1/quickAlign_rRNA/sample1.bam ./results/sample1/quickAlign_rRNA/sample1.sam                                                                                                 
                             
                                                                                                                                                                                                                                    
/metagenomica/exec/samtools-0.1.19/samtools sort  ./results/sample1/quickAlign_rRNA/sample1.bam ./results/sample1/quickAlign_rRNA/sample1_sorted                                                                                                                                           
                                                                                                                                                                                                                                     
/metagenomica/exec/samtools-0.1.19/samtools index  ./results/sample1/quickAlign_rRNA/sample1_sorted.bam                                                                                                                                                                                
/metagenomica/exec/samtools-0.1.19/samtools view -u -f 12 -F 256 ./results/sample1/quickAlign_rRNA/sample1.bam > ./results/sample1/quickAlign_rRNA/sample1_unmapped.bam

/share/apps/genomics/bedtools-2.25.0/bin//bamToFastq -i  ./results/sample1/quickAlign_rRNA/sample1_unmapped.bam -fq ./results/sample1/quickAlign_rRNA/sample1_1.fq -fq2 ./results/sample1/quickAlign_rRNA/sample1_2.fq                                                                                                                                                                                                                                                 

/metagenomica/exec/seqtk-master/seqtk seq -a ./results/sample1/quickAlign_rRNA/sample1_1.fq > ./results/sample1/quickAlign_rRNA/sample1_1.fasta 
/metagenomica/exec/seqtk-master/seqtk seq -a ./results/sample1/quickAlign_rRNA/sample1_2.fq > ./results/sample1/quickAlign_rRNA/sample1_2.fasta 
                     
cat  ./results/sample1/quickAlign_rRNA/sample1_1.fasta  ./results/sample1/quickAlign_rRNA/sample1_2.fasta > ./results/sample1/quickAlign_rRNA/sample1_paired.fasta 


 awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print $0}' ./results/sample1/quickAlign_rRNA/sample1_paired.fasta > ./results/sample1/quickAlign_rRNA/sample1_paired_filtered.txt

                                                                                                                                                                           

                                                                                                                                                                                               
/share/apps/genomics/blast-2.9.0/bin/blastn -db .//DBs/silva_rRNA/new_db/SILVA_128_SSU_LSU_UniVec -query ./results/sample1/quickAlign_rRNA/sample1_paired.fasta  -outfmt 6  -num_alignments 1   -evalue 0.1  -culling_limit 1 -num_threads 6  > ./results/sample1/blastn_rRNA/sample1_paired.ncbiBLASTn



                                                                                                                                                                                             

                         
                                                                                                                                                                                               
####### Make fasta files for non overlapping reads (paired end)                                                                                                                                
                                                                                                                                                                                               
	file1=./results/sample1/blastn_rRNA/sample1_paired.ncbiBLASTn
	file1tmp=./results/sample1/blastn_rRNA/sample1_paired_ncbi_tmp.txt

	file2=./results/sample1/quickAlign_rRNA/sample1_paired_filtered.txt             
	file2tmp=./results/sample1/quickAlign_rRNA/sample1_paired_filtered_tmp.txt                                                                                                              
                                                                                                                                                                                               
        filtered=./results/sample1/blastn_rRNA/sample1_paired_filtered_rRNA.txt                                                                                                                                

                   
        awk  -v x="/"  'BEGIN {FS=OFS="\t"} {split($1, a, x);  print a[1], $1, $3}' ${file1} | uniq > ${file1tmp}
	awk  -v x="/"  'BEGIN {FS=OFS="\t"} {split($1, a, x);  print a[1], $1, $2}' ${file2} | uniq > ${file2tmp}                                                                                                                                                                                                                                                                                     
        echo $file1tmp                                                                                                                                                                    
        echo $file2tmp                                                                                                                                                                                                                                                                                                                                                                                 
        awk  'BEGIN {FS=OFS="\t"} NR==FNR{a[$1]=$1;next} a[$1]!=$1{print $2,$3}' ${file1tmp} ${file2tmp} > ${filtered}                                                     

        wc -l ${filtered} > ./results/sample1/blastn_rRNA/NM_paired.number

        awk -F"\t" '{print ">"$1"\n"$2}' ${filtered} > ./results/sample1/blastn_rRNA/NM_paired_sample1_filtered_rRNA.fasta







# ##############megaBLAST for nucleotides

/share/apps/genomics/blast-2.9.0/bin/blastn -db .//DBs/refseq_nucleotide_28March2020/custom_nucleotide_refseq_nounplaced_July2019_March2020.fasta -query ./results/sample1/blastn_rRNA/NM_paired_sample1_filtered_rRNA.fasta -outfmt 6 -max_target_seqs 10 -max_hsps 1 -num_threads 6  >  ./results/sample1/nucleotide/sample1_megaBLAST_norRNA.tab





/metagenomica/exec/diamond blastx -d .//DBs/refseq_protein_28March2020/protein_march2020.refseq.faa.dmnd -q ./results/sample1/blastn_rRNA/NM_paired_sample1_filtered_rRNA.fasta   -f 6  -k 10 -p 6 -sensitive -e 1  > ./results/sample1/protein/sample1_diamond.tab






# ### remove files that take up too much storage space and won't be used downstream

# ### from bowtie folders
 rm   ./results/sample1/quickAlign_human/sample1.sam

 rm  ./results/sample1/quickAlign_rRNA/sample1.sam

# ##from QC folder
# rm ./results/sample1/QC/NM_paired_sample1.txt
# rm ./results/sample1/QC/NM_paired_sample1_tmp.txt

# ## from blastn folders
# rm ./results/sample1/blastn/sample1_paired_ncbi_tmp.txt                  


 

                                                                                                                                                                                         
grep 'Processed reads:' ./results/sample1/trim/sample1_1.fastq.gz_trimming_report.txt | awk -F":" 'BEGIN {OFS="\t"}{print "Raw reads", $2*2}' > ./results/sample1/summary/sample1_raw.count   

grep 'Trimmed reads:' ./results/sample1/trim/sample1_1.fastq.gz_trimming_report.txt | awk -F":" 'BEGIN {OFS="\t"}{print "Trimmed reads", $2*2}' > ./results/sample1/summary/sample1_trimmed.count   
                                                                                                                                                                                                   
wc -l  ./results/sample1/quickAlign_human/sample1_1.fq |  awk  'BEGIN {OFS="\t"}{ print "Non human (bowtie2)", ($1/4)*2}' > ./results/sample1/summary/sample1_nonhuman.count           

                                                                                                                                                                                              
wc -l  ./results/sample1/quickAlign_rRNA/sample1_1.fq |  awk 'BEGIN {OFS="\t"}{ print "Non rRNA (bowtie2)", ($1/4)*2}' > ./results/sample1/summary/sample1_nonrRNA.count              


wc -l  /sample1_paired.fasta | awk 'BEGIN {OFS="\t"}{ print "after QC", /2}' > ./results/sample1/summary/sample1_QC.count                                                                                                                                                                                                                                                                                                                                                                                                             
awk '{s+=$1} END {print s}' ./results/sample1/blastn/NM_paired.number | awk -F":" 'BEGIN {OFS="\t"}{ print "non human (blastn)", $1}'  > ./results/sample1/summary/sample1_nonhuman_blastn.count            
                                                                                                                                                                                                      
awk '{s+=$1} END {print s}' ./results/sample1/blastn_rRNA/NM_paired.number | awk -F":" 'BEGIN {OFS="\t"}{ print "non rRNA (blastn)", $1}'  > ./results/sample1/summary/sample1_nonrRNA_blastn.count                 
                                                                                                                                                                                                      
awk -F"\t" '{print $1}' ./results/sample1/protein/sample1_diamond.tab | sort | uniq | wc -l | awk -F" " '{print "Reads mapping to prot DB", $1}' > ./results/sample1/summary/sample1_maptoprot.count     
                                                                                                                                                                                                     
                                                          
awk -F"\t" '{print $1}' ./results/sample1/nucleotide/sample1_megaBLAST.tab | sort | uniq | wc -l | awk -F" " '{print "Reads mapping to nucl DB", $1}' > ./results/sample1/summary/sample1_maptonucl.count                                                                                                                                                                                                      
                                                                                                                                                                                                     
                                                                 
cat ./results/sample1/summary/sample1_raw.count ./results/sample1/summary/sample1_trimmed.count ./results/sample1/summary/sample1_nonhuman.count ./results/sample1/summary/sample1_nonrRNA.count  ./results/sample1/summary/sample1_QC.count ./results/sample1/summary/sample1_nonhuman_blastn.count ./results/sample1/summary/sample1_nonrRNA_blastn.count ./results/sample1/summary/sample1_maptoprot.count ./results/sample1/summary/sample1_maptonucl.count  > ./results/sample1/summary/sample1_summary.count.new                                                                                                                                                         



 
