source /share/apps/source_files/metaMix.source
results=/cluster/scratch8b/Morfopoulou_encephalitis/
mainDir=/blah/blah/blah

   for sample in `cat sampleNames.tab`; do

       echo $sample
       	mkdir ${results}/${sample}/metamix_nucleotide_norRNA
       	mkdir ${results}/${sample}/metamix_nucleotide_norRNA/cluster
       	mkdir ${results}/${sample}/metamix_nucleotide_norRNA/cluster/error
       	mkdir ${results}/${sample}/metamix_nucleotide_norRNA/cluster/out

       	awk 'BEGIN{FS=OFS="\t"} {print $1, length($2)}' ${results}/${sample}/blastn_rRNA/${sample}_paired_filtered_rRNA.txt > ${results}/${sample}/metamix_nucleotide_norRNA/read_lengths.tab


      echo "
library(metaMix)

step1<-generative.prob.nucl(blast.output.file=\"${results}/${sample}/nucleotide/${sample}_megaBLAST_norRNA.tab\",  outDir=\"./\",  blast.default=TRUE, read.length.file=\"${results}/${sample}/metamix_nucleotide_norRNA/read_lengths.tab\", contig.weight.file=1,  accession.taxon.file=\"${mainDir}/DBs/refseq_nucleotide_28March2020/custom_nucleotide_refseq_nounplaced_July2019_March2020_nucl_gb.accession2taxid\", genomeLength=\"${mainDir}/DBs/refseq_nucleotide_28March2020//custom_nucleotide_refseq_nounplaced_July2019_March2020_nuclLength.tab\")

step2<-reduce.space(step1=\"step1.RData\", taxon.name.map=\"${mainDir}/DBs/refseq_nucleotide_28March2020/taxonomy/names.dmp\")    

step3<-parallel.temper.nucl(step2=\"step2.RData\",  readSupport=10)

step4<-bayes.model.aver(step2=\"step2.RData\", step3=\"step3.RData\",  taxon.name.map=\"/${mainDir}/DBs/refseq_nucleotide_28March2020/taxonomy/names.dmp\")    

" > ${results}/${sample}/metamix_nucleotide_norRNA/Rsubmit.R



     echo "


#!/bin/bash
#$ -S /bin/bash
#$ -o ${results}/${sample}/metamix_nucleotide_norRNA/cluster/out
#$ -e ${resuts}/${sample}/metamix_nucleotide_norRNA/cluster/error
#$ -wd ${results}/${sample}/metamix_nucleotide_norRNA 
#$ -pe smp 6
#$ -l tmem=7.1G,h_vmem=7.1G
#$ -l h_rt=72:00:00
#$ -V
#$ -R y

source /share/apps/source_files/metaMix.source

mpirun --mca btl_tcp_if_include 128.41.97.0/21 -np 1 /share/apps/genomics/metaMix/R-3.5.2/bin/R --slave CMD BATCH --no-save --no-restore  Rsubmit.R step.out

  " > ${results}/${sample}/metamix_nucleotide_norRNA/submit_metamix_${sample}.sh

      qsub ${results}/${sample}/metamix_nucleotide_norRNA/submit_metamix_${sample}.sh


 done


