source /share/apps/source_files/metaMix.source
results=/cluster/scratch8b/Morfopoulou_encephalitis/
mainDir=/blah/blah/blah

   for sample in `cat sampleNames.tab`; do

       echo $sample
       mkdir ${results}/${sample}/metamix_protein_
       mkdir ${results}/${sample}/metamix_protein/cluster
       mkdir ${results}/${sample}/metamix_protein/cluster/error
       mkdir ${results}/${sample}/metamix_protein/cluster/out

       awk 'BEGIN{FS=OFS="\t"} {print $1, length($2)}' ${results}/${sample}/blastn_rRNA/${sample}_paired_filtered_rRNA.txt > ${results}/${sample}/metamix_protein/read_lengths.tab


      echo "
library(metaMixFast)

step1<-generative.prob(blast.output.file=\"${results}/${sample}/protein/${sample}_diamond.tab\",  outDir=\"./\",  blast.default=TRUE, read.length.file=\"${results}/${sample}/metamix_protein/read_lengths.tab\", contig.weight.file=1, gi.or.acc=\"acc\", accession.taxon.file=\"${mainDir}/DBs/refseq_protein_28March2020/protein_march2020.refseq_prot_gb.accession2taxid\")

step2<-reduce.space(step1=\"step1.RData\", taxon.name.map=\"${mainDir}/DBs/${refseq_protein_28March2020/taxonomy/names.dmp\")             

step3<-parallel.temper(step2=\"step2.RData\", readSupport=10)

step4<-bayes.model.aver(step2=\"step2.RData\", step3=\"step3.RData\", taxon.name.map==\"${mainDir}/DBs/refseq_protein_28March2020/taxonomy/names.dmp\")



" > ${results}/${sample}/metamix_protein/Rsubmit.R



     echo "


#!/bin/bash
#$ -S /bin/bash
#$ -o /cluster/scratch8b/Morfopoulou_encephalitis/${sample}/metamix_protein/cluster/out
#$ -e /cluster/scratch8b/Morfopoulou_encephalitis/${sample}/metamix_protein/cluster/error
#$ -wd /cluster/scratch8b/Morfopoulou_encephalitis/${sample}/metamix_protein 
#$ -pe smp 6
#$ -l tmem=6.1G,h_vmem=6.1G
#$ -l h_rt=48:00:00
#$ -V
#$ -R y



source /share/apps/source_files/metaMix.source

mpirun --mca btl_tcp_if_include 128.41.97.0/21 -np 1 /share/apps/genomics/metaMix/R-3.5.2/bin/R --slave CMD BATCH --no-save --no-restore  Rsubmit.R step.out

  " > ${results}/${sample}/metamix_protein/submit_metamix_${sample}.sh

      qsub ${results}/${sample}/metamix_protein/submit_metamix_${sample}.sh



 done


















