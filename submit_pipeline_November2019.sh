mainDir=/blah/blah/blah     #### This needs to be changed to the directory where the main pipeline script is as well as the DBs


pipeline=${mainDir}/scripts/pipeline_RNAseq_November2019_pchuckle.sh   ####Main pipeline script

###metagenomic databases  ###
 
referenceHuman=${mainDir}/DBs/human_db_new/bowtie/                            ####### bowtie2 index for human reference

referencerRNA=${mainDir}/DBs/silva_rRNA/new_db/bowtie/                            ####### bowtie2 index for rRNA reference

dbHuman=${mainDir}/DBs/human_db_new/human_GRCh38.p9_additions                     ####### blast database Human                                                      
                 
dbrRNA=${mainDir}/DBs/silva_rRNA/new_db/SILVA_128_SSU_LSU_UniVec            ###########blast DB rRNA
dbunivec=${mainDir}/DBs/UniVec/UniVec                                      ##### blast DB vectors

dbcustom=${mainDir}/DBs/refseq_protein_28March2020/protein_march2020.refseq.faa.dmnd   ###DIAMOND index for protein DB

dbcustomnucl=${mainDir}/DBs/refseq_nucleotide_28March2020/custom_nucleotide_refseq_nounplaced_July2019_March2020.fasta  ###megablast index for nucleotide DB 

delim="/"



results=/cluster/scratch8b/Morfopoulou_encephalitis/   ##### results directory

if [ ! -e $results ]; then mkdir $results; fi    #### if results dir doesn't exist, create dir


 for sample in `cat sampleNames.tab`; do  
     echo $sample
     mkdir $results/$sample

 done



for sample in `cat sampleNames.tab`; do


outDir=${results}/${sample}                                                                                     ####output directory for specific sample

inputFiles="2  /blah/blah/${sample}*1*.fastq.gz  /blah/blah/${sample}*2*.fastq.gz"     ####input fastqs, usually in backed up directories

scriptStep1=${outDir}/cluster/submission/step1.sh    #### script that will be submitted to the job queue


##STEP1                                                                                                                                                                                              

bash ${pipeline} --script ${scriptStep1} --inputFiles ${inputFiles}  --outDir ${outDir} --sample ${sample} --delimiter ${delim}  --reference ${referenceHuman} --referencerRNA ${referencerRNA} --db ${dbHuman} --dbrRNA ${dbrRNA} --dbunivec ${dbunivec} --dbprot ${dbcustom} --dbnucl ${dbcustomnucl} --step1 TRUE




done





