#### This needs to be changed to the directory where the main pipeline script is as well as the DBs; it also look like you need to change the inputFiles variable (fastq files)

# Main pipeline script
pipeline=/app/pipeline_RNAseq_November2019_pchuckle.sh

# metagenomic databases
referenceHuman=${mainDir}/DBs/human_db_new/bowtie/                                                                      # bowtie2 index for human reference
referencerRNA=${mainDir}/DBs/silva_rRNA/new_db/bowtie/                                                                  # bowtie2 index for rRNA reference
dbHuman=${mainDir}/DBs/human_db_new/human_GRCh38.p9_additions                                                           # blast database Human                                                      
dbrRNA=${mainDir}/DBs/silva_rRNA/new_db/SILVA_128_SSU_LSU_UniVec                                                        # blast DB rRNA
dbunivec=${mainDir}/DBs/UniVec/UniVec                                                                                   # blast DB vectors
dbcustom=${mainDir}/DBs/refseq_protein_28March2020/protein_march2020.refseq.faa.dmnd                                    # DIAMOND index for protein DB
dbcustomnucl=${mainDir}/DBs/refseq_nucleotide_28March2020/custom_nucleotide_refseq_nounplaced_July2019_March2020.fasta  # megablast index for nucleotide DB 

# create results dir
results=/output
if [ ! -e $results ]; then mkdir $results; fi

samples=$(seq 1 1)

# make per-sample results dir
for sample in $samples; do
     echo $sample
     mkdir $results/$sample
done

for sample in $samples; do
    outDir=${results}/${sample}
    inputFiles="2 /input/${sample}*1*.fastq.gz /input/${sample}*2*.fastq.gz" # input fastqs, usually in backed aup directories
    scriptStep1=/output/generated_scripts.sh

    #generate and submit script to job queue 
    bash ${pipeline} --script ${scriptStep1} \
                     --inputFiles ${inputFiles} \
                     --outDir ${outDir} \
                     --sample ${sample} \
                     --delimiter ${delim} \
                     --reference ${referenceHuman} \
                     --referencerRNA ${referencerRNA} \
                     --db ${dbHuman} \
                     --dbrRNA ${dbrRNA} \
                     --dbunivec ${dbunivec} \
                     --dbprot ${dbcustom} \
                     --dbnucl ${dbcustomnucl} \
                     --step1 TRUE
done





