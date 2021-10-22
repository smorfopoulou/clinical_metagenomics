########## everything else below should be automated (if all executables are in the same directory)
diamond=diamond
trimG=trim_galore
seqtk=seqtk
samtools=samtools
prinseq_script=prinseq-lite.pl
bedtools=bedtools
bowtie=bowtie2
blastn=blastn
  
############ default values
script=${output}/cluster/submission/default.sh
step1=FALSE

until [ -z "$1" ]; do
	# use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--inputFiles)
	    shift
	    i=0
	    for fileloc in $@; do 
		inputFiles[ $i ]=$fileloc
		((i=i+1))
	    done;;
	--outDir )
	    shift
	    output=$1;;
	--delimiter )
	    shift
	    delimiter=$1;;
	--script)
	    shift
	    script=$1;;
	 --reference)
            shift
            reference=$1;;
	 --referencerRNA)
            shift
            referencerRNA=$1;;
	 --db )
	   shift
	   db=$1;;
  	 --dbunivec )
	   shift
	   dbunivec=$1;;
	 --dbrRNA )
	   shift
	   dbrRNA=$1;;
	 --dbnucl)
	    shift
	    dbnucl=$1;;
	 --dbprot)
	    shift
	    dbprot=$1;;
         --sample )
           shift
           sample=$1;;
	 --step1 )
	   shift
	   step1=$1;;
	-* )
	    echo "Unrecognized option: $1"
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done

#creating all the output folders
echo -e "Output folder: $output\n"
clusterDir=${output}/cluster
cluster_out=${clusterDir}/out
cluster_error=${clusterDir}/error
cluster_submission=${clusterDir}/submission
output_qc=${output}/QC
output_trim=${output}/trim
output_quickAlign_human=${output}/quickAlign_human
output_quickAlign_rRNA=${output}/quickAlign_rRNA
output_blastn=${output}/blastn
output_univec_blastn=${output}/blastn_univec
output_rRNA=${output}/blastn_rRNA
output_protein=${output}/protein
output_nucleotide=${output}/nucleotide
output_summary=${output}/summary

myFolders="$output $clusterDir $cluster_out $cluster_error $cluster_submission $output_trim $output_quickAlign_human $output_quickAlign_rRNA $output_qc $output_blastn $output_univec_blastn $output_rRNA $output_protein $output_nucleotide $output_summary "
for folder in $myFolders; do
    if [ ! -e $folder ]; then 
	echo "Creating $folder"
	mkdir $folder
    fi
done

# Now writing the script
echo "Output script:  $script"
echo "
#!/bin/bash
#$ -S /bin/bash
date ##to measure the duration
hostname
" > $script

# STEP1 (Trim, bowtie2, prinseq, megablastn to remove host & rRNA, Diamond for protein, megaBLAST for nucleotide)
if [[ "$step1" == "TRUE" ]]; then
    nfiles=${inputFiles[0]}
    if [[ "$nfiles" != "2" ]]; then
        echo "You specified $nfiles input files".
        echo "Error: currently the input data MUST be paired end."
        exit;
    fi
    seq1=${inputFiles[ 1 ]}
    seq2=${inputFiles[ 2 ]}
	
    # check that raw data files & reference exist
    for file in $seq1 $seq2 $reference; do
        ls -lh $file
        if [ ! -e "$file" ]; then
            echo "Error, file $file does not exist"
            exit
        fi
    done

    echo "#1) Trim, bowtie2, QC, remove human with blast, univec for nucleotide/rRNA for protein, DIAMOND and megaBLAST to protein/nucleotide DBs" >> $script

#change the time depending on the size of dataset

    echo "
    #!/bin/bash -l
    #$ -S /bin/bash
    #$ -o $cluster_out
    #$ -e $cluster_error
    #$ -l h_rt=48:00:00
    #$ -nl tmem=34.9,h_vmem=34.9G
    #$ -R y
    #$ -N step1
    #$ -wd  ${output}
    #$ -V
    " >> $script

    echo -e  "#1) Trim adapters from reads and basic quality control"  >> $script
    echo "
    ${trimG} -q 20 --length 25 --paired  $seq1  $seq2  -o ${output_trim}
    mv ${output_trim}/*val_1.fq.gz ${output_trim}/${sample}_1_val_1.fq.gz
    mv ${output_trim}/*val_2.fq.gz ${output_trim}/${sample}_2_val_2.fq.gz
    " >> $script
     
    echo "#2) Bowtie2 - human"  >> $script
    echo "
    ${bowtie} -p 6 -x ${reference}/human_GRCh38.p9_addition -1 ${output_trim}/${sample}_1_val_1.fq.gz -2 ${output_trim}/${sample}_2_val_2.fq.gz -S ${output_quickAlign_human}/${sample}.sam
    " >> $script

    echo -e "#3) Select pairs of reads that are both non mapping" >> $script
    echo "
    ${samtools} view -bS  -o ${output_quickAlign_human}/${sample}.bam ${output_quickAlign_human}/${sample}.sam
    ${samtools} sort  ${output_quickAlign_human}/${sample}.bam ${output_quickAlign_human}/${sample}_sorted
    ${samtools} index  ${output_quickAlign_human}/${sample}_sorted.bam
    ${samtools} view -u -f 12 -F 256 ${output_quickAlign_human}/${sample}.bam > ${output_quickAlign_human}/${sample}_unmapped.bam
    ${bedtools} bamtofastq -i  ${output_quickAlign_human}/${sample}_unmapped.bam -fq ${output_quickAlign_human}/${sample}_1.fq -fq2 ${output_quickAlign_human}/${sample}_2.fq
    " >> $script

    echo "#6) PrinSeq reads" >> $script
    echo "
    export PERL5LIB=\${PERL5LIB}:${Software}/bioperl-live
    perl ${prinseq_script} -fastq ${output_quickAlign_human}/${sample}_1.fq -fastq2 ${output_quickAlign_human}/${sample}_2.fq  -min_qual_mean 20 -out_good ${output_qc}/Sample${sample}  -trim_qual_right     10 -trim_qual_left 10  -lc_method dust -lc_threshold 7 -no_qual_header -qual_noscale -graph_data ${output_qc}/Sample${sample}_quality.plot  -out_bad null                                                    
    " >> $script

    echo -e "#7) Submit blastn job against human reference"  >> $script
    echo "
    ${seqtk} seq -a ${output_qc}/Sample${sample}_1.fastq > ${output_qc}/Sample${sample}_1.fasta
    ${seqtk} seq -a ${output_qc}/Sample${sample}_2.fastq > ${output_qc}/Sample${sample}_2.fasta
    cat ${output_qc}/Sample${sample}_1.fasta ${output_qc}/Sample${sample}_2.fasta > ${output_qc}/${sample}_paired.fasta
    rm ${output_qc}/Sample${sample}_1.fasta  ${output_qc}/Sample${sample}_2.fasta
    " >> $script

    echo "
    awk 'BEGIN{RS=\">\"}NR>1{sub(\"\n\",\"\t\"); gsub(\"\n\",\"\"); print \$0}' ${output_qc}/${sample}_paired.fasta > ${output_qc}/NM_paired_${sample}.txt
    " >> $script

    echo "
    $blastn -db $db -query  ${output_qc}/${sample}_paired.fasta  -outfmt 6  -num_alignments 1   -evalue 1  -culling_limit 1 -num_threads 6  > ${output_blastn}/${sample}_paired.ncbiBLASTn
     " >> $script

    echo -e "#8) filter out blastn hits and keep filtered dataset"  >> $script
    echo "
    # Make fasta files for non overlapping reads (paired end)
    file1=${output_blastn}/${sample}_paired.ncbiBLASTn
    file1tmp=${output_blastn}/${sample}_paired_ncbi_tmp.txt
    file2=${output_qc}/NM_paired_${sample}.txt
    file2tmp=${output_qc}/NM_paired_${sample}_tmp.txt
    filtered=${output_blastn}/${sample}_paired_filtered.txt
    awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$3}' \${file1} | uniq > \${file1tmp}
    awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${file2} | uniq > \${file2tmp}
    echo \$file1tmp
    echo \$file2tmp
    awk  'BEGIN {FS=OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$2,\$3}' \${file1tmp} \${file2tmp} > \${filtered}                                                     
    wc -l \${filtered} > ${output_blastn}/NM_paired.number
    awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filtered} > ${output_blastn}/NM_paired_${sample}_filtered.fasta
    " >> $script

    # blast against univec
    echo "
    $blastn -db $dbunivec -query  ${output_blastn}/NM_paired_${sample}_filtered.fasta  -outfmt 6  -num_alignments 1   -evalue 1  -culling_limit 1 -num_threads 6  > ${output_univec_blastn}/${sample}_paired.ncbiBLASTn
    " >> $script

    echo -e "#9) filter out blastn hits and keep filtered dataset"  >> $script
    echo "
    ####### Make fasta files for non overlapping reads (paired end)
    file1=${output_univec_blastn}/${sample}_paired.ncbiBLASTn
    file1tmp=${output_univec_blastn}/${sample}_paired_ncbi_tmp.txt
    file2=${output_blastn}/${sample}_paired_filtered.txt
    file2tmp=${output_blastn}/NM_paired_${sample}_tmp.txt
    filtered=${output_univec_blastn}/${sample}_paired_filtered.txt
    awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$3}' \${file1} | uniq > \${file1tmp}
    awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${file2} | uniq > \${file2tmp}
    echo \$file1tmp
    echo \$file2tmp
    awk  'BEGIN {FS=OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$2,\$3}' \${file1tmp} \${file2tmp} > \${filtered}
    wc -l \${filtered} > ${output_univec_blastn}/NM_paired.number
    awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filtered} > ${output_univec_blastn}/NM_paired_${sample}_filtered.fasta
    " >> $script

    echo "
    #megaBLAST for nucleotides
    #$blastn -db $dbnucl -query ${output_univec_blastn}/NM_paired_${sample}_filtered.fasta -outfmt 6 -max_target_seqs 10 -max_hsps 1 -num_threads 6  >  ${output_nucleotide}/${sample}_megaBLAST.tab
    " >> $script

    # And now continue processing suitable for protein search
    # First create fastqs from blastn human filtered fasta to feed into bowtie2
    echo "
    filtered=${output_blastn}/${sample}_paired_filtered.txt
    awk  'BEGIN {FS=OFS=\"\\t\"} {print \$1}'  \${filtered}  > ${output_blastn}/ids.txt
    ${seqtk} subseq ${output_qc}/Sample${sample}_1.fastq ${output_blastn}/ids.txt > ${output_blastn}/${sample}_filtered_1.fastq
    ${seqtk} subseq ${output_qc}/Sample${sample}_2.fastq ${output_blastn}/ids.txt > ${output_blastn}/${sample}_filtered_2.fastq
    " >> $script

    echo "#9) Bowtie2 - rRNA"  >> $script
    echo "
    ${bowtie} -p 6 -x ${referencerRNA}/SILVA_128_SSU_LSU_UniVec -1 ${output_blastn}/${sample}_filtered_1.fastq -2 ${output_blastn}/${sample}_filtered_2.fastq -S ${output_quickAlign_rRNA}/${sample}.sam
    " >> $script

    echo -e "#4) Select pairs of reads that are both non mapping"  >> $script
    echo "
    ${samtools} view -bS  -o ${output_quickAlign_rRNA}/${sample}.bam ${output_quickAlign_rRNA}/${sample}.sam
    ${samtools} sort  ${output_quickAlign_rRNA}/${sample}.bam ${output_quickAlign_rRNA}/${sample}_sorted
    ${samtools} index  ${output_quickAlign_rRNA}/${sample}_sorted.bam
    ${samtools} view -u -f 12 -F 256 ${output_quickAlign_rRNA}/${sample}.bam > ${output_quickAlign_rRNA}/${sample}_unmapped.bam
    ${bedtools} bamtofastq -i  ${output_quickAlign_rRNA}/${sample}_unmapped.bam -fq ${output_quickAlign_rRNA}/${sample}_1.fq -fq2 ${output_quickAlign_rRNA}/${sample}_2.fq
    ${seqtk} seq -a ${output_quickAlign_rRNA}/${sample}_1.fq > ${output_quickAlign_rRNA}/${sample}_1.fasta
    ${seqtk} seq -a ${output_quickAlign_rRNA}/${sample}_2.fq > ${output_quickAlign_rRNA}/${sample}_2.fasta
    cat  ${output_quickAlign_rRNA}/${sample}_1.fasta  ${output_quickAlign_rRNA}/${sample}_2.fasta > ${output_quickAlign_rRNA}/${sample}_paired.fasta
    awk 'BEGIN{RS=\">\"}NR>1{sub(\"\n\",\"\t\"); gsub(\"\n\",\"\"); print \$0}' ${output_quickAlign_rRNA}/${sample}_paired.fasta > ${output_quickAlign_rRNA}/${sample}_paired_filtered.txt
    " >> $script

    # megaBLAST to filter out rRNA
    #blastn against ribosomal RNA (silva)
    echo -e "#9) Submit blastn job against rRNA reference"  >> $script
    echo "
    $blastn -db $dbrRNA -query ${output_quickAlign_rRNA}/${sample}_paired.fasta  -outfmt 6  -num_alignments 1   -evalue 0.1  -culling_limit 1 -num_threads 6  > ${output_rRNA}/${sample}_paired.ncbiBLASTn
    " >> $script

    echo -e "#10) filter out blastn hits and keep filtered dataset"  >> $script
    echo "
    # Make fasta files for non overlapping reads (paired end)
    file1=${output_rRNA}/${sample}_paired.ncbiBLASTn
    file1tmp=${output_rRNA}/${sample}_paired_ncbi_tmp.txt
    file2=${output_quickAlign_rRNA}/${sample}_paired_filtered.txt
    file2tmp=${output_quickAlign_rRNA}/${sample}_paired_filtered_tmp.txt
    filtered=${output_rRNA}/${sample}_paired_filtered_rRNA.txt
    awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$3}' \${file1} | uniq > \${file1tmp}
    awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${file2} | uniq > \${file2tmp}
    echo \$file1tmp
    echo \$file2tmp
    awk  'BEGIN {FS=OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$2,\$3}' \${file1tmp} \${file2tmp} > \${filtered}
    wc -l \${filtered} > ${output_rRNA}/NM_paired.number
    awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filtered} > ${output_rRNA}/NM_paired_${sample}_filtered_rRNA.fasta
    " >> $script

    echo -e "#nucleotide similarity"  >> $script
    echo "
    # megaBLAST for nucleotides
    $blastn -db $dbnucl -query ${output_rRNA}/NM_paired_${sample}_filtered_rRNA.fasta -outfmt 6 -max_target_seqs 10 -max_hsps 1 -num_threads 6  >  ${output_nucleotide}/${sample}_megaBLAST_norRNA.tab
    " >> $script

    echo -e "#11) Submit DIAMOND for protein"  >> $script

    echo "
    $diamond blastx -d $dbprot -q ${output_rRNA}/NM_paired_${sample}_filtered_rRNA.fasta   -f 6  -k 10 -p 6 -sensitive -e 1  > ${output_protein}/${sample}_diamond.tab
    " >> $script

    # echo -e "12) Remove large files that won't be reused i.e sam "
    echo "
    # remove files that take up too much storage space and won't be used downstream
    # from bowtie folders
    rm ${output_quickAlign_human}/${sample}.sam
    rm  ${output_quickAlign_rRNA}/${sample}.sam
    # from QC folder
    # rm ${output_qc}/NM_paired_${sample}.txt
    # rm ${output_qc}/NM_paired_${sample}_tmp.txt
    # from blastn folders
    # rm ${output_blastn}/${sample}_paired_ncbi_tmp.txt
     " >> $script

    # make reads stats summary
    echo -e "#13) Create reads stats summary"  >> $script

    echo "
    grep 'Processed reads:' ${output_trim}/${sample}_1.fastq.gz_trimming_report.txt | awk -F\":\" 'BEGIN {OFS=\"\t\"}{print \"Raw reads\", \$2*2}' > ${output_summary}/${sample}_raw.count
    grep 'Trimmed reads:' ${output_trim}/${sample}_1.fastq.gz_trimming_report.txt | awk -F\":\" 'BEGIN {OFS=\"\t\"}{print \"Trimmed reads\", \$2*2}' > ${output_summary}/${sample}_trimmed.count
    wc -l  ${output_quickAlign_human}/${sample}_1.fq |  awk  'BEGIN {OFS=\"\t\"}{ print \"Non human (bowtie2)\", (\$1/4)*2}' > ${output_summary}/${sample}_nonhuman.count
    wc -l  ${output_quickAlign_rRNA}/${sample}_1.fq |  awk 'BEGIN {OFS=\"\t\"}{ print \"Non rRNA (bowtie2)\", (\$1/4)*2}' > ${output_summary}/${sample}_nonrRNA.count
    wc -l  ${output_QC}/${sample}_paired.fasta | awk 'BEGIN {OFS=\"\t\"}{ print \"after QC\", $1/2}' > ${output_summary}/${sample}_QC.count
    awk '{s+=\$1} END {print s}' ${output_blastn}/NM_paired.number | awk -F\":\" 'BEGIN {OFS=\"\t\"}{ print \"non human (blastn)\", \$1}'  > ${output_summary}/${sample}_nonhuman_blastn.count
    awk '{s+=\$1} END {print s}' ${output_rRNA}/NM_paired.number | awk -F\":\" 'BEGIN {OFS=\"\t\"}{ print \"non rRNA (blastn)\", \$1}'  > ${output_summary}/${sample}_nonrRNA_blastn.count
    awk -F\"\t\" '{print \$1}' ${output_protein}/${sample}_diamond.tab | sort | uniq | wc -l | awk -F\" \" '{print \"Reads mapping to prot DB\", \$1}' > ${output_summary}/${sample}_maptoprot.count
    awk -F\"\t\" '{print \$1}' ${output_nucleotide}/${sample}_megaBLAST.tab | sort | uniq | wc -l | awk -F\" \" '{print \"Reads mapping to nucl DB\", \$1}' > ${output_summary}/${sample}_maptonucl.count
    cat ${output_summary}/${sample}_raw.count ${output_summary}/${sample}_trimmed.count ${output_summary}/${sample}_nonhuman.count ${output_summary}/${sample}_nonrRNA.count  ${output_summary}/${sample}_QC.count ${output_summary}/${sample}_nonhuman_blastn.count ${output_summary}/${sample}_nonrRNA_blastn.count ${output_summary}/${sample}_maptoprot.count ${output_summary}/${sample}_maptonucl.count  > ${output_summary}/${sample}_summary.count.new
    " >> $script

    #qsub  $script # Jrayner commented out qsub
fi
