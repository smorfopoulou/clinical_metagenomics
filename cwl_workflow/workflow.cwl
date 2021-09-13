class: Workflow
cwlVersion: v1.0
label: Encephalitis pipeline prototype for GOSH.
doc: This is a run of machineprofile null of pipeline Encephalitis_pipeline
id: Encephalitis_pipeline
requirements:
- class: InlineJavascriptRequirement
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
inputs:
  blastn_2__evalue:
    type: float?
  blastn__culling_limit:
    type: int?
  diamond__max_target_seqs:
    type: int?
  trim_galore__quality:
    type: int
  prinseqlite__trim_qual_left:
    type: int
  blastn_2__num_alignments:
    type: int?
  tar_tool_4__extractdirectory:
    type: string
  nuclLength:
    type: File
  Refseq_Blast_DB:
    type: File
  trim_galore__length:
    type: int
  Read2:
    type: File
  diamond__evalue:
    type: float
  blastn_3__evalue:
    type: float?
  accession2taxid:
    type: File
  Human_Bowtie_Reference:
    type: File
  tar_tool_2__extractdirectory:
    type: string
  tar_tool_5__extractdirectory:
    type: string
  prinseqlite__lc-threshold:
    type: int
  awk__program:
    type: string
  diamond__output_fmt:
    type: int
  blastn__num_alignments:
    type: int?
  RNA_Blast_DB:
    type: File
  RNA_Bowtie_Reference:
    type: File
  prinseqlite__lc-method:
    type: string
  Read1_Fastq:
    type: File
  blastn_4__max_hsps:
    type: int?
  blastn_4__max_target_seqs:
    type: int?
  blastn_3__num_alignments:
    type: int?
  blastn__evalue:
    type: float?
  Univec_Blast_DB:
    type: File
  Refseq_Diamond_DB:
    type: File
  blastn_3__culling_limit:
    type: int?
  names:
    type: File
  tar_tool__extractdirectory:
    type: string
  Human_Blast_DB:
    type: File
  blastn_2__culling_limit:
    type: int?
  prinseqlite__trim_qual_right:
    type: int
  tar_tool_6__extractdirectory:
    type: string
  tar_tool_3__extractdirectory:
    type: string
  prinseqlite__min_qual_mean:
    type: int
outputs:
  awk__result:
    outputSource: awk/result
    type: File
  metamix__presentSpecies_assignedReads:
    outputSource: metamix/presentSpecies_assignedReads
    type: File
  metamix__logLikelihood_traceplot_all:
    outputSource: metamix/logLikelihood_traceplot_all
    type: File
  blastn__result:
    outputSource: blastn/result
    type: File
  filter_blastn_2__count:
    outputSource: filter_blastn_2/count
    type: File
  filter_blastn_2__paired_filtered:
    outputSource: filter_blastn_2/paired_filtered
    type: File
  diamond__result:
    outputSource: diamond/result
    type: File
  prinseqlite__fastq_output_2:
    outputSource: prinseqlite/fastq_output_2
    type: File
  prinseqlite__plot:
    outputSource: prinseqlite/plot
    type: File
  prinseqlite__fastq_output_1:
    outputSource: prinseqlite/fastq_output_1
    type: File
  filter_blastn_2__filtered_blastn_hits_fasta:
    outputSource: filter_blastn_2/filtered_blastn_hits_fasta
    type: File
  metamix__histograms_cdf:
    outputSource: metamix/histograms_cdf
    type: File
  trim_galore__fastq_R2:
    outputSource: trim_galore/fastq_R2
    type: File
  metamix__allreads_classified:
    outputSource: metamix/allreads_classified
    type: File
  trim_galore__fastq_R1:
    outputSource: trim_galore/fastq_R1
    type: File
  encephalitis_summary_report_tool__report:
    outputSource: encephalitis_summary_report_tool/report
    type: File
  metamix__logLikelihood_traceplot_40:
    outputSource: metamix/logLikelihood_traceplot_40
    type: File
  blastn_2__result:
    outputSource: blastn_2/result
    type: File
steps:
  prinseqlite:
    run: prinseqlite.cwl
    in:
      trim_qual_right:
        source: prinseqlite__trim_qual_right
        default:
        - 10
      trim_qual_left:
        source: prinseqlite__trim_qual_left
        default:
        - 10
      fastq2: bam2fastq/fastq_r2
      lc-method:
        source: prinseqlite__lc-method
        default:
        - dust
      fastq1: bam2fastq/fastq_r1
      min_qual_mean:
        source: prinseqlite__min_qual_mean
        default:
        - 20
      lc-threshold:
        source: prinseqlite__lc-threshold
        default:
        - 7
    out:
    - fastq_output_1
    - plot
    - error_log
    - fastq_output_2
  tar_tool_6:
    run: untar_reference.cwl
    in:
      extractdirectory:
        source: tar_tool_6__extractdirectory
        default:
        - custom_nucleotide_refseq_nounplaced_July2019_March2020.fasta
      tarfile: Refseq_Blast_DB
    out:
    - extacted_tar_dir
    - error_log
  tar_tool_5:
    run: untar_reference.cwl
    in:
      extractdirectory:
        source: tar_tool_5__extractdirectory
        default:
        - SILVA_128_SSU_LSU_UniVec
      tarfile: RNA_Blast_DB
    out:
    - extacted_tar_dir
    - error_log
  metamix:
    run: metamix.cwl
    in:
      accession2taxid: accession2taxid
      taxonomy_names: names
      read_lengths: awk/result
      nuclLength: nuclLength
      blast_output: blastn_4/result
    out:
    - allreads_classified
    - histograms_cdf
    - logLikelihood_traceplot_40
    - logLikelihood_traceplot_all
    - presentSpecies_assignedReads
    - step2_species_annotated
  tar_tool_2:
    run: untar_reference.cwl
    in:
      extractdirectory:
        source: tar_tool_2__extractdirectory
        default:
        - human_GRCh38.p9_additions
      tarfile: Human_Blast_DB
    out:
    - extacted_tar_dir
    - error_log
  tar_tool_4:
    run: untar_reference.cwl
    in:
      extractdirectory:
        source: tar_tool_4__extractdirectory
        default:
        - SILVA_128_SSU_LSU_UniVec
      tarfile: RNA_Bowtie_Reference
    out:
    - extacted_tar_dir
    - error_log
  tar_tool_3:
    run: untar_reference.cwl
    in:
      extractdirectory:
        source: tar_tool_3__extractdirectory
        default:
        - UniVec
      tarfile: Univec_Blast_DB
    out:
    - extacted_tar_dir
    - error_log
  encephalitis_summary_report_tool:
    run: summary_stats.cwl
    in:
      NM_paired: filter_blastn/count
      diamond: diamond/result
      fastq_2: bam2fastq/fastq_r2
      paired_fasta: converFastqsToFastasAndConcat/fasta
      RNA_NM_paired: filter_blastn_3/count
      fastq_1: bam2fastq/fastq_r1
      trimming_fastq_1: trim_galore/report_1
      blast: blastn_4/result
    out:
    - report
    - error_log
    - summary_script
  tar_tool:
    run: untar_reference.cwl
    in:
      extractdirectory:
        source: tar_tool__extractdirectory
        default:
        - human_GRCh38.p9_addition
      tarfile: Human_Bowtie_Reference
    out:
    - extacted_tar_dir
    - error_log
  blastn_2:
    run: blastn.cwl
    in:
      evalue:
        source: blastn_2__evalue
        default:
        - 1.0
      culling_limit:
        source: blastn_2__culling_limit
        default:
        - 1
      query_fasta: filter_blastn/filtered_blastn_hits_fasta
      num_alignments:
        source: blastn_2__num_alignments
        default:
        - 1
      db: tar_tool_3/extacted_tar_dir
    out:
    - result
    - error_log
  filterMapReads_2:
    run: filter_mapped_reads.cwl
    in:
      bai: samtools_index_2/index
      bam: samtools_sort_2/sorted_bam
    out:
    - filtered_bam
    - error_log
  blastn_4:
    run: blastn.cwl
    in:
      max_hsps:
        source: blastn_4__max_hsps
        default:
        - 1
      max_target_seqs:
        source: blastn_4__max_target_seqs
        default:
        - 10
      query_fasta: filter_blastn_3/filtered_blastn_hits_fasta
      db: tar_tool_6/extacted_tar_dir
    out:
    - result
    - error_log
  trim_galore:
    run: trim_galore.cwl
    in:
      fastqs: Read1_Fastq
      fastq2: Read2
      length:
        source: trim_galore__length
        default:
        - 25
      quality:
        source: trim_galore__quality
        default:
        - 20
    out:
    - fastq_R1
    - fastq_R2
    - report_1
    - report_2
    - error_log
    - out_log
  blastn_3:
    run: blastn.cwl
    in:
      evalue:
        source: blastn_3__evalue
        default:
        - 0.1
      culling_limit:
        source: blastn_3__culling_limit
        default:
        - 1
      query_fasta: converFastqsToFastasAndConcat_2/fasta
      num_alignments:
        source: blastn_3__num_alignments
        default:
        - 1
      db: tar_tool_5/extacted_tar_dir
    out:
    - result
    - error_log
  samtools_sort:
    run: samtools_sort.cwl
    in:
      bam: sam2bam/bam
    out:
    - sorted_bam
    - error_log
  sam2bam_2:
    run: sam2bam.cwl
    in:
      sam: bowtie2_2/sam
    out:
    - bam
  filterMapReads:
    run: filter_mapped_reads.cwl
    in:
      bai: samtools_index/index
      bam: samtools_sort/sorted_bam
    out:
    - filtered_bam
    - error_log
  bam2fastq:
    run: bam2fastq.cwl
    in:
      bam: filterMapReads/filtered_bam
    out:
    - fastq_r1
    - fastq_r2
    - error_log
  converFastqsToFastasAndConcat_2:
    run: seqtk_seq.cwl
    in:
      fastq_2: bam2fastq_2/fastq_r2
      fastq_1: bam2fastq_2/fastq_r1
    out:
    - fasta
    - NM
    - error_log
    - script
  bowtie2_2:
    run: bowtie.cwl
    in:
      reference: tar_tool_4/extacted_tar_dir
      fastq_R1: seqtk_subseq/filtered_fastq1
      fastq_R2: seqtk_subseq/filtered_fastq2
    out:
    - sam
    - error_log
  bowtie2:
    run: bowtie.cwl
    in:
      reference: tar_tool/extacted_tar_dir
      fastq_R1: trim_galore/fastq_R1
      fastq_R2: trim_galore/fastq_R2
    out:
    - sam
    - error_log
  converFastqsToFastasAndConcat:
    run: seqtk_seq.cwl
    in:
      fastq_2: prinseqlite/fastq_output_2
      fastq_1: prinseqlite/fastq_output_1
    out:
    - fasta
    - NM
    - error_log
    - script
  samtools_index_2:
    run: samtools_index.cwl
    in:
      bam_sorted: samtools_sort_2/sorted_bam
    out:
    - index
    - error_log
  blastn:
    run: blastn.cwl
    in:
      evalue:
        source: blastn__evalue
        default:
        - 1.0
      culling_limit:
        source: blastn__culling_limit
        default:
        - 1
      query_fasta: converFastqsToFastasAndConcat/fasta
      num_alignments:
        source: blastn__num_alignments
        default:
        - 1
      db: tar_tool_2/extacted_tar_dir
    out:
    - result
    - error_log
  bam2fastq_2:
    run: bam2fastq.cwl
    in:
      bam: filterMapReads_2/filtered_bam
    out:
    - fastq_r1
    - fastq_r2
    - error_log
  sam2bam:
    run: sam2bam.cwl
    in:
      sam: bowtie2/sam
    out:
    - bam
  filter_blastn_3:
    run: filter_blastn.cwl
    in:
      blastn_result: blastn_3/result
      blast_query: converFastqsToFastasAndConcat_2/NM
    out:
    - filtered_blastn_hits_fasta
    - count
    - ids
    - error_log
    - filter_blastn_script
    - paired_filtered
  diamond:
    run: diamond.cwl
    in:
      output_fmt:
        source: diamond__output_fmt
        default:
        - 6
      evalue:
        source: diamond__evalue
        default:
        - 1.0
      max_target_seqs:
        source: diamond__max_target_seqs
        default:
        - 10
      query_fasta: filter_blastn_3/filtered_blastn_hits_fasta
      db: Refseq_Diamond_DB
    out:
    - result
    - error_log
  filter_blastn_2:
    run: filter_blastn.cwl
    in:
      blastn_result: blastn_2/result
      blast_query: filter_blastn/paired_filtered
    out:
    - filtered_blastn_hits_fasta
    - count
    - ids
    - error_log
    - filter_blastn_script
    - paired_filtered
  awk:
    run: awk.cwl
    in:
      target_files:
        linkMerge: merge_flattened
        source:
        - filter_blastn_3/paired_filtered
      program:
        source: awk__program
        default:
        - BEGIN{FS=OFS="\t"} {print $1, length($2)}
    out:
    - result
  seqtk_subseq:
    run: seqtk_subseq.cwl
    in:
      fastq_2: prinseqlite/fastq_output_2
      ids: filter_blastn/ids
      fastq_1: prinseqlite/fastq_output_1
    out:
    - filtered_fastq1
    - error_log
    - filtered_fastq2
  filter_blastn:
    run: filter_blastn.cwl
    in:
      blastn_result: blastn/result
      blast_query: converFastqsToFastasAndConcat/NM
    out:
    - filtered_blastn_hits_fasta
    - count
    - ids
    - error_log
    - filter_blastn_script
    - paired_filtered
  samtools_index:
    run: samtools_index.cwl
    in:
      bam_sorted: samtools_sort/sorted_bam
    out:
    - index
    - error_log
  samtools_sort_2:
    run: samtools_sort.cwl
    in:
      bam: sam2bam_2/bam
    out:
    - sorted_bam
    - error_log
