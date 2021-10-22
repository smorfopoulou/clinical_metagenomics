
    library(metaMix)
    step1<-generative.prob.nucl(blast.output.file="/output/1/nucleotide/1_megaBLAST_norRNA.tab",  outDir="./",  blast.default=TRUE, read.length.file="/output/1/metamix_nucleotide_norRNA/read_lengths.tab", contig.weight.file=1,  accession.taxon.file="//DBs/refseq_nucleotide_28March2020/custom_nucleotide_refseq_nounplaced_July2019_March2020_nucl_gb.accession2taxid", genomeLength="//DBs/refseq_nucleotide_28March2020//custom_nucleotide_refseq_nounplaced_July2019_March2020_nuclLength.tab")
step2<-reduce.space(step1="step1.RData", taxon.name.map="//DBs/refseq_nucleotide_28March2020/taxonomy/names.dmp")    
step3<-parallel.temper.nucl(step2="step2.RData",  readSupport=10)
step4<-bayes.model.aver(step2="step2.RData", step3="step3.RData",  taxon.name.map="///DBs/refseq_nucleotide_28March2020/taxonomy/names.dmp")
     
