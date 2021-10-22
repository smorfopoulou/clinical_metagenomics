
#!/bin/bash
#$ -S /bin/bash
#$ -o /output/1/metamix_nucleotide_norRNA/cluster/out
#$ -e /1/metamix_nucleotide_norRNA/cluster/error
#$ -wd /output/1/metamix_nucleotide_norRNA 
#$ -pe smp 6
#$ -l tmem=7.1G,h_vmem=7.1G
#$ -l h_rt=72:00:00
#$ -V
#$ -R y
source /share/apps/source_files/metaMix.source
mpirun --mca btl_tcp_if_include 128.41.97.0/21 -np 1 /share/apps/genomics/metaMix/R-3.5.2/bin/R --slave CMD BATCH --no-save --no-restore  Rsubmit.R step.out

