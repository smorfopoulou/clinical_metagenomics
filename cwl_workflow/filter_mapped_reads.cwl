#!/usr/bin/env cwl-runner

# (Re)generated by BlueBee Platform

$namespaces:
  bb: http://bluebee.com/cwl/
  ilmn-tes: http://platform.illumina.com/rdf/iap/
cwlVersion: cwl:v1.0
class: CommandLineTool
bb:toolVersion: '1'
id: filterMappedReads
requirements:
- class: ShellCommandRequirement
hints:
- class: ResourceRequirement
  ilmn-tes:resources:
    tier: standard
    type: standardHiCpu
    size: small
- class: DockerRequirement
  dockerPull: clinical_metagenomics_preprocessing:latest
label: filterMapReads
stderr: stderr.txt
inputs:
  bam:
    type: File
    secondaryFiles: # ToDo: remove secondary file (I think this is left over from old tool version?)
    - ''
    inputBinding:
      position: 100
  bai:
    type: File
outputs:
  filtered_bam:
    type: File
    outputBinding:
      glob:
      - '*.filtered.bam'
  error_log:
    type: stderr
arguments:
- position: 0
  valueFrom: -u
- position: 3
  prefix: -f
  valueFrom: '12'
- position: 4
  prefix: -F
  valueFrom: '256'
- position: 1
  prefix: -@
  valueFrom: \$(nproc --all)
  shellQuote: false  
- position: 2
  prefix: -o
  valueFrom: $(inputs.bam.nameroot).filtered.bam
baseCommand:
- samtools
- view