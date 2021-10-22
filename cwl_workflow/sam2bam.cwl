#!/usr/bin/env cwl-runner

# (Re)generated by BlueBee Platform

$namespaces:
  bb: http://bluebee.com/cwl/
  ilmn-tes: http://platform.illumina.com/rdf/iap/
cwlVersion: cwl:v1.0
class: CommandLineTool
bb:toolVersion: '1'
id: samtools view
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
label: sam2bam
inputs:
  sam:
    type: File
    streamable: true
    inputBinding:
      position: 100
outputs:
  bam:
    type: File
    streamable: true
    outputBinding:
      glob:
      - '*.bam'
arguments:
- position: 3
  prefix: -o
  valueFrom: $(inputs.sam.nameroot).bam
- position: 2
  prefix: -@
  valueFrom: \$(nproc --all)
  shellQuote: false    
baseCommand:
- samtools
- view
