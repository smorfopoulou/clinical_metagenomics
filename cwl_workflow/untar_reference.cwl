#!/usr/bin/env cwl-runner

# (Re)generated by BlueBee Platform

$namespaces:
  bb: http://bluebee.com/cwl/
  ilmn-tes: http://platform.illumina.com/rdf/iap/
cwlVersion: cwl:v1.0
class: CommandLineTool
bb:toolVersion: '1'
id: untar_reference
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: clinical_metagenomics_preprocessing:latest
- class: ResourceRequirement
  ilmn-tes:resources:
    tier: standard
    type: standardHiCpu
    size: small
label: untar_reference
stderr: stderr.txt
inputs:
  tarfile:
    type: File
    inputBinding:
      position: 0
      prefix: -xf
  extractdirectory:
    type: string
    inputBinding:
      position: 1
outputs:
  extacted_tar_dir:
    type: Directory
    outputBinding:
      glob:
      - $(inputs.extractdirectory)
  error_log:
    type: stderr
baseCommand:
- tar
