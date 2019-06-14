#!/usr/bin/env cwl-runner
class: CommandLineTool
id: run-immclassifier
label: run-immclassifier
cwlVersion: v1.0

requirements:
  - class: DockerRequirement
    dockerPull: guoxindi/r-garnett
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.synapse_config)
        entryname: .synapseConfig

baseCommand: [Rscript, /usr/local/bin/run-garnett.R]

inputs:
  synapse_config:
    type: File
  input_path:
    type: File
    inputBinding:
      position: 1
      prefix: --input
  output_path:
    type: string
    inputBinding:
      position: 2
      prefix: --output
  marker_path:
    type: File
    inputBinding:
      position: 3
      prefix: --marker

outputs:
  predictions:
    type: File
    outputBinding:
      glob: $(inputs.output_path)
