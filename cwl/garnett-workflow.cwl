class: Workflow
label: imm-class-workflow
id: imm-class-workflow
cwlVersion: v1.0

inputs:
  input_id:
    type: string
  output_name:
    type: string
  marker_id:
    type: string
  synapse_config:
    type: File
  folder_id:
    type: string

outputs:
  prediction:
    type: File

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement

steps:
  get-input-file:
    run: https://raw.githubusercontent.com/Sage-Bionetworks/synapse-client-cwl-tools/master/synapse-get-tool.cwl
    in:
      synapseid: input_id
      synapse_config: synapse_config
    out: 
      input_file
  get-marker-file:
    run: https://raw.githubusercontent.com/Sage-Bionetworks/synapse-client-cwl-tools/master/synapse-get-tool.cwl
    in:
      synapseid: marker_id
      synapse_config: synapse_config
    out:
      marker_file
  run-garnett:
    run: run-garnett.cwl
    in:
      synapse_config: synapse_config
      input_path: get-input-file/input_file
      output_path: output_name
      marker_path: get-marker_file/marker_file
    out:
      output_file
  store-output-file:
    run: https://raw.githubusercontent.com/Sage-Bionetworks/synapse-client-cwl-tools/master/synapse-store-tool.cwl
    in: 
      synapse_config: synapse_config
      file_to_store: run-garnett/output_file
      parentid: folder_id
    out: 
      []
