#!/usr/bin/env cwl-runner
# view the workflow with https://view.commonwl.org/workflows

cwlVersion: v1.0
class: Workflow

outputs:
  KPI_CO2_emission:
    type: float
    doc: "C02 emission of the global structure"

inputs:
  mix_wz:
    type: float
    default: 1.0


steps:
  mix_design_performance_prediction:
    run: mix_design_performance_prediction.cwl
    in:
      mix_wz: mix_wz
    out: [CO2]

  structural_simulation:
    run: structural_simulation.cwl
    in :
            CO2: mix_design_performance_prediction/CO2
    out: [KPI_CO2_emission]

