#!/usr/bin/env cwl-runner
# view the workflow with https://view.commonwl.org/workflows

cwlVersion: v1.0
class: Workflow

outputs
  KPI_CO2_emission:
    type: float
    doc: "C02 emission of the global structure"
  KPI_t_removal_of_formwork:
    type: float
    doc: "time at which the formwork can be removed (stress under dead loads)"
  KPI_load_bearing_capacity_28:
    type: float
    doc: "Ratio between the maximum_stress and the strength after 28d at the most critical position"

inputs:
  mix_wz:
    type: float
    default: 1.0

steps:

  mix_design_performance:
    run: mix_design_performance_prediction.cwl
    in:
      domain_size: mix_wz
    out: [  hydration_A,
            stiffness_evolution_H,
            strength_evolution_S,
            shrinkage_ABC,
            CO2_per_m3,
            density,
            max_t,
            delta_t,
            geometry_A,
            load_A
         ]

  structural_simulation:
    run: structural_simulation.cwl
    in:
      domain_size: domain_size
    out: [KPI_CO2_emission, KPI_t_removal_of_formwork, KPI_load_bearing_capacity_28]

