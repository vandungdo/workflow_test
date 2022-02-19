#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  Prediction of mix design properties
baseCommand: ["python structural.py"]
arguments: ["hydration_A", "$(inputs.param_hydration)",
            "stiffness_evolution_H", "$(inputs.param_stiffness_evolution_H)",
            "strength_evolution_S", "$(inputs.param_strength_evolution_S)",
            "shrinkage_ABC", "$(inputs.param_shrinkage_ABC)",
            "CO2_per_m3", "$(inputs.CO2_per_m3)",
            "density", "$(inputs.density)",
            "max_t", "$(inputs.max_t)",
            "delta_t", "$(inputs.max_t)",
            "geometry_A", "$(inputs.geometry_A)",
            "load_A", "$(inputs.load_A)",
            ]

inputs:
  mix_wz:
    type: float
    doc: "water cement ratio of the mix"

outputs:
  hydration_A:
    type: float
    doc: "parameter A of the hydration model (possibly add more here)"
  stiffness_evolution_H:
    type: float
    doc: "parameter H of the stiffness evolution model (possibly add more here)"
  strength_evolution_S:
    type: float
    doc: "parameter S of the strength evolution model (possibly add more here)"
  shrinkage_ABC:
    type: float
    doc: "parameter ABC of the shrinkage model (possibly add more here)"
  CO2_per_m3:
    type: float
    doc: "CO2 emission per m³ of concrete (make sure to upscale from mortar to concrete)"
  density:
    type: float
    doc: "density of concrete"
  max_t:
    type: float
    doc: "maximum time for the computation"
  delta_t:
    type: float
    doc: "initial timestep for the computation"
  geometry_A:
    type: float
    doc: "parameter A of the geometry (e.g. length/height/thickness)"
  load_A:
    type: float
    doc: "point load in the center (define coordinate systems)"