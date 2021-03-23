from pathlib import Path
from doit.task import clean_targets
import sys

#DOIT_CONFIG = {"default_tasks": ["test"]}
PYTHON_EXE = sys.executable


def task_generate_artifical_data():
    """
    Create the artifical data and store them in csv file.
    """
    script = Path(__file__).parents[0] / 'generate_artifical_data.py'
    
    metadata_files = Path(__file__).parent.glob('*_meta.yaml')
    data_files = Path(__file__).parent.glob('*_data.csv')

    return {
        'actions': [f'{PYTHON_EXE} {script}'],
        'file_dep': [script, *metadata_files],
        "targets": [*data_files],
        "verbosity": 2,  # show stdout
        #"setup": ['tbd'],
        "clean": [clean_targets]
    }

def task_diffusion_model():
    """
    Simulation of diffusion.
    """
    script = Path(__file__).parents[0] / 'diffusion_model.py'
    
    modeldata_files = Path(__file__).parent.glob('*_model.yaml')
    model_files = Path(__file__).parent.glob('*_model.csv')

    return {
        'actions': [f'{PYTHON_EXE} {script}'],
        'file_dep': [script, *modeldata_files],
        "targets": [*model_files],
        "verbosity": 2,  # show stdout
        #"setup": ['tbd'],
        "clean": [clean_targets]
    }