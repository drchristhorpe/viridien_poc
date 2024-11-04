from typing import List, Dict, Tuple

import json
import os
import datetime
import click
from rich.console import Console

from functions import load_config, load_prediction_list, load_allele_sequences, load_b2m_sequence, make_filepath, create_combined_sequence

def write_log_file(log:Dict, filepath:str):
    with open(filepath, 'w') as f:
        json.dump(log, f, indent=4)

def get_prediction_path(experiment_number: str) -> str:
    predictions_path = None
    if experiment_number is None:
        pre
    else:
        # or we'll use the experiment number to create the output folder
        predictions_path = f"experiments/{experiment_number}"
    return predictions_path


def load_prediction_data(config:Dict, structure_set:str) -> Tuple[List[Dict], Dict, str]:
    # first we'll load the list of predictions to run from the CSV file 
    structures_to_predict = load_prediction_list(config)

    # then we'll load the allele sequences for the alleles in the list of predictions
    allele_sequences = load_allele_sequences(structures_to_predict, config)

    # next we'll load the canonical sequence of the B2M gene
    b2m_seq = load_b2m_sequence(config)

    if structure_set == 'partial':
        # we'll filter the list of predictions to only include the first 10
        structures_to_predict = structures_to_predict[:10]

    elif structure_set == 'single':
        # we'll filter the list of predictions to only include the first one
        structures_to_predict = structures_to_predict[:1]
    
    return structures_to_predict, allele_sequences, b2m_seq


def generate_colabfold_command():
    pass


@click.command()
@click.option("--environment", default='local', help="The name of the environment, can either be local or poc.")
@click.option("--structure_set", default='full', help="The structure set to use for the predictions, can either be full, partial or single.")
@click.option("--experiment_number", default=None, help="The experiment number to use for the predictions.")
@click.option("--gpu_number", default='all', help="The GPU number to use for the predictions.")
@click.option("--gpu_count", default=None, help="The number of GPUs to use for the predictions.")
def run_predictions(environment, structure_set, experiment_number, gpu_number, gpu_count):

    config = load_config(environment)

    print ('CLI options selected:')
    print (f'Environment: {environment}')
    print (f'Structure set: {structure_set}')
    print (f'Experiment number: {experiment_number}')
    print (f'GPU number: {gpu_number}')
    print (f'GPU count: {gpu_count}')

    experiment_a3m = None
    a3m_tmp_filepath = 'outputs/tmp'

    if not os.path.exists(a3m_tmp_filepath):
        os.makedirs(a3m_tmp_filepath)


    if experiment_number is None:
        exit()

    else:

        
        experiment_a3m_filepath = f"inputs/experiments/experiment{experiment_number}.a3m"
        
        if not os.path.exists(experiment_a3m_filepath):
            print (f"Experiment a3m file does not exist at {experiment_a3m_filepath}")
            exit()

        experiment_folder = f"outputs/experiments/{experiment_number}"
        if not os.path.exists(experiment_folder):
            os.makedirs(experiment_folder)
        
        # we'll load the JSON log file for the experiment
        experiment_log_filepath = f"{experiment_folder}/log.json"

        if os.path.exists(experiment_log_filepath):
            with open(experiment_log_filepath, 'r') as f:
                experiment_log = json.load(f)
        else:
            experiment_log = {}
            write_log_file(experiment_log, experiment_log_filepath)

        
        
    
    console = Console()
    if config is not None:

        with console.status("Running predictions...", spinner="dots"):
            
            structures_to_predict, allele_sequences, b2m_seq = load_prediction_data(config, structure_set)


            # now we'll iterate over the list of predictions and run the predictions
            for structure in structures_to_predict:

                with open(experiment_a3m_filepath, 'r') as f:
                    experiment_a3m = f.read()

                # first we'll generate a a3m file for the prediction if it doesn't exist, these won't change, so we'll only generate them once
                pdb_code = structure['pdb_code']
                
                if pdb_code in experiment_log:
                    if experiment_log[pdb_code]['status'] == 'done':
                        console.print(f"[bold green]Predictions already exist for {pdb_code}[/bold green]")
                        continue
                else:
                    experiment_log[pdb_code] = {'status': 'preparing'}
                    write_log_file(experiment_log, experiment_log_filepath)

                # we'll create the local and in container filepath for the a3m file
                local_a3m_filepath = f"{a3m_tmp_filepath}/{pdb_code}_{experiment_number}.a3m"
                docker_a3m_filepath = f"/work/{config['OUTPUT_FOLDER']}/tmp/{pdb_code}_{experiment_number}.a3m"

                if not os.path.exists(local_a3m_filepath):
                    with open(local_a3m_filepath, 'w') as f:
                        combined_sequence = create_combined_sequence(allele_sequences, structure, b2m_seq, 274)
                        prediction_string = f">101\t102\t103\n{combined_sequence}"
                        this_a3m = experiment_a3m.replace('###', prediction_string)

                        f.write(this_a3m)
                
                # we'll check if the predictions have already been run for this PDB code
                item_path = f"{config['OUTPUT_FOLDER']}/experiments/{experiment_number}/{pdb_code}"
                print (f'Item path: {item_path}')
                local_output_folder = f"{config['PROJECT_FOLDER']}/{item_path}"
                print (f'Local output folder: {local_output_folder}')
                docker_output_folder = f"/work/{item_path}"
                print (f'Docker output folder: {docker_output_folder}')
                predictions_exist = False

                if os.path.exists(local_output_folder):
                    # each colabfold run will generate 26 files, so if we have 26 files we can assume the predictions are done for that PDB code
                    if len(os.listdir(local_output_folder)) == 26:
                        predictions_exist = True
                else:
                    os.makedirs(local_output_folder)

                # if the predictions don't exist we'll run them
                if not predictions_exist:
                    experiment_log[pdb_code]['status'] = 'running'
                    start_time = datetime.datetime.now()
                    experiment_log[pdb_code]['start_time'] = start_time.isoformat()
                    write_log_file(experiment_log, experiment_log_filepath)

                    colabfold_command = f"docker run --user $(id -u) -ti --rm --gpus={gpu_number} -v {config['AF_WEIGHTS_FOLDER']}:/cache:rw -v $(pwd):/work:rw {config['CONTAINER_IMAGE']} colabfold_batch {config['COLABFOLD_OPTIONS']} {docker_a3m_filepath} {docker_output_folder}"
                    os.system(colabfold_command)

                    experiment_log[pdb_code]['status'] = 'done'
                    end_time = datetime.datetime.now()
                    experiment_log[pdb_code]['end_time'] = end_time.isoformat()
                    experiment_log[pdb_code]['elapsed_time'] = (end_time - start_time).total_seconds()
                    write_log_file(experiment_log, experiment_log_filepath)

                    #console.print(colabfold_command)
                else:
                    console.print(f"[bold green]Predictions already exist for {pdb_code}[/bold green]")
                    
    else:
        console.print("[bold red]Cannot run. There was an error loading the configuration, please check you have filled in the config file.[/bold red]")
    pass




if __name__ == "__main__":
    run_predictions()

