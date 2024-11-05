from typing import List, Dict, Tuple

import json
import os
import datetime
import click
from rich.console import Console

from functions import load_config, load_prediction_list, load_allele_sequences, load_b2m_sequence, make_filepath, create_combined_sequence


def write_log_file(log:Dict, filepath:str, testing:bool) -> None:
    """ 
    This function writes a log file to the specified filepath.

    Args:
        log (Dict): The log dictionary to write to the file.
        filepath (str): The filepath to write the log file to.
        testing (bool): Whether we are in testing mode or not.

    """
    if not testing:
        with open(filepath, 'w') as f:
            json.dump(log, f, indent=4)
    pass


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



@click.command()
@click.option("--environment", default='local', help="The name of the environment, can either be local or poc.")
@click.option("--structure_set", default='full', help="The structure set to use for the predictions, can either be full, partial or single.")
@click.option("--experiment_number", default=None, help="The experiment number to use for the predictions.")
@click.option("--gpu_number", default='all', help="The GPU number to use for the predictions.")
@click.option("--testing", default=None, help="Whether we are in testing mode or not.")

def run_predictions(environment, structure_set, experiment_number, gpu_number, testing):

    # First we'll load the configuration file for the chosen environment
    config = load_config(environment)
    print (config)

    console = Console()

    # We'll set the experiment a3m to None for now
    experiment_a3m = None

    # We'll check the input parameter for whether we're in testing mode or not
    if testing is not None:
        testing = True
    else:
        testing = False

    # This code only runs MSA experiments, so we'll check if an experiment number has been provided
    if experiment_number is None:
        print ("No experiment number provided, please provide an experiment number and try again.")
        exit()

    # We'll set the a3m filepath for the experiment, and if that experiment a3m file doesn't exist we'll exit
    experiment_a3m_filepath = f"inputs/experiments/experiment{experiment_number}.a3m"
    if not os.path.exists(experiment_a3m_filepath):
        print (f"Experiment a3m file does not exist at {experiment_a3m_filepath}")
        exit()


    # We'll print out the CLI options selected
    print ('CLI options selected:')
    print (f'Environment: {environment}')
    print (f'Structure set: {structure_set}')
    print (f'Experiment number: {experiment_number}')
    print (f'GPU number: {gpu_number}')
    print (f'Testing: {testing}')


    # We'll create the output folder structure if it doesn't exist
    a3m_tmp_filepath = 'outputs/tmp'
    if not os.path.exists(a3m_tmp_filepath):
        os.makedirs(a3m_tmp_filepath)

    experiments_folder = 'outputs/experiments'
    if not os.path.exists(experiments_folder):
        os.makedirs(experiments_folder)

    experiment_folder = f"outputs/experiments/{experiment_number}"
    if not os.path.exists(experiment_folder):
        os.makedirs(experiment_folder)


    # if the log file exists we'll load it, otherwise we'll create an empty dictionary
    experiment_log_filepath = f"{experiment_folder}/log.json"
    if os.path.exists(experiment_log_filepath):
        with open(experiment_log_filepath, 'r') as f:
            experiment_log = json.load(f)
    else:
        experiment_log = {}
        # and then we'll write the empty dictionary to the log file
        write_log_file(experiment_log, experiment_log_filepath, testing)

    # finally for the set up we'll load the prediction data
    structures_to_predict, allele_sequences, b2m_seq = load_prediction_data(config, structure_set)

    if config is not None and len(structures_to_predict) > 0:

        with console.status("Running predictions...", spinner="dots"):
            
            # now we'll iterate over the list of predictions and run the predictions
            for structure in structures_to_predict:

                # we'll load the experiment a3m file, we do this for each prediction as we need to create one per prediction with the concatenated sequence as the first sequence in the alignment
                with open(experiment_a3m_filepath, 'r') as f:
                    experiment_a3m = f.read()

                # first we'll generate a a3m file for the prediction if it doesn't exist
                pdb_code = structure['pdb_code']
                
                # we'll check if the predictions have already been run for this PDB code
                if pdb_code in experiment_log:
                    if experiment_log[pdb_code]['status'] == 'done':
                        console.print(f"[bold green]Predictions already exist for {pdb_code}[/bold green]")
                        continue
                else:
                    experiment_log[pdb_code] = {'status': 'preparing'}
                    write_log_file(experiment_log, experiment_log_filepath, testing)

                # we'll create the local and in container filepath for the a3m file
                local_a3m_filepath = f"{a3m_tmp_filepath}/{pdb_code}_{experiment_number}.a3m"
                docker_a3m_filepath = f"/work/{config['OUTPUT_FOLDER']}/tmp/{pdb_code}_{experiment_number}.a3m"

                # we'll create the a3m file if it doesn't exist
                if not os.path.exists(local_a3m_filepath):
                    with open(local_a3m_filepath, 'w') as f:
                        combined_sequence = create_combined_sequence(allele_sequences, structure, b2m_seq, 274)
                        prediction_string = f">101\t102\t103\n{combined_sequence}"
                        this_a3m = experiment_a3m.replace('###', prediction_string)

                        f.write(this_a3m)
                
                # we'll create the in container filepath for the output folder and the local output folder
                item_path = f"{config['OUTPUT_FOLDER']}/experiments/{experiment_number}/{pdb_code}"
                local_output_folder = f"{config['PROJECT_FOLDER']}/{item_path}"
                docker_output_folder = f"/work/{item_path}"

                # we'll set the predictions_exist flag to False
                predictions_exist = False

                # we'll check if the predictions have already been run for this PDB code
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
                    write_log_file(experiment_log, experiment_log_filepath, testing)

                    if gpu_number == 'all':
                        gpu_field = "--gpus=all"
                    else:
                        gpu_field = f"--gpus=\"device={gpu_number}\""
                    colabfold_command = f"docker run --user $(id -u) -ti --rm {gpu_field} -v {config['AF_WEIGHTS_FOLDER']}:/cache:rw -v $(pwd):/work:rw {config['CONTAINER_IMAGE']} colabfold_batch {config['COLABFOLD_OPTIONS']} {docker_a3m_filepath} {docker_output_folder}"
                    
                    if testing:
                        console.print(colabfold_command)
                    else:
                        os.system(colabfold_command)

                        experiment_log[pdb_code]['status'] = 'done'
                        end_time = datetime.datetime.now()
                        experiment_log[pdb_code]['end_time'] = end_time.isoformat()
                        experiment_log[pdb_code]['elapsed_time'] = (end_time - start_time).total_seconds()
                        write_log_file(experiment_log, experiment_log_filepath, testing)

                else:
                    console.print(f"[bold green]Predictions already exist for {pdb_code}[/bold green]")
                    
    else:
        console.print("[bold red]Cannot run. There was an error loading the configuration, please check you have filled in the config file.[/bold red]")
    pass




if __name__ == "__main__":
    run_predictions()

