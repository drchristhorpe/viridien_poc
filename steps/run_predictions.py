from typing import List, Dict

import json
import os
import datetime
import click
from rich.console import Console

from functions import load_config, load_prediction_list, load_allele_sequences, load_b2m_sequence, make_filepath, create_tmp_fasta_file






@click.command()
@click.option("--environment", default='local', help="The name of the environment, can either be local or poc.")
def run_predictions(environment):
    config = load_config(environment)

    console = Console()
    if config is not None:
        
        
        # we'll create a timestamp to use for the output folder, in the final PoC this will be a unique identifier for the conditions of the run
        timestamp = datetime.datetime.now().strftime("%Y%m%d")
        
        # we'll create the output folder structure if it doesn't exist
        predictions_path = make_filepath(config, 'output', 'predictions', timestamp)
        if not os.path.exists(predictions_path):
            os.makedirs(predictions_path)


        with console.status("Running predictions...", spinner="dots"):
            # first we'll load the list of predictions to run from the CSV file
            structures_to_predict = load_prediction_list(config)

            # then we'll load the allele sequences for the alleles in the list of predictions
            allele_sequences = load_allele_sequences(structures_to_predict, config)
            
            # next we'll load the canonical sequence of the B2M gene
            b2m_seq = load_b2m_sequence(config)

            # now we'll iterate over the list of predictions and run the predictions
            for structure in structures_to_predict:

                # first we'll generate a fasta file for the prediction if it doesn't exist, these won't change, so we'll only generate them once
                pdb_code = structure['pdb_code']
                tmp_fasta_file = make_filepath(config, 'output', 'fasta', f"{pdb_code}.fasta")
                if not os.path.exists(tmp_fasta_file):
                    fasta_str = create_tmp_fasta_file(allele_sequences, structure, b2m_seq, tmp_fasta_file)

                
                # we'll create the in container filepath for the fasta file
                docker_fasta_file = f"/work/{config['OUTPUT_FOLDER']}/fasta/{pdb_code}.fasta"
                
                # we'll check if the predictions have already been run for this PDB code
                item_path = f"{config['OUTPUT_FOLDER']}/predictions/{timestamp}/{pdb_code}"
                local_output_folder = f"{config['PROJECT_FOLDER']}/{item_path}"
                docker_output_folder = f"/work/{item_path}"
                predictions_exist = False

                if os.path.exists(local_output_folder):
                    # each colabfold run will generate 26 files, so if we have 26 files we can assume the predictions are done for that PDB code
                    if len(os.listdir(local_output_folder)) == 26:
                        predictions_exist = True

                # if the predictions don't exist we'll run them
                if not predictions_exist:
                    colabfold_command = f"docker run --user $(id -u) -ti --rm --gpus=all -v {config['AF_WEIGHTS_FOLDER']}:/cache:rw -v $(pwd):/work:rw {config['CONTAINER_IMAGE']} colabfold_batch {config['COLABFOLD_OPTIONS']} {docker_fasta_file} {docker_output_folder}"
                    os.system(colabfold_command)
                    
    else:
        console.print("[bold red]Cannot run. There was an error loading the configuration, please check you have filled in the config file.[/bold red]")
    pass




if __name__ == "__main__":
    run_predictions()

