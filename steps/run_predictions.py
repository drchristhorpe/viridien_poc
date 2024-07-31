from typing import List, Dict

import json
import os
import datetime


from functions import load_config, load_prediction_list, load_allele_sequences, load_b2m_sequence, make_filepath, create_tmp_fasta_file

from rich.console import Console

console = Console()



def run_predictions(config:Dict):
    timestamp = datetime.datetime.now().strftime("%Y%m%d%H")
    
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

        for structure in structures_to_predict:

            pdb_code = structure['pdb_code']

            tmp_fasta_file = make_filepath(config, 'output', 'fasta', f"{pdb_code}.fasta")

            if not os.path.exists(tmp_fasta_file):
                fasta_str = create_tmp_fasta_file(allele_sequences, structure, b2m_seq, tmp_fasta_file)

            docker_fasta_file = f"/work/{config['OUTPUT_FOLDER']}/fasta/{pdb_code}.fasta"
            docker_output_folder = f"/work/{config['OUTPUT_FOLDER']}/predictions/{timestamp}/{pdb_code}"



            
            colabfold_command = f"docker run --user $(id -u) -ti --rm --gpus=all -v {config['AF_WEIGHTS_FOLDER']}:/cache:rw -v $(pwd):/work:rw {config['CONTAINER_URL']} colabfold_batch {config['COLABFOLD_OPTIONS']} {docker_fasta_file} {docker_output_folder}"
            print (colabfold_command)
            os.system(colabfold_command)
            break
    pass




if __name__ == "__main__":
    config = load_config()
    run_predictions(config)

