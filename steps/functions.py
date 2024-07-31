from typing import Dict, List

import toml
import csv
import json


def make_filepath(config:Dict, in_or_out:str, foldername:str, filename:str) -> str:
    """
    This function will generate a filepath based on the input parameters.

    Args:
        config (Dict): A dictionary details of input/ouput paths and project location
        in_or_out (str): A string to specify if the file is an input or output file
        foldername (str): The name of the folder where the file is located
        filename (str): The name of the file

    Returns:
        filepath (str): The path to the file
    """
    if in_or_out == 'input':
        filepath = f"{config['PROJECT_FOLDER']}/{config['INPUT_FOLDER']}/{foldername}/{filename}"
    elif in_or_out == 'output':
        filepath = f"{config['PROJECT_FOLDER']}/{config['OUTPUT_FOLDER']}/{foldername}/{filename}"
    else:
        raise ValueError("in_or_out must be either 'input' or 'output'")
        filepath = None
    return filepath



def deslugify_allele_slug(allele_slug:str) -> str:
    """
    This function takes an allele slug and returns a deslugified allele number.

    Args:
        allele_slug (str): The allele slug to deslugify e.g. hla_a_24_02

    Returns:
        str: The deslugified allele number e.g. HLA-A*24:02
    """
    # split the allele slug on underscores
    allele_number_components = allele_slug.split("_")
    # join the components with the correct characters
    allele_number = f"{allele_number_components[0]}-{allele_number_components[1]}*{allele_number_components[2]}:{allele_number_components[3]}"
    # return the allele number
    return allele_number.upper()


def load_config() -> Dict:
    """
    This function will load the configuration file from the config.toml file.

    Args:
        None

    Returns:
        config (Dict): A dictionary with the following keys
            tmp_filepath (str): The path to the temporary directory to save the fasta file
            mhcflurry_models (str): The path to the mhcflurry models
            netmhcpan_path (str): The path to the netmhcpan executable
    """
    with open("config.toml", "r") as config_file:
        config = toml.load(config_file)
    return config


def load_prediction_list(config:Dict) -> List:
    """
    This function will load the prediction list from the hla_class_i.csv file which is either hand generated or from a datasette query.
    
    Args:
        config (Dict): A dictionary details of input/ouput paths and project location
    
    Returns:
        structures_to_predict (List): A list of dictionaries with the following keys
            allele_slug (str): The allele_slug to predict
            locus (str): The locus of the allele
            peptide (str): The peptide to predict
            pdb_code (str): The PDB code to predict
            resolution (str): The resolution of the PDB structure
    """
    with open(make_filepath(config, 'input', 'complexes', 'hla_class_i.csv'), "r") as allele_list_file:
        structures_to_predict = list(csv.DictReader(allele_list_file))
    return structures_to_predict


def load_allele_sequences(predictions_to_run:List[Dict], config:Dict) -> dict:
    """
    This function will generate a dictionary of allele sequences for the alleles in the predictions_to_run list.

    Args:
        predictions_to_run (List[Dict]): A list of dictionaries with the following keys
            allele_slug (str): The allele_slug to predict
            allele_number (str): The allele number to predict
            peptide (str): The peptide to predict
            pdb_code (str): The PDB code to predict

    Returns:
        allele_seq_dict (dict): A dictionary with the allele as the key and the sequence as the value
    """
    allele_list = sorted(list(set([prediction['allele_slug'] for prediction in predictions_to_run])))
    loci = ['hla_a', 'hla_b', 'hla_c']
    allele_seq_dict = {}
    locus_dict = {}
    for locus in loci:
        with open(make_filepath(config, 'input', 'sequences', f"{locus}.json"), "r") as locus_file:
            locus_dict[locus] = json.load(locus_file)
    for allele in allele_list:
        allele_locus = allele[:5]
        allele_seq_dict[allele] = locus_dict[allele_locus][allele]['canonical_sequence']
    return allele_seq_dict


def load_b2m_sequence(config:Dict) -> str:
    """
    This function will load the canonical sequence of the B2M gene from the human_b2m.json file.

    Args:
        e

    Returns:
        b2m_seq (str): The canonical sequence of the B2M gene
    """
    with open(make_filepath(config, 'input', 'sequences', 'human_b2m.json'), "r") as b2m_file:
        b2m_seq = json.load(b2m_file)['canonical_sequence']
    return b2m_seq


def create_tmp_fasta_file(allele_sequences:Dict, prediction:Dict, b2m_seq:str, filename:str) -> str:
    """
    This function will create a temporary fasta file with the concatenated sequence for the prediction in the format:
        mhci_allele_sequence:b2m_sequence:peptide
    Args:
        allele_sequences (Dict): A dictionary with the allele as the key and the sequence as the value
        prediction (Dict): A dictionary with the following keys
            allele_slug (str): The allele_slug to predict
            allele_number (str): The allele number to predict
            peptide (str): The peptide to predict
            pdb_code (str): The PDB code to predict
        b2m_seq (str): The canonical sequence of the B2M gene
        filename (str): The name of the file to save the fasta file
    Returns:
        fasta_file (str): The string of the FASTA file
    """
    prediction_sequence = f"{allele_sequences[prediction['allele_slug']].replace('-','')}:{b2m_seq}:{prediction['peptide_sequence']}"
    fasta_file = f">{prediction['pdb_code']} | {deslugify_allele_slug(prediction['allele_slug'])}:{prediction['peptide_sequence']}\n{prediction_sequence}\n"
    with open(filename, "w") as fasta:
        fasta.write(fasta_file)
    return fasta_file