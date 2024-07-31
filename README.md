# Viridien peptide:MHC AlphaFold prediction Proof of Concept
Proof of Concept computing experiment around AlphaFold tuning for peptide:MHC prediction using custom MSAs


## Setting up the environment

It's been designed to be as self-contained and lightweight as possible. You'll need to create a python virtual environment and install the packages contained in the 'requirements.txt' file

## Installing relevant files

The Colabfold container relies upon the latest AlphaFold weights. These can be obtained by following the step in the instructions for running [Colabfold in Docker](https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker)

## Updating the 'config.toml' file

Once you have downloaded the AlphaFold weights, you'll need to edit the 'poc.toml' file to provide the path to the cache folder containing the weights and the project folder containing the code. You'll also need to provide the name of the container image. The 'local.toml' gives an example of local development. You could also change the values in this file to run the software locally. 

## Running a simple set of predictions

For running locally just type

'python steps/run_predictions.py'

To run in the PoC environment

'python steps/run_predictions.py --environment poc'