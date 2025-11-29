from typing import List, Dict

import os




# colabfold_batch custommsa.a3m output_directory

colabfold_command = f"docker run --user $(id -u) -ti --rm --gpus=all -v /home/chris/science/cache:/cache:rw -v $(pwd):/work:rw ghcr.io/sokrypton/colabfold:1.5.5-cuda11.8.0 colabfold_batch --num-recycle 3 --random-seed 42 --amber --use-gpu-relax work/outputs/experiments/experiment9.a3m work/outputs/test"

print (colabfold_command)
#colabfold_command = f"docker run --user $(id -u) -ti --rm --gpus=all -v {config['AF_WEIGHTS_FOLDER']}:/cache:rw -v $(pwd):/work:rw {config['CONTAINER_IMAGE']} colabfold_batch {config['COLABFOLD_OPTIONS']} {docker_fasta_file} {docker_output_folder}"
os.system(colabfold_command)