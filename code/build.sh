# This script allows for the creation of environments on different platforms ####

# Setup ####
conda remove -p ./env --all

# Create a conda environment from the exported configuration file ####
conda env create -f config.yml -p ./env
