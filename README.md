# DeepSurf - SHREC22
This repo contains a version of DeepSurf used in a SHREC 2022 track (http://shrec.ge.imati.cnr.it/shrec22_protein/) and aims at the replicability of the results presented in the corresponding publication ('SHREC 2022: Protein-ligand binding site recognition'). For the original DeepSurf method, please refer to the https://github.com/stemylonas/DeepSurf repo.

Installation
---------------

1) Python 3 and CUDA 9 are required 
2) Install DMS from http://www.cgl.ucsf.edu/Overview/software.html
3) Install openbabel (version 2.4.1 originally used)
4) Download trained models from https://drive.google.com/file/d/1nIBoD3_5nuMqgRGx4G1OHZwLsiUjb7JG/view?usp=sharing
5) Install python dependencies (requirements.txt)
6) Execute 'lds/compile' to compile the custom LDS-module. If you have g++>=5, add -D_GLIBCXX_USE_CXX11_ABI=0 to the g++ commands.


Usage example
---------------

To replicate the DeepSurf results f 'SHREC 2022: Protein-ligand binding site recognition', download the corresponding dataset and execute the following command, where the tesr directory refers to the fodler with the testing data.

```
python multi_predict.py -inp test_directory -mp model_path -o output_directory --f 1
```

