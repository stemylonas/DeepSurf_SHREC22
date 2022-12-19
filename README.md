# DeepSurf - SHREC22
This repo contains a version of DeepSurf used in a SHREC 2022 track (http://shrec.ge.imati.cnr.it/shrec22_protein/) and aims at the replicability of the results presented in the corresponding publication ('SHREC 2022: Protein-ligand binding site recognition'). For the original DeepSurf method or the instalaltion process, please refer to the https://github.com/stemylonas/DeepSurf repo. 

Usage example
---------------

To replicate the DeepSurf results of 'SHREC 2022: Protein-ligand binding site recognition', download the corresponding dataset and execute the following command, where the test directory refers to the folder with the testing data.

```
python multi_predict.py -inp test_directory -mp model_path -o output_directory --f 1
```

