<p align="center">
  <img src="docs/NIAPU_LOGO.svg" alt="NIAPU logo" width=30%>
</p>

![C](https://img.shields.io/badge/C-00599C?style=for-the-badge&logo=c&logoColor=white) ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![Perl](https://img.shields.io/badge/Perl-39457E?style=for-the-badge&logo=perl&logoColor=white)

# NIAPU

This is the official implementation for [**NIAPU: network-informed adaptive positive-unlabeled learning for disease genes identification**](https://doi.org/10.1093/bioinformatics/btac848).

NIAPU is formed by two main components: the computation of the NeDBIT (Network diffusion and biology-informed topological) features and the usage of APU (Adaptive Positive-Unlabelled label propagation).

## How to use

* Step 1: **NeDBIT feature computation**

Use the file ```preprocess.pl``` (you should have Perl installed) which takes as input 1) ```ppi_file``` in the format **gene1** **gene2** containing the PPI links (a Biogrid-derived file is provided in the data folder) and 2) seed_genes_scores in the format **gene** **score** and outputs a) **out_links** which contains PPI links as ordinal numbers (for the subsequent computations) and **out_genes**, contating all the PPI genes (score 0 for non seed genes).

```
perl preprocess.pl ppi_file seed_genes_scores out_links out_genes
```

Then, after compiling nedbit_features_calculator.c, execute it as

```
nedbit_features_calculator out_links out_genes nedbit_features
```

which outputs **nedbit_features** file with the newly computed features.

* Step 2: **Adaptive PU labelling**

After having compiled ```apu_label_propagation.c```, execute it. It takes as input **nedbit_features**, a 0/1 Boolean value indicating the presence of an header in the features file, the name of the output gene ranking and two float values indicating the quantile threshold for the removal of weak links and for the Reliable Negative computation, respectively.

```
./apu_label_propagation nedbit_features HEADER_PRESENCE output_gene_ranking 0.05 0.2 
```

Of course, one is free to use their own features with the APU labelling system, provided the correct format.

* Step 3: **ML classification**

Once you obtained the ranking along with the pseudo-labels, you can use those to train ML algorithms. We provide as an example a python script ```classification.py``` that reproduces the choices made in the paper and uses three ML algorithm (MLP, RF and SVM) to perform classification. This script is ready-to-use even without performing the steps above if one wants to use the provided disease file and features. It will look for APU scores in the folder ```../data/APU_scores/```. Filename format is specified in the file. However, one can decide to use their favorite classifiers and data. 

The feature computation code and the label propagation algorithm are written in C for efficiency reasons; we improved the previous R implementation obtaining a great speedup of the execution time.

## Citation

If you use our work, please cite our paper 😊

Paola Stolfi, Andrea Mastropietro, Giuseppe Pasculli, Paolo Tieri, Davide Vergni, NIAPU: Network-Informed Adaptive Positive-Unlabeled learning for disease gene identification, Bioinformatics, 2023;, btac848, https://doi.org/10.1093/bioinformatics/btac848
