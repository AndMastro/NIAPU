# NIAPU
NIAPU: network-informed adaptive positive-unlabelled learning for disease genes identification

NIAPU is formed by two main components: the computation of the NeDBIT (Network diffusion and biology-informed topological) features and the usage of APU (Adaptive Positive-Unlabelled label propagation).

## How to use

* Step 1: **NeDBIT feature computation**

Use the file ```preprocess.pl``` which takes as input 1) ```ppi_file``` in the format **gene1** **gene2** containing the PPI links and 2) seed_genes_scores in the format **gene** **score** and outputs a) **out_links** which contains PPI links as ordinal numbers (for the subsequent computations) and **out_genes**, contating all the PPI genes (score 0 for non seed genes).

```
preprocess.pl ppi_file seed_genes_scores out_links out_genes
```

Then, after compiling nedbit_features_calculator.c, execute it as

```
nedbit_features_calculator out_links out_genes nedbit_features
```

which outputs **nedbit_features** file with the newly computed features.

* Step 2: **Adaptive PU labelling**

Use the provided ```APU.R``` file taking as input **nedbit_features** to obtain a gene ranking.

* Step 3: ML classification

Once you obtained the ranking, you can assign pseudo-labels to unlabbeled genes from it (details are in the paper below). We provide as an example a python script that reproduces the choices made in the paper and uses three ML algorithm (MLP, RF and SVM) to perform classification. However, once can decide to use their favorite classifiers.

