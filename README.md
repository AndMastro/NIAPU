# NIAPU
NIAPU: network-informed adaptive positive-unlabelled learning for disease genes identification


NIAPU is formed by two main components: the computation of the NeDBIT (Network diffusion and biology-informed topological) features and the usage of APU (Adaptive Positive-Unlabelled label propagation).

In order to use NIAPU, perform the following steps.

* Step 1: **NeDBIT feature computation**

Use the file ```preprocess.pl``` which takes as input 1) ```ppi_file``` in the format **gene1** **gene2** containing the PPI links and 2) seed_genes_scores in the format **gene** **score** and outputs a) **out_links** which contains PPI links as ordinal numbers (for the subsequent computations) and **out_genes**, contating all the PPI genes (score 0 for non seed genes).

```
preprocess.pl ppi_file seed_genes_scores out_links out_genes
```
