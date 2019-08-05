# SDG example datasets

This is a list of all the datasets used on the examples.

## *E. coli* paired-end and PacBio

This is a simple dataset containing a subsample of a test *E. coli* K12 2x300bp MiSeq run, and the PacBio data from [Reducing assembly complexity of microbial genomes with single-molecule sequencing (Koren et al.)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r101)

This dataset can be downloaded from [https://opendata.earlham.ac.uk/opendata/data/sdg_datasets/ecoli/](https://opendata.earlham.ac.uk/opendata/data/sdg_datasets/ecoli/).

## PseudoSeq.jl trio dataset

The synthetic genome creation and sequencing package [Pseudoseq.jl](https://github.com/bioinfologics/Pseudoseq.jl) was used to create two parental genomes.

Chromosomes 4 and 5 of the reference genome of the yeast strain S288C were used as a template to create a diploid, 2 chromosome genome for each parent. The proportion of heterozygous sites in each parent's genome was set to ~1\%.

Crossing over and recombination between the homologous pairs of chromosomes was simulated by dividing them into an number of chunks of equal base pairs in length, and then swapping those chunks at random. Then, for each homologous pair of chromosomes, the child inherited one chromosome from the first parent at random, and one chromosome from the second parent at random.\\

These three simulated genomes were then pseudo-sequenced using Pseudoseq.jl, generating simulated read datasets. Simulated paired end reads were generated for each genome, using an average fragment length of 700bp and a read length of 250bp, and an expected coverage of 70x with error rate was set to .1\%.

This dataset can be downloaded from https://opendata.earlham.ac.uk/opendata/data/sdg_datasets/trio/.

