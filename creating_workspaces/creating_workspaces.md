# Creating SDG workspaces

SDG can be used to explor GFA graphs created by any assembler, alongside raw data, and/or to create its own graphs. For all of these, a workspace needs to be created first, and data needs to mapped into the graph.

On this page we present code and examples of how to do precisely that. It can be done from the command line or thorugh the python interface.



## SDG Workspace from an existing GFA file and PE+LMP data

A SequenceDistanceGraph can be created by reading the graph topology, the sequences and the link distances from a GFA or GFA2 file. Because of the SDG representation, the following restrictions must be met:

* Sequences must be given just once (not RC).
* Overlaps must represent perfect matches.
* S records with length but no sequence in GFA1 are considered gaps and intepreted as positive distances.

ABySS and w2rap-contigger output graphs that can be read into SDG, and in general most short-read assemblers should produce output that is compatible with SDG.

As an example, we start by generating a `k=71` assembly of our [e. coli test dataset](../datasets/datasets.md#E.-coli-paired-end-and-PacBio) with abyss:



```shell
abyss-pe in='../ecoli_pe_r1.fastq ../ecoli_pe_r2.fastq' k=71 name=ecoli_abyss graph=gfa2
```



