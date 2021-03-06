## Creating SDG workspaces

SDG can be used to explor GFA graphs created by any assembler, alongside raw data, and/or to create its own graphs. For all of these, a workspace needs to be created first, and data needs to mapped into the graph.

On this page we present code and examples of how to do precisely that. It can be done from the command line or thorugh the python interface.



### Raw reads to datastores

Before you can include them in a WorkSpace, your raw reads need to be put into datastores. The easiest way to do this is thorugh the `sdg-datastore` CLI tool.

To create the two datastores for the [e. coli test dataset](../datasets/datasets.md#E.-coli-paired-end-and-PacBio), simply run:

```shell
sdg-datastore make -t paired -o ecoli_pe -1 ../ecoli_pe_r1.fastq -2 ../ecoli_pe_r2.fastq
sdg-datastore make -t long -o ecoli_pb -L ../ecoli_pb_all.fastq
```

This will create the files `ecoli_pe.prseq` with the paired end data and `ecoli_pb.loseq` with the PacBio data.

### SDG Workspace from an existing GFA file

A SequenceDistanceGraph can be created by reading the graph topology, the sequences and the link distances from a GFA or GFA2 file. Because of the SDG representation, the following restrictions must be met:

* Sequences must be given just once (not RC).
* Overlaps must represent perfect matches.
* S records with length but no sequence in GFA1 are considered gaps and intepreted as positive distances.

ABySS and w2rap-contigger output graphs that can be read into SDG, and in general most short-read assemblers should produce output that is compatible with SDG.

As an example, we start by generating a `k=71` assembly of our [e. coli test dataset](../datasets/datasets.md#E.-coli-paired-end-and-PacBio) with abyss:

```shell
abyss-pe in='../ecoli_pe_r1.fastq ../ecoli_pe_r2.fastq' k=71 name=ecoli_abyss graph=gfa2
```

Now we will create an SDG workspace from the final ABySS file, add the paired and long reads, map them, and create  and add a KmerCounter with a count for the PE reads. This is a typical way to import a dataset to work on this problem.

```shell
sdg-workspace make -g ecoli_abyss-8.gfa2 -p ecoli_pe.prseq -L ecoli_pb.loseq -o sdg_ecoli_abyss
sdg-kmercounter make -w sdg_ecoli_abyss.sdgws -o sdg_ecoli_abyss_kcounts -d ecoli_pe.prseq -n PE
sdg-workspace add_counter -w sdg_ecoli_abyss.sdgws -n kcounts -c sdg_ecoli_abyss_kcounts.sdgkc -o sdg_ecoli_abyss
sdg-mapper -w sdg_ecoli_abyss.sdgws -o sdg_ecoli_abyss

```

The last two commands, when adding the KmerCounter to the WorkSpace and mapping the reads, use an output name that is the same of the input workspace. This effectively replaces the file of the old workspace with the new one. Since SDG by design does not update WorkSpaces on disk, this is the appropiate way to keep a single file when building the WorkSpace in stages.