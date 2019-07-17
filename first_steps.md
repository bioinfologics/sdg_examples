# Workspace

The workspace is the foundation of bsg, a workspace contains everything you need to assemble a genome. To create a workspace from scratch, first import bsg from the installation directory

```python
import pysdg as SDG
```

then simply instanciate the `WorkSpace()` class. 

```python
ws = SDG.WorkSpace()
```

This is the the object where everything is going to be stored, the different components of the workspace can be accessed using dot notation.

The main components of the workspace are the graph (ws.sg), the counted kmers object (ws.kci) and the datastores and mappers (ws.)

![image-20181109112347970](/Users/ggarcia/Library/Application Support/typora-user-images/image-20181109112347970.png)

a workspace can be saved to disk and loaded from disk, this allows to preserve the sessions. See (functions to dump and restore to disk). In bsg the processes that are applied to a workspace and it's component live in memory, so if you want to persist the work that you have done in your graph you'll need to save a copy of the workspace to the disk, this can be done using.

```python
>>> ws.dump_to_disk("./persisted_workspece.bsgws")
```

and the session can be restored by loading the workspace instead of rebuilding it from scratch every time

```
## Create a new Workspace
>>> restored_ws = bsg.WorkSpace()
>>> restored_ws.load_from_disk("./persisted_workspece.bsgws")
```

Now restored_ws contains all the same elements as ws. 

Please note that some of the components that are derived from the data ( like graph kmer indexes) are not persisted in the session and need to be recomputed after loading the workspace from disk.

#### ws.sg

Next we need to add a graph as the base for the rest of the information, bsg as the name implies is graph based, so the graph is the 

in bsg sequences are stored in their canonical form in the nodes. Links store the topological relation between nodes:

![image-20181107213601519](/Users/ggarcia/Library/Application Support/typora-user-images/image-20181107213601519.png)


To import a graph in an empty workspace you can load it directly from a gfa v1.0 graph:

```python
>>> ws.sg.load_from_gfa('./graph.gfa')
```

Once the graph is loaded all the informatino can be accessed accessing the graph in the workspace via dot notation. The nodes are stored inside the graph in a vector of Node objects, nodes start from index 1, index 0 is reserved for the NULL node. The same is valid for the links, links[0] is reserved for the null link.

```
>>> print("Number of nodes: %s" %(len(ws.sg.nodes)))
>>> print("Number of links: %s" %(sum([len(ws.sg.links[n]) for n in range(1, len(ws.sg.nodes))]), ))
Number of nodes: 964121
Number of links: 2826446
```


Information in the nodes and links can also be accessed using dot notation ( for more information about the methods and attributes check the [Node](https://bioinfologics.github.io/bsg/class_node.html) and [Link](https://bioinfologics.github.io/bsg/class_link.html) documentation), for example if you need to acess the sequence of some nodes:

```python
>>> for node in ws.sg.nodes[:10]:
>>>     print(node.sequence)

GATTCTCAATTAACAAAATTAAAAATAAAATATATTCAAATTAAAATAAAATACTGGCTCCTGTTTCTCACTCTGTTCCTCCGATCTTATTTTAATTAATTAAAAAAAATAATTAATTCTTAATTCTAATTAAAGGGGGAATTCAAGAATTCCCTTAATTCAAGAAGAATTCCCCTTTGAATTCCCCAAGGAGAATATTA
AAAGTTAACTGTTAACTTTGTTCAAACAAAGAATAGCTTGCTAATCCCCCCCCCCTATTTCTTACCCAACTCTGTTGCAACATCATCGATGATTATGGCTTGGAAATAGAGTGTTCTCATTGTCATACCTCATATTCATCTTCATCTTCATGCCATTTGGTACCCCGTCTTAACGGTTTATCCGTATTCTGTCAAGTAGGATACTCTATTAGAGGAGCGCTTCCTATTTAATCATCCACTCATTGATTTAGTTTCGAT
CTGCAATTCACCCTGTTTTCCCATTGTTGGGCGACCCAATATAAACAAAATATTAACATTAACATCCTGTGAATGTTATTAAATAATAAGTAATAAATAATATAATCTTCATTAGAATGTATCTACTTAATGTTAGGTTCTTAAAAGCAAGCACTTAAACTATTAAAGTAACACAGGGATGGTTACCCTACTATTAGTTA
AGTATTAAAAATACTAATTATAATTAAATGAAGTGTGATCCGTTAATAAATAAATATATTAAAATGAGGGTGTTATTCCGTTCAACCAAATTGAGATGAAGTGAGATGAGAATGAGATTAGTCAGAAATCACCAAATATCAACCATTTAACCCCCCTTTCAACATTAACATCCTTAAATTATTTAATTTAAGGATTATAA
CCCCCCTTGCTTTCGCTTCCCTTCCAACTGCTCTTCACTCACTCACTTAATCAACGTTGATAAGTACAATTAGAATTAATGTGTTTGGTTCAACCATTCCCATCAAAATTACCAACCAACCAACCAACCTACCAACCAACCAACCAACCAACCAACCATCCTCAATTGGATCGAGAAGCGATAACTGCTTTAAATCAGAAATACTCTATTTCTTCAAGACAGCTGGAATA
ACTCTCCCCAAGTAGTACTCCGTTAAAAATACTCATCATGTCGTGACTGTGTGGACCTATTCTCATATCCCCTCTCATTCGCTTCGATGTCCGCTCTGATGTTTGGTAAACCATAATAAGTAGTAGGATGTTATTAGAATACTTAATAGGTCGCTGCATGTGAGGCGGTAAGATCTACTACCTCATAGTGACATCTGTCC
CCTCTAGCACAACTAATAGTAATTATGAGCTAGGTAGTTACTCGATAACGGTAAGAGTGACAGGTATGACTGAGGGGGTTGAAGGACCGTTTCCGAGGAGGTATAAACAGTAAGGGTATGAGCGAGAACGTTAAGGAAGCGTTAAGGAGTAGGGGTGGAGGTAGTGGGGGGTTGAATAAGATTTAAGACAGAGAAGTTATAGTATTAGAATCGTAGGTCTTATTTAAGGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTAGTGTGTTTGTAAAATATGGGATAAGGATAGCGATATTATTAGAGATAAGGTTGACACGGTGAGCGAGATAAGGGGTTGTGAGGGTTGTTGAAGTGAGA
CCTCTCCCCAAGTAGTACTCCGTTAAAAATACTCATCATGTCGTGACTGTGTGGACCTATTCTCATATCCCCTCTCATTCGCTTCGATGTCCGCTCTGATGTTTGGTAAACCATAATAAGTAGTAGGATGTTATTAGAATACTTAATAGGTCGCTGCATGTGAGGCGGTAAGATCTACTACCTCATAGTGACATCTGTCA
GATAGCATGAATGCTGAAAGTGTAACACAAAAATATGGTCCTTGGCTATTGTATATATATTATCAAGTCTTCCATGACCTTGCCCTGAAGCCAAACTGACTGACTGACTGACTTACCCTGAAGTCCAACTGACTGTGACTCACGTGAATCTCATAGGTAAACGATTGAAACGTTATAAAAAAAAATATAATCCCAGCCCAGAA
```


Or if you need to get the 101th node in the collection and check some properties you have to acess the 101th position in the nodes collection of the graph and check if it's sotred in the canonical orientation or not:

```
>>> node101 = ws.sg.nodes[101]
>>> print ("Sequence: " %(node101.sequence, ))
AGTAGGGGGATATAAAGATA ... GTGTTGAAACGCACGCGAGG

>>> print ("Is canonical: " %(node101.is_canonical(), ))
True
```

Most of the collections in the graph are indexed according to the node collection, so the 101th position in the links vector correspondes with the links of the 101th node.

```
>>> print("Number of links in node 101: %s" %(len(ws.sg.links[101]), ))
Number of links in node 101: 7

>>> print(ws.sg.links[101])
(<-101 -> 192099>, <-101 -> 299295>, <-101 -> 305649>, <-101 -> 648660>, <101 -> 392194>, <101 -> 410676>, <101 -> 865355>)
```


in most cases acessing the links/nodes directly is not necesarry because the graph object has functions to navigate the graph, for example for the 101th node we can get the forward facing links:

```
>>> print("Fowrard links:")
>>> for fwl in ws.sg.get_fw_links(101):
>>>     print(fwl)
Fowrard links:
-101 -> 192099
-101 -> 299295
-101 -> 305649
-101 -> 648660
```

Or the backwards facing links:

```
>>> print("Backwards links:")
>>> for bwl in ws.sg.get_bw_links(101):
>>>     print(bwl)
Backwards links:
101 -> 392194
101 -> 410676
101 -> 865355
```


The same navigation functions are available to traverse the graph using nodes only, so you can check which nodes are in the neighbourhood of the 101th node using:

```
>>> ## Both fw and bw nodes at the same time (neighbours)
>>> for nn in ws.sg.get_neighbour_nodes(101):
>>>     print(nn)
299295
192099
648660
865355
392194
305649
410676
```


An check which nodes are forward from 101th or backwards from 101th:

```
>>> ## Get forward nodes for the 101th node
>>> print ("Forward nodes: ")
>>> for fwn in ws.sg.get_fw_nodes(101):
>>>     print(fwn)
>>> print("")
>>> print ("Backwards nodes: ")
>>> for bwn in ws.sg.get_bw_nodes(101):
>>>     print(bwn)
Forward nodes:
192099
299295
305649
648660

Backwards nodes:
392194
410676
865355
```

This is just a small subset of the available functinos, for a complete reference check the documentation ( https://bioinfologics.github.io/bsg/ )

## Short reads (kci)

BSG can store the kmer count for the grph Kmers, the count of the graph kmers in the graph are stored in the graph, but also can store the kmer count of the graph kmers in other datasets. For example one important kmer count that is usually counted is the kmer count of each graph kmer in the illumina pcr free data that was used to create the contigs graph. This count is directly related with the number of copies of a particular kmer in the original genome. By comparing the number of times that a kmer appears in the illumina pe dataset en the number of times that that same kmers apperas in the produced graph reveals if the produced assembly captured the correct copy number for the assembled sequence. we call that relation KCI (kmer compression index).

Kmer counts in the Workspace are stored as KmerCount object vectors over the kmers existing in the graph. First we need to create a kmer index for the graph, this will take the sequences of the nodes, transform that into kmers and store the kmers in a graph kmers hash, at the same time will create a kmer count of the kmers in the graph.
This can be done as follows:

```
>>> ws.kci.index_graph()
```
this will popuate the `ws.kci.graph_kmers` collection of `KmerCount` this is needed to store the count of the kmers in the graph and to count the other kmer sources as well. Each KmerCount stores a kmer `KmerCount.kmer` and a count `KmerCount.count` of that kmer in the graph.

```
>>> ws.kci.graph_kmers[1000].kmer
19459
>>> ws.kci.graph_kmers[1000].count
19
```

Kmer counts for additinoal sources are stored in the `ws.kci.read_counts` vector. The first dimension of this vector is the each of the counted datasets, the second dimension is the the count for all the kmers of the graph in that particular dataset, the kmers keep the same order as in the graph, so the kmer 1001th in the graph Khmers vector is the same as the 1001th kmer in the read_counts vector but with the corresponding count. 

To add a new count to the `ws`, first add an empty container for the new count

```
>>> ws.kci.start_new_count()
```

and then count directly from the fastq file, this will accumulate the count to the last available container, in this case the new empty one.

```
>>> strVec = bsg.vectorString(["./pe_10M_R1.fastq","./pe_10M_R2.fastq",])
>>> ws.kci.add_counts_from_file(strVec)
```

- Accessing the kmer coverage

The count of the kmers in the graph and the kmer collections can be accessed by quering the collection of kmers in the graph and in the reads.

```
>>> kmer_position = 100
>>> print("Kmer sequence: %s" %(ws.kci.graph_kmers[kmer_position].kmer))
>>> print("Freq in the graph: %s" %(ws.kci.graph_kmers[kmer_position].count))
>>> print("Freq in the reads: %s" %(ws.kci.read_counts[0][kmer_position]))

Kmer sequence: 142
Freq in the graph: 3
Freq in the reads: 3
```

The kmer count can also be recovered for entire nodes using `ws.kci.compute_node_coverage_profil(<sequence>, , <collection_index>)` this will return the kmer count for each of the kmers in the provided sequence. The function returns 3 vectors of length `sequence_length-K+1`, the first vector is the kmer count for the kmers in the `reads_count[collection_index]` the second vector is the unique kmers profile, that is the kmers that are unique in the graph will have a 1 and the ones that are not unique qill have 0, and the 3rd collection is graph_kmer profile, is the count of the kmer in the graph.

For example take the kmer coverage profile of the node 7538 from the graph and count the kmers.

```
>>> nd = ws.sg.nodes[7538]
>>> kmers = [x for x in ws.kci.compute_node_coverage_profile(nd.sequence, 0)]

>>> plt.figure(figsize=(20, 5))
>>> plt.plot(kmers[0], 'r')
>>> plt.plot(kmers[1], 'g')
>>> plt.plot(kmers[2], 'b')
>>> plt.legend(["Reads", "Unique", "Graph"])
```

![image-20181121181048236](/Users/ggarcia/Library/Application Support/typora-user-images/image-20181121181048236.png)

- Calculate KCI

KCI is the kmer count in the reads and divided by the frequency of the first peak so the expressed values now express the expected copy number for the kmers. An average KCI can be calculated for each node by combining the kcis of all nodes.

First compute all the stats for the read kmer collection, this will calculate the mode of the unique content distirbution. The obtained mode can be checked acessing the internal value, it's a good idea to check this value before calculating the KCI for all the nodes.

```
>>> ws.kci.compute_compression_stats()
>>> ws.kci.uniq_mode
2
```

With this information the KCI for the all the nodes can be calculated

```
>>> ws.kci.compute_all_nodes_kci()
```

This finction will calculate the average KCI for all the nodes in the graph and store the results in `ws.kci.nodes_depth`, Only the kmers with graph frequency < 10 are used in the calculation, the max frequency can be changed.

Now the collection `ws.sg.nodes_depth` contains all the KCI values for the nodes in the graph in the same order as in the nodes collections.

To get the distribution of node KCI values a histogram can be created.

```
>>> plt.hist([x for x in ws.kci.nodes_depth], bins=50)
```

![image-20181121231946478](/Users/ggarcia/Library/Application Support/typora-user-images/image-20181121231946478.png)

Nodes with `kci~1` are nodes that corresponde to the heterozygous part of the genome while othre frequencies corresponds to other parts.

## Short reads (mp)

- create datastore
- add datastore to ws
- map reads
- access mappings

## Linked reads

The reads need to be included in the workspace as datastores. To create a datastore we need to create a datastore first and then import that datastore in the workspace and the last step is to map the reads to the graph.

```bash
## Create datastore Object
lrds = bsg.LinkedReadsDatastore()

## Import the reads to the datastore
lrds.build_from_fastq('./subssample_R1.fastq', './subssample_R2.fastq', '10xds.10xds', bsg.LinkedReadsFormat_UCDavis)

## Append the datastore to the current workspace
ws.getLinkedReadDatastores().append(lrds)
```

now the datastore is stored in the linked reads collection in the first position, more datastores can be added if necesarry

```
ws.linked_read_datastores[0]
```

Next step is to map the reads to the graph, the reads are going to be mapped using a unique Khmers approach, a read is considered mapping if and only if all the reads in the read map uniquely to the same contig. every other case will be discarded by the mapper and the read is going to be marked as non-mapped.

To map the reads first we need to create an index of unique Khmers (this can be done at k-63 as well), then add a read mapper to the workspace and finally call remap_all to remal all read collections in the workspace (this functino will trigger a remap in all datatypes).

```
ws.create_index()
ws.create_63mer_index()
lrm = bsg.LinkedReadMapper(ws.sg, ws.linked_read_datastores[0], ws.uniqueKmerIndex, ws.unique63merIndex)
ws.getLinkedReadMappers().append(lrm)
ws.remap_all()
```

The result of this step is the transformation functions between the rad data and the graph that we are going to need to execute the untangling steps. This collections are 

```
ws.linked_read_mappers[0].reads_in_node[1]
```

will give you the reads of the first mapped linked library in the node id 1.

You can also get some qc statistics for the datastore like tag ocupancy (interesting for icing datasets)

```bash
lrds.dump_tag_occupancy_histogram('tag_ocupancy_histogram.hist')
```

Also you can get get the tags for a particular read or all reads with a particular tag

```
lrds.get_read_tag(100001)
for r in lrds.get_tag_reads(122385086):
    print(ws.r)
```

## Long reads

- create datastore

  

  ```bash
  ## Create datastore Object
  lrds = bsg.LinkedReadsDatastore()
  
  ## Import the reads to the datastore
  lrds.build_from_fastq('./subssample_R1.fastq', './subssample_R2.fastq', '10xds.10xds', bsg.LinkedReadsFormat_UCDavis)
  
  ## Append the datastore to the current workspace
  ws.getLinkedReadDatastores().append(lrds)
  ```

  

  

- add datastore to ws

- map reads

- access mappings

## Basic untangling or genome resolution

- select by kci and show the idea of a backbone selectino for haplotype phasing
- connect the backbones using some data
- Produce lines, save the graph output, modify the lines to show that can be done.

