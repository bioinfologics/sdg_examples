# SDG first steps walkthrough

## Workspace

The workspace is the container of a project, contains everything needed work with a  genome. 

To create a workspace from scratch, first import pysdg.

```python
import pysdg as SDG
```

instanciate a `WorkSpace()` and save the instance in a variable.

```python
ws = SDG.WorkSpace()
```

This is the the object where everything is going to be stored, the different components of the workspace can be accessed using dot notation.

The main components of the workspace are the graph (ws.sdg), the datastores with their corresponding mappers and the KmerCounters (ws.kci) and most frequently all this data is accessed usgin NodeViews().

![image-20190730114736552](./workspace-diagram.png)

a workspace can be saved to disk and loaded from disk, this allows to preserve the sessions. See (functions to dump and restore to disk). In bsg the processes that are applied to a workspace and it's component live in memory, so if you want to persist the work that you have done in your graph you'll need to save a copy of the workspace to the disk, this can be done using.

```python
>>> ws.dump_to_disk("./persisted_workspece.sdgws")
```

and the session can be restored by loading the workspace instead of rebuilding it from scratch every time

```python
## Create a new Workspace
>>> restored_ws = bsg.WorkSpace()
>>> restored_ws.load_from_disk("./persisted_workspece.sdgws")
```

Now restored_ws contains all the same elements as ws. 

Please note that some of the components that are derived from the data ( like graph kmer indexes) are not persisted in the session and need to be recomputed after loading the workspace from disk.

### ws.sdg

sdg as the name implies is graph based, each workspace has a main graph where all the data is referenced to. sdg has to types of graphs sequence digraphs (sdg) and distance graphs (dg). Sequence digraphs are composed of nodes and edges, the nodes store the sequences and the edges represent distances between nodes (distance is the d in sdg). Negative distances indicate overlap and positive distances indicate N gaps. SDG sequences are stored in their canonical form in the nodes. Distance graphs (dg) only store distances (links) and a reference to the nodes in the main sdg, this type of graph is used to store the results of calculations.

![image-20181107213601519](/Users/ggarcia/Library/Application Support/typora-user-images/image-20181107213601519.png)



To import a graph in an empty workspace you can load it directly from a gfa file:

```python
>>> ws.sdg.load_from_gfa('./graph.gfa')
```

Once the graph is loaded all the information can be accessed directly from the `workspace` using dot notation or better using `Nodeviws()` to access the information from the perspective of a node. 

Now that the main graph is in the workspace other data can be ingested by the workspace. SDG supports generic datatypes like paired reads, long reads and linked reads. 

Raw data is converted to DataStores, this is done so SDG can store an indexed version of the data in disk. This allows the framework to have quick access to the entire dataset without the need to have the entire dataset in memory, the result is that workspaces can be bigger than the available ram. Datastores are type specific, each available data type is stored according to it's main characteristic.

Datastores can be created using the python API or the command line utilities. For example, to create a paired end datastore using the API.

```python
>>> SDG.PairedReadsDatastore_build_from_fastq("prds.prds", "./pe-reads_R1.fastq", "./pe-reads_R2.fastq", "pe_reads")
```

the `.PairedReadsDatastore_build_from_fastq()` method created a paired read datastore form fastq files and stores the datastore in the disk ready to be attached to a workspace.

```python
>>> ws.add_paired_reads_datastore("./prds.prds", "PE")
```

`.add_paired_reads_datastore()` attaches the datastore to the workspace and names it `pe` in this case. now the dataset is attached to the workspace. The attached datastores can be listed using the .`list_<datatype>_datastores()` method, in this case:

```python
>> ws.list_paired_reads_datastores()
('PE',)
```

an alternative way of creating datastores is using the command line interfase, in this case `sdg-datastore` 

```bash
$ sdg-datastore make -t paired -o cli-prds.prds -1 ./pe-reads_R1.fastq -2 ./pe-reads_R2.fastq
```

this will create a datastore that can be attached in the same way as  to the ws

```python
>>> ws.add_paired_reads_datastore("./cli-prds.prds", "cli-PE")
>>> ws.list_paired_reads_datastores()
('PE', 'cli-PE',)
```

in a similar way linked reads, long reads datastores and kmer counters can be attached to the workspace 

```python
## Create linked reads ds and attach it to the ws
>>> SDG.LinkedReadsDatastore_build_from_fastq("./lirds.lirds", "li_reads", "./child/child-link-reads_R1.fastq.gz", "./child/child-link-reads_R2.fastq.gz", SDG.LinkedReadsFormat_raw, readsize=250, chunksize=10000000)
>>> ws.add_linked_reads_datastore("./lirds.lirds", "LI")
>>> ws.list_linked_reads_datastores()
('LI',)

## Create long reads ds and attach it to the ws
>>> SDG.LongReadsDatastore.build_from_fastq("./lords.lords", "longreads_pe", "./child/child-long-reads.fastq")
>>> ws.add_long_reads_datastore('./lords.lords', 'LO')
>>> ws.list_long_reads_datastores()
('LO',)

## Create kmer counter collection in the ws
>>> ws.add_kmer_counter("kmers-PE-k27", 27)
>>> ws.get_kmer_counter("kmers-PE-k27").add_count("PE", ["./pe-reads_R1.fastq", "./pe-reads_R2.fastq"])
>>> ws.list_kmer_counters()
('kmers-PE-k27', )
>>> ws.get_kmer_counter("kmers_noncanonical").list_names()
("PE",)
```

now all data is attacched to the workspace and can be usad from within. 

The next thing to do is to indicate SDG how the data in the datastores related to the graph, this is how we want to map the data from the datastores in the nodes. 

Each datastore in the ws gets a mapper object within, the mapper is used to map the data to the nodes of the graph. 

The mapper is in charge of mapping the data back to the nodes, the result of a mapping the data is a set of `Mappings` or a collection of 

Each data type has a particular way of mapping the data back to the node and produces different mappings. For example to map long reads from the datastore:

```python
## set the mapper k to 15
ws.get_long_reads_datastore("long_pe").mapper.k=15
## execute the read mapping
ws.get_long_reads_datastore("long_pe").mapper.map_reads()
```

alternatively a pointer to the mapper can be constructed to make the code more readable, the mapper is pointed by the variable but is still contained within the datastore in the mapper.

```python
## pointer to the long read mapper
lorm = ws.get_long_reads_datastore("long_pe").mapper

## set k and map the reads
lorm.mapper.k=15
lorm.map_reads()
```

the resulting mappings are stared inside the mapper in the datastore, for example if you want to check wich nodes were mapped to the 1000th read you can user the `get_raw-mappings_from_read()` function but the best way to access the mappings is using `NodeViews` (see next section).

```python
>>> for m in lorm.get_raw_mappings_from_read(1000):
>>> .... print (m)
LongReadMapping 1000 (0:8302) -> -429 (5451:13753)  1169 hits
LongReadMapping 1000 (8530:14690) -> -4990 (190:6340)  963 hits
LongReadMapping 1000 (15831:20460) -> 2904 (8:4637)  611 hits
LongReadMapping 1000 (24032:25491) -> -4169 (14:1473)  213 hits
LongReadMapping 1000 (24032:25491) -> -2393 (14:1473)  259 hits
LongReadMapping 1000 (27748:29666) -> 2983 (21:1939)  240 hits
LongReadMapping 1000 (29531:30795) -> 4596 (12:1276)  215 hits
LongReadMapping 1000 (29531:30795) -> 4773 (12:1276)  127 hits
```

The same mappings can be applied to the rest of the datastores

```python
## Map the rest of the ds in the workspace
>>> ws.get_paired_reads_datastore("PE").mapper.map_reads()
>>> ws.get_linked_reads_datastore("LI").mapper.map_reads()
```

Now that all the data is atteched to the ws and the reads alre mapped it is a good time to save the workspace.

```python
ws.dump_to_disk("./persisted.reads.mapped.sdgws")
```

Now all the work is stored in disk and can be loaded when needed.



### NodeViews

NodeViews are used to access the workspace data from a node perspective, to create a nodeview of the 10th node of the main sdg graph in the workspace:

```python
>>> nv = ws.sdg.get_nodeview(10)
```

using nodeviews you can access node properties, navigate the graph, access datastores and mappings. For example to see central node properties:

```python
## Get the node id
>>> print(nv.node_id())
10

## get the node sequence
>>> print(nv.sequence())
'AAAAAAAAAACTGTTAATTTTTTTTTTTACTGATTTACAAATTCCTTGTACCTCTATCTCAGCCCATGCATATATTGGTTGTTTTACTGCTTGATACCACTAGAAACCTCTAGTTGTTTAA'

## get the node size
>>> print(nv.size())
121
```

also topology of the graph from the central node perspective, for example, if you want to get the node neighbours

```python
>>> for n in nv.netx():
>>> ....print(n)
<< NodeDistanceView: -60bp to 2 >>
<< NodeDistanceView: -60bp to -3 >>

>>> for p in nv.netx():
>>> ....print(p)
<< NodeDistanceView: ...bp to ... >>
<< NodeDistanceView: ...bp to ... >>
```

Mappings can also be queried from the node view, 

```python
## Get the names of the available linked reads datastore
>>> ws.list_linked_reads_datastores()
('LI',)

## get the id of the first 3 reads mapped to the node
>>> for m in nv.get_linked_mappings('LI')[:3]:
>>> ....print(m)
12576
23882
30795

## get the first 3 mappings from the datastore
>>> for m in nv.get_linked_mappings('LI')[:3]:
>>> ....print(m)
<pysdg.pysdg.ReadMapping; proxy of <Swig Object of type 'ReadMapping *' at 0x107bb5f90> >
<pysdg.pysdg.ReadMapping; proxy of <Swig Object of type 'ReadMapping *' at 0x107bb55d0> >
<pysdg.pysdg.ReadMapping; proxy of <Swig Object of type 'ReadMapping *' at 0x107bb5f90> >
```

mappings contain the information to go from nodes to reads and vice versa, the information inside the mappings can be queried individually

```python
### Get the first mappings
>>> mapping = nv.get_linked_mappings('LI')[0]
## Check the node of the mapping
>>> mapping.node
10

## Check the read id
>>> mapping.read_id
12576

## Check the first and last mapping position in the node for the mapping
>>> mapping.first_pos
641
>>> mapping.last_pos
638

## Check the orientation of the mapping (True=opposite to node, False=Same as node)
>>> mapping.rev
Out[27]: True
```

Each data type has a slightly different mapping format.

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

