# ORTHOSCOPE
Our web servise is available from 
http://orthoscope.jp  
[Japanese instruction](http://www.fish-evol.com/orthoscope_ji.html)

---

## Mode
![mode](images/mode.jpg)

---

## Case Studies in Inoue and Satoh (2018)
Query seqeunces from genes with known function.

| Actinopterygii | Vertebrata | Deuterostomia | Protostomia |
:---: | :---: | :---: | :---:
| PLCB1* | ALDH1A* | Brachyury | Brachyury |
| [Queries][t1-1] | [Queries][t1-2] | [Queries][t1-3] | [Queries][t1-4] |
| [Result][t1-5] | [Result][t1-6] | [Result][t1-7] | [Result][t1-8] |

[t1-1]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/ActinopterygianPLCB1.fas
[t1-2]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/VertebrateALDH1A.fas
[t1-3]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/DeuterostomeBra.fas
[t1-4]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/ProtostomeBra.fas
[t1-5]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/ActinopterygianPLCB1.zip
[t1-6]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/VertebrateALDH1A.zip
[t1-7]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/DeuterostomeBra.zip
[t1-8]:https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/ProtostomeBra.zip

#         
Query sequence collectoin from assemble database*

1. Download Coregonus lavaretus TSA file ([GFIG00000000.1](https://www.ncbi.nlm.nih.gov/nuccore/GFIG00000000.1)) form NCBI.
2. Translate raw sequences into amino acid and coding sequences using [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki).
```
./TransDecoder.LongOrfs -t GFIG01.1.fsa_nt
```
3. Make blast databases using [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
```
makeblastdb -in longest_orfs.pep -dbtype prot -parse_seqids 
makeblastdb -in longest_orfs.cds -dbtype nucl -parse_seqids
```
4. BLASTP seaech against amino acid database.
```
blastp -query query.txt -db longest_orfs.pep -num_alignments 10 -evalue 1e-12 -out 010_out.txt
```
5. Retrieve blast top hit sequences from coding sequence file using sequence id.
```
blastdbcmd -db longest_orfs.cds -dbtype nucl -entry_batch queryIDs.txt -out 020_out.txt
```

<br />
<br />

---
## Focal Group
![analisis group](images/analysisGroup.jpg)

<br />
<br />  

---
## Upload Files
Coding sequence

![file format](images/UplodFile.jpg)

Case 1: Query seqeunce is present in the ORTHOSCOPE database

![registered sequence search](images/example1.jpg)

Case 2: Query seqeunce is not present in the ORTHOSCOPE database

![unregistered sequence search](images/yourOwnSequence.jpg)

<br />
<br />  

---
## Species Tree Hypothesis

![SpeciesTree](images/SpeciesTree.jpg)

Our hypothetical species tree (newick) can be downloaded from [here](https://fish-evol.unit.oist.jp/orthoscope/examples/SpeciesTreeHypothesis.tre).


| [Metazoa][treeA] | [Hexapoda][treeB] | [Vertebrata][treeC] | [Aves][treeD] | [Actinopterygii][treeE] |
:---: | :---: | :---: | :---: | :---:

[treeA]:https://github.com/jun-inoue/orthoscope/raw/master/images/SpeciesTree_Metazoa.pdf
[treeB]:https://github.com/jun-inoue/orthoscope/raw/master/images/SpeciesTree_Hexapoda.pdf
[treeD]:https://github.com/jun-inoue/orthoscope/raw/master/images/SpeciesTree_Aves.pdf
[treeC]:https://github.com/jun-inoue/orthoscope/raw/master/images/SpeciesTree_Vertebrata.pdf
[treeE]:https://github.com/jun-inoue/orthoscope/raw/master/images/SpeciesTree_Actinopterygii.pdf


Phylogenetic relationships without references follow the [NCBI Taxonomy Common Tree](https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi).

Newick formats can be modifed using [TreeGraph2](http://treegraph.bioinfweb.info/).

![treegraph2](images/treeGraph2.jpg)

<br />
<br />  

---
## Sequence Collection
![sequence collection](images/BlastEvalue.jpg)

<br />
<br />  

---
## Alignment
![sequence alignment](images/Aligned-site_rate.jpg)

<br />
<br />  

---
## Tree Search
Dataset

![codon mode](images/dataset.jpg)

<br />
Rearrangement BS value threshold 

![branch rearrangement](images/rearrangeBS.jpg)

NJ analysis is conducted using the software package [Ape](https://cran.r-project.org/web/packages/ape/ape.pdf) in R (coding) and [FastME](http://www.atgc-montpellier.fr/fastme/) (amino acid). Rearrangement analysis is done using a method implemented in [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/).

<br />
<br />  

---
## Genome Taxon Sampling

Feasibility of completion

Number of hits to report per genome | Number of species
:---: | :---:
3 | <50
5 | <40 
10 | <30 

<br />
<br />  

---
## Tree Estimation of Orthogroup Members using Additional Sequences
The script is specialized for a Macintosh use. Windows users need some modifications.
[Example](https://github.com/jun-inoue/orthoscope/raw/master/tarfiles/DeuterostomeBra_2ndAnalysis.zip).

### Installing Dependencies

Estimation of the small tree requires some dependencies to be installed and in the system path.
<br />  


#### RAxML:

Available here: [https://github.com/stamatak/standard-RAxML](https://github.com/stamatak/standard-RAxML)

Download the the latest release and extract it.
Cd into the extracted directry (e.g., standard-RAxML-8.2.12), compile the PThreads version, and copy the executable to a directory in your system path, e.g.:
```
cd standard-RAxML-8.2.12
make -f Makefile.SSE3.PTHREADS.gcc
cp raxmlHPC-PTHREADS-SSE3 ~/bin
```
Add the directory containing the directory to your PATH variable. e.g.:
```
export PATH=$PATH:~/bin
``` 
<br />  

#### Mafft:
Available here: [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)

After compilation, set your PATH [following this site](https://mafft.cbrc.jp/alignment/software/add_path.html).  

<br />   


#### trimAl v1.2 (Official release):
Available here: [http://trimal.cgenomics.org/downloads](http://trimal.cgenomics.org/downloads).  
Cd to trimAl/source, type make, and copy the executable.
```
make
cp trimal ~/bin
```  
<br />  


#### PAL2NAL: 
Available here: [http://www.bork.embl.de/pal2nal/#Download](http://www.bork.embl.de/pal2nal/#Download)
Change the permission of perl script and copy it.
```
chmod 755 pal2nal.pl
cp pal2nal.pl ~/bin
```  
<br />  


#### Ape in R:
R is availab here [R](https://cran.ism.ac.jp). By installing R, [rscript](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html) will be installed automatically. 
[APE in R](http://ape-package.ird.fr) can be installed from the R console.
```
install.packages("ape")
```
<br />
<br />  

### Tree Estimation
1. Select an appropriate outgroup and orthogroup members and save 010_candidates_nucl.txt file. The outgroup sequence should be placed at the top of alignment. Additional sequences can be included.

[![query sequences](images/treeSearchWithOrthologs.jpg)](images/treeSearchWithOrthologsL.jpg)

2. Decompress 100_2ndTree.tar.gz file.
3. cd into 100_2ndTree file.
4. Run the pipeline.
```
./100_estimate2ndTree.py
```
5. ML tree is saved in 200_RAxMLtree_Exc3rd.pdf automatically.

[![ML tree](images/200_RAxMLtree_Exc3rd.jpg)](images/200_RAxMLtree_Exc3rdL.jpg)

<br />  


### Duplicated Node Estimation
Using [Notung](http://www.cs.cmu.edu/~durand/Notung/), duplicated nodes can be identified. Here, we will analyze the gene tree of orthogroup members.

1. Double click the downloaded .jar file (here, Notung-2.9.jar).  
2. Save the species tree (newick format) as a new file (here, speciesTree.tre), from 000_summary.txt file.  
3. Open the species tree file, speciesTree.tre (File > Open Gene Tree), from Notung.  
4. Open the gene tree file, RAxML_bootstrap.txt (File > Open Gene Tree).  
5. Set "Edge Weight THreshold" (here 70) from “Edit Values button“. This value corresponds to
“Rearrangement BS value threshold” in ORTHOSCOPE.  
6. From "Rearrange" tab in the bottum, select "Prefix of the general label".  
7. Push "Reconcile" button.  
8. Duplicated nodes are shown with "D".  

[![Rearranged tree](images/Notung_rearrangement1.jpg)](images/Notung_rearrangement1L.jpg)

<br />
<br />  
 


---
## Supported Browsers
Chrome | Firefox | Safari | IE
:---: | :---: | :---: | :---:
Supported | Supported | 11.0 or later | Not supported
<br />
<br />  

---
## History

Date | Version | Revision
--- | --- | ---
21 Dec. 2018 | [Version 1.0.1.](http://fish-evol.unit.oist.jp/orthoscope11/) | Released. In the rearranged gene tree, nodes identified as speciation events were marked with "D".
18 Dec. 2018 | Version 1.0.1.beta | Some xenacoelomorph, platyhelminth, priapulid, avian data were newly added.
10 July 2018 | [Version 1.0](http://fish-evol.unit.oist.jp/orthoscope100) | Published in Inoue and Satoh (2018).

<br />
<br />  

---
## Database
Available from [here](https://zenodo.org/record/1452077#.W7xEfS_ANsM)
(10.5281/zenodo.1452077). 10 October 2018.

ORTHOSCOPE employs a genome-scale, protein-coding gene database (coding and amino acid sequence datasets) constructed for each species. In order to count numbers of orthologs in each species, only the longest sequence is used, when transcript variants exist for single locus.
<br />
<br />  

---
## Citation
Inoue J. and Satoh N. ORTHOSCOPE: an automatic web tool of analytical pipeline for ortholog identification using a species tree. in press. [MBE](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msy226/5229933).

---
Previous versions: 


Email: [_jun.inoue_ AT _oist.jp_](http://www.geocities.jp/ancientfishtree/index_eng.html)
