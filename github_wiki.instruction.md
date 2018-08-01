# ORTHOSCOPE: instruction
This page is an instruction for [ORTHOSCOPE](https://fish-evol.unit.oist.jp/yurai/index.html).

## Two analyses (mode) of ORTHOSCOPE
![mode](https://fish-evol.unit.oist.jp/yurai/examples/mode.jpg)

## Query seqences: Inoue et al.
Query seqeunces from genes with known function.

| Actinopterygii | Vertebrata | Deuterostomia | Protostomia |
:---: | :---: | :---: | :---:
| PLCB1* | ALDH1A | Bra | Bra |
| [Queries][t1-1] | [Queries][t1-2] | [Queries][t1-3] | [Queries][t1-4] |
| [Result][t1-5] | [Result][t1-6] | [Result][t1-7] | [Result][t1-8] |

[t1-1]:http://fish-evol.unit.oist.jp/yurai/examples/ActinopterygianPLCB1.txt.tar.gz
[t1-2]:http://fish-evol.unit.oist.jp/yurai/examples/VertebrateALDH1A.txt.tar.gz
[t1-3]:http://fish-evol.unit.oist.jp/yurai/examples/DeuterostomeBra.txt.tar.gz
[t1-4]:http://fish-evol.unit.oist.jp/yurai/examples/ProtostomeBra.txt.tar.gz
[t1-5]:http://fish-evol.unit.oist.jp/yurai/examples/ActinopterygianPLCB1.tar.gz
[t1-6]:http://fish-evol.unit.oist.jp/yurai/examples/VertebrateALDH1A.tar.gz
[t1-7]:http://fish-evol.unit.oist.jp/yurai/examples/DeuterostomeBra.tar.gz
[t1-8]:http://fish-evol.unit.oist.jp/yurai/examples/ProtostomeBra.tar.gz

Query sequences from assemble database*

1. Download Coregonus lavaretus TSA file ([GFIG00000000.1](https://www.ncbi.nlm.nih.gov/nuccore/GFIG00000000.1)) form NCBI.
2. Translate raw sequences into amino acid and cDNA sequences using [TrandDecoder](https://github.com/TransDecoder/TransDecoder/wiki).
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
5. Retrieve blast top hit sequences from cDNA file using seq id.
```
blastdbcmd -db longest_orfs.cds -dbtype nucl -entry_batch queryIDs.txt -out 020_out.txt
```

## Query seqences: example

|   | Actinopterygii | Vertebrata | Deuterostomia | Protostomia |
| :---: | :---: | :---: | :---: | :---:
| Bra | [3queries.txt][t2-1] | [6queries.txt][t2-2] | [6queries.txt][t2-3] | [3queries.txt][t2-4] |
| Actin | [4queries.txt][t2-5] | [4querist.txt][t2-6] | [7queries.txt][t2-7] | [9queries.txt][t2-8] |
| MHC | [4queries.txt][t2-9] | [10queries.txt][t2-10] | [6queries.txt][t2-11] | [2queries.txt][t2-12] |

[t2-1]:https://fish-evol.unit.oist.jp/yurai/examples/Actinopterygii_Bra3queries.txt.tar.gz
[t2-2]:https://fish-evol.unit.oist.jp/yurai/examples/Vertebrata_Bra6queries.txt.tar.gz
[t2-3]:https://fish-evol.unit.oist.jp/yurai/examples/Deuterostomia_Bra6queries.txt.tar.gz
[t2-4]:https://fish-evol.unit.oist.jp/yurai/examples/Protostomia_Bra3queries.txt.tar.gz
[t2-5]:https://fish-evol.unit.oist.jp/yurai/examples/Actinopterygii_Actin_4queries.txt.tar.gz
[t2-6]:https://fish-evol.unit.oist.jp/yurai/examples/Vertebrate_Actin_4queries.txt.tar.gz
[t2-7]:https://fish-evol.unit.oist.jp/yurai/examples/Deuterostomia_Actin_7queries.txt.tar.gz
[t2-8]:https://fish-evol.unit.oist.jp/yurai/examples/Protostomia_Actin_9queries.txt.tar.gz
[t2-9]:https://fish-evol.unit.oist.jp/yurai/examples/Actinopterygii_MHC_4queries.txt.tar.gz
[t2-10]:https://fish-evol.unit.oist.jp/yurai/examples/Vertebrate_MHC_10queries.txt.tar.gz
[t2-11]:https://fish-evol.unit.oist.jp/yurai/examples/Deuterostomia_MHC6spp.txt.tar.gz
[t2-12]:https://fish-evol.unit.oist.jp/yurai/examples/Protostomia_MHC.txt.tar.gz

## Analysis group
![analisis group](https://fish-evol.unit.oist.jp/yurai/examples/analysisGroup.jpg)

## Upload files
Contents

![file format](https://fish-evol.unit.oist.jp/yurai/examples/UplodFile.jpg)

Case 1: Query seqeunce is present in ORTHOSCOPE database

![registered sequence search](https://fish-evol.unit.oist.jp/yurai/examples/example1.jpg)

Case 2: Query seqeunce is not present in ORTHOSCOPE database

![unregistered sequence search](https://fish-evol.unit.oist.jp/yurai/examples/yourOwnSequence.jpg)

---
Species tree hypothesis

![species tree](https://fish-evol.unit.oist.jp/yurai/examples/Metazoa.tre.jpg)

The tree file can be modifed using [TreeGraph2](http://treegraph.bioinfweb.info/).

![treegraph2](https://fish-evol.unit.oist.jp/yurai/examples/treeGraph2.jpg)

## Sequence collection
![sequence collection](https://fish-evol.unit.oist.jp/yurai/examples/BlastEvalue.jpg)

## Alignment
![sequence alignment](https://fish-evol.unit.oist.jp/yurai/examples/Aligned-site_rate.jpg)

## Tree search
Dataset

![codon mode](https://fish-evol.unit.oist.jp/yurai/examples/dataset.jpg)

Rearrangement BS value threshold 

![branch rearrangement](https://fish-evol.unit.oist.jp/yurai/examples/rearrangeBS.jpg)

NJ analysis is conducted using the software package [Ape](https://cran.r-project.org/web/packages/ape/ape.pdf) in R (cDNA) and [FastME](http://www.atgc-montpellier.fr/fastme/) (amino acid). Rearrangement analysis is done using a method implemented in [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/).

## Tree estimation using identified orthologs
Mac only
1. Select only orthologs and save 010_candidates_nucl.txt file

![query sequences](https://fish-evol.unit.oist.jp/yurai/examples/treeSearchWithOrthologsL.jpg)

2. Decompress 100_2ndTree.tar.gz file
3. cd into 100_2ndTree file
4. Run the pipeline
```
./100_estimate2ndTree.py
```
5. ML tree is saved is 200_RAxMLtree_Exc3rd.pdf.

![ML tree](https://fish-evol.unit.oist.jp/yurai/examples/200_RAxMLtree_Exc3rd.jpg)

## Supported browsers
Chrome | Firefox | Safari | IE
:---: | :---: | :---: | :---:
OK | OK | 11.0 or later | Not supported

## History
10 July 2018 	Version 1.0.

## Citation
Inoue J. and Satoh N. ORTHOSCOPE: an automatic web tool of analytical pipeline for ortholog identification using a species tree. in prep.

---
Previous versions: [_jun.inoue_ AT _oist.jp_](http://www.geocities.jp/ancientfishtree/index_eng.html)