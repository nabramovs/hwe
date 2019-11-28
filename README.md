# Hardy-Weinberg Equilibrium in the Large Scale Genomic Sequencing Era.
Analysis of deviations from Hardy-Weinberg equilibrium in a large population database (gnomAD).

### Project Description

Hardy-Weinberg Equilibrium (HWE) is used to estimate the number of homozygous and heterozygous variant carriers based on its allele frequency in populations that are not evolving. Previously, deviation from HWE in large population databases were investigated to detect genotyping errors, which can result in extreme heterozygote excess (HetExc). However, HetExc might also be a sign of natural selection since recessive disease causing variants are expected to occur less frequently in a homozygous state in the general population, but might reach high allele frequency, especially if they are advantageous, in a heterozygote state. We developed a filtering strategy to detect these variants and applied it on genome data from 137,842 individuals. We found that the main limitations of this approach were quality of genotype calls and insufficient population sizes, whereas population structure and high level of inbreeding could reduce  sensitivity, but not precision, in certain populations. Nevertheless, we identified 365 HetExc variants in 326 genes, most of which were specific to African/African American populations (~84.7%). Although the majority of them were not associated with known diseases, or were classified as “benign”, they were enriched in genes associated with autosomal recessive diseases. The resulting dataset also contained two known recessive disease causing variants with evidence of heterozygote advantage in the genes HBB and CFTR. Finally, we provide in silico evidence of a novel heterozygote advantageous variant in the CHD6 gene (involved in influenza virus replication). We anticipate that our approach will allow the detection of rare recessive disease causing variants in the future. 

### Usage

1) Install mongoDB, see installation instructions here:
https://docs.mongodb.com/manual/installation/

2) Install required python modules:
```
pip install -r requirements.txt
```
3) Unzip the source datasets:
```
unzip source_data.zip
```
4) Import datasets into the local database: 
```
python import_data.py
```
5) Perform the analysis:
```
python hwe_analysis.py
```
6) Produce all the figures and tables:
```
python figures.py
```

**IMPORTANT NOTE**: Before running the above scripts, check their main methods and uncomment all the functions required to perform the analysis. Functions which work with gnomAD API require significant amount of time to run (up to a couple of days), but had to be executed only once, as they store results in the local database. See comments in the main methods for additional information.

### Code Description

- **import_data.py** - imports data from various datasets into the local database

- **hwe_analysis.py** - performs the main analysis

- **figures.py** - produces all the figures, tables and reports statistics

- **hardy_weinberg.py** - contains the method used to calculate Hardy-Weinberg Equilibrium statistics 

- **gnomad_api.py** - contains methods used to obtain data from gnomAD API

- **common.py** - contains methods commonly used by other scripts

- **csv_reader.py** - custom csv reader, used to import data into local database
