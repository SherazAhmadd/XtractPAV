# XtractPAV  

A pipeline for extracting Presence-Absence Variations (PAVs) in eukaryotic and prokaryotic genomes across multiple assemblies.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
  - [Using Conda](#installation)
  - [Using Debian](#installation)
- [Parameters](#parameters)
- [Input & Output Files](#input--output-files)
- [Run XtractPAV with Test Data](#run-xtractpav-with-test-data)
- [Reference](#reference)
- [Contact Us](#contact-us)


## Introduction

**PAV (Present/Absence Variation)** is a form of structural genetic variation where specific genes or genomic regions are completely present in some individuals or species but entirely absent in others. Unlike small-scale variations like SNPs, PAVs involve larger segments of DNA and can significantly impact phenotype by altering gene content. Common in both eukaryotic and prokaryotic genomes, PAVs contribute to genetic diversity, adaptation, evolution, and traits such as disease resistance or pathogenicity. They often result from gene duplication, deletion, or horizontal gene transfer, and their analysis is crucial in fields like evolutionary biology, agriculture, microbial genomics, and personalized medicine.


**XtractPAV** is an automated pipeline, designed to extract **Presence/Absence Variations (PAVs)** from genomic datasets. The pipeline utilizes **Mummer4** for the comparative analysis of genomes and incorporates custom Python scripts for the extraction of raw PAVs. the pipeline applies a set of coverage and similarity criteria  to identify authentic PAVs. **XtractPAV** is capable of identifying **genic PAVs** and provides annotations for these variations. The pipeline generates a final report, which includes interactive visualizations illustrating the distribution and characteristics of the PAVs across the genomes analyzed.

## Features
- Extracts the authentic PAVs and allows users to set a minimum length threshold for the PAVs.
- Support for multiple genomes, enabling the tracking of variations across diverse samples.
- Compatibility with various genome file formats, accepting extensions such as .fa, .fasta, .fna, or any other common genome file suffix.
- Usability across a wide range of genomes, including both eukaryotic and prokaryotic organisms.
- Detection of genic PAVs, accompanied by functional annotations for the identified variations.
- Generation of a comprehensive final report, detailing the PAVs, and providing interactive visualizations that illustrate their distribution and characteristics.
## Installation

Using *XtractPAV* is very easy. Simply clone the repository or use the wget and uncompressed it.

```
git clone https://github.com/SherazAhmadd/XtractPAV
cd XtractPAV/XtractPAV-pipeline
# Add the bin to PATH
export PATH=$PATH:/path/to/XtractPAV
```


- #### using conda
    To set up the environment using Conda, run the following:

    ```
    conda env create -f XtractPAV-Dependencies.yml
    conda activate XtractPAV  
    ```

- #### Using Ubuntu/Debian-based systems
    To install the required tools and dependencies on Ubuntu/Debian-based systems, follow these steps:

    1. **MUMmer4**:
       MUMmer is available [**here**](http://mummer.sourceforge.net/), with **MUMmer-4.0.0** being the version used in this study. It is important to have a **GCC compiler** (**g++ version ≥ 4.7**) for proper installation.

       ```
       sudo apt install mummer4
       
       or

        wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
        tar -xvzf mummer-4.0.0beta2.tar.gz
        ./configure --prefix=/path/to/installation
        make
        make install
       # Add MUMmer tools to your PATH
        export PATH=/path/to/installation/:$PATH
       ```

    2. **BLAST+**:
       You can install BLAST+ by following [**here**](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). The latest version is used in this study.

       For Ubuntu/Debian-based systems, run:

       ```
       sudo apt install ncbi-blast+

       or

        wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
        tar zxvf ncbi-blast-2.7.1+-x64-linux.tar.gz
       # Add Blast+ tools to your PATH
        export PATH=/path/to/blast+/bin:$PATH
       ```
    3. **Bedtools**    
      Bedtools - the swiss army knife for genome arithmetic. You can download via [**here**](https://bedtools.readthedocs.io/en/latest/). XtractPAV is not version specific to bedtool however version 2.27.1 or the latest version is preferred.
         ```
         sudo apt install bedtools

       or

        tar -zxvf bedtools-2.27.1.tar.gz
        cd bedtools2
        make
       # Add Bedtools tools to your PATH
         export PATH=/path/to/bedtools/bin:$PATH

         ```
    4. **Python dependencies**:
       Install the necessary Python packages:

       ```
       pip install biopython==1.78
       pip install plotly==6.0.1

       ```
XtractPAV currently only supports Linux system due to the software dependencies
## Parameters

The following parameters can be used to customize the execution of the **XtractPAV** pipeline:

| Parameter       | Description                                                                                     | Example Usage                          |
|------------------|-------------------------------------------------------------------------------------------------|----------------------------------------|
| `--rf`          | Path to the reference genome file.                                                              | `--rf reference_genome.fa`             |
| `--ra`          | Path to the reference GFF3 annotation file.                                                     | `--ra reference_annotation.gff3`       |
| `--qf`          | Comma-separated list of query genome files.                                                     | `--qf "query1.fa,query2.fa"`             |
| `--qa`          | Comma-separated list of query GFF3 annotation files.                                            | `--qa "query1.gff3,query2.gff3"`         |
| `--cov`         | Minimum coverage threshold for filtering PAVs.                                                  | `--cov 0.8`                             |
| `--sim`         | Minimum similarity percentage for filtering PAVs.                                               | `--sim 90.0`                             |
| `--len`         | Minimum length threshold for PAVs.                                                              | `--len 100`                            |
| `--thr`         | Number of threads to use for parallel processing.                                               | `--thr 8`                              |
| `--help`        | Display help information and usage instructions.                                                | `--help`                               |
| `--version`     | Display the version of the XtractPAV pipeline.                                                    | `--version`                            |
|

*Note:* Coverage and similarity thresholds must be in Float format (e.g., `0.8` for 80% coverage or `90.0` for 90% similarity). The length threshold must be an integer (e.g., `100` for a minimum length of 100 base pairs). 

### Example Command
To run **XtractPAV** with custom parameters:
```
XtractPAV.sh --rf reference_genome.fa --ra reference_annotation.gff3 --qf query1.fa,query2.fa --qa query1.gff3,query2.gff3 --cov 0.8 --sim 90.0 --len 100 --thr 8
```
# Input & Output Files
## Input Files

To run **XtractPAV**, you need at least one reference genome file and one or more query genome files, along with their corresponding annotation files.

#### Genome Sequence Files
The genome sequence files should be in FASTA format with the following structure:

```
>chromosome1
ATCGATCG...
```

File extensions such as `.fa`, `.fasta`, `.fna`, or other common genome file suffixes are supported. The prefix of the sequence file name will be used to generate temporary files, so it is recommended to use descriptive names (e.g., `species1.fa` or `sample1.fasta`).

#### Annotation Files
Annotation files should be in GFF3 format and include information about genes, mRNA, exons, and CDS. Below is an example of a valid GFF3 file:

```
ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=GENE1
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=GENE1.1
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00001
ctg123 . CDS             1201  1500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=GENE1_protein
```

Annotations in GFF format containing `gene` information lines are also supported by **XtractPAV**.


## Output Files

The **XtractPAV** pipeline generates the following output files:

####  PAV Sequences (FASTA Format)
This file contains the extracted PAV sequences in FASTA format. Each sequence header includes the chromosome name and coordinates of the PAV. Below is an example:

```
>NC_010067.1__1__7443
ATCGATCG...
>NC_010067.1__70167__86803
GCTAGCTA...
>NC_010067.1__150997__157358
TGCATGCA...
```

####  Genic PAVs (GFF3 Format)
This file provides information about genic PAVs and their functional annotations. Only `gene` entries are included. Below is an example:

```
NC_010067.1 . gene 1       7443    . + . ID=gene00001;Name=GENE1
NC_010067.1 . gene 70167   86803   . + . ID=gene00002;Name=GENE2
NC_010067.1 . gene 150997  157358  . + . ID=gene00003;Name=GENE3
```

####  Final Report (HTML Format)
The final report includes interactive visualizations and a detailed table of PAVs. Below is an outline of the report:

- **Figures**:
  1. **PAV Length Distribution**  
     ![PAV Length Distribution](https://github.com/SherazAhmadd/XtractPAV/blob/main/PAV_length.png)
  2. **Number of PAVs per Chromosome**  
     ![Number of PAVs per Chromosome](https://github.com/SherazAhmadd/XtractPAV/blob/main/PAVs_per_chromosome.png)
  3. **Percentage of PAVs per Chromosome**  
     ![Percentage of PAVs per Chromosome](https://github.com/SherazAhmadd/XtractPAV/blob/main/PAVS_Distribution.png)

- **PAV Details Table**:
  | Number of PAV | Chromosome       | Start   | End     | PAV Length |
  |---------------|------------------|---------|---------|------------|
  | PAV1          | NC_010067.1      | 1       | 7443    | 7443       |
  | PAV2          | NC_010067.1      | 70167   | 86803   | 16637      |
  | PAV3          | NC_010067.1      | 150997  | 157358  | 6362       |


# Run XtractPAV with Test Data
#### Test Data
The test data is available in the `sample_genomes` directory. You can use the provided test data to run **XtractPAV** and verify its functionality.
#### Run XtractPAV
To run **XtractPAV** with the test data, use the following command:

```
XtractPAV.sh --rf S_Entrica_LT2.fna --ra S_Entrica_LT2.gff --qf S_Agona_SL483.fna --qa S_Agona_SL483.gff --cov 0.9 --sim 95.0 --len 100 --thr 8
```

# Reference
Rana Sheraz Ahmad, Muhammad Tahir ul Qamar, Heng Li. XtractPAV: An Automated Pipeline for Identifying Presence–Absence Variations Across Multiple Genomes Genomes. *Bioinformatics*, [https://doi.org/10.1093/bioinformatics/example](https://doi.org/10.1093/bioinformatics/example)


# Contact Us
For any questions or issues, please contact us at:
- **Rana Sheraz Ahmad**: [ranasheraz.202101902@gcuf.edu.pk](mailto:ranasheraz.202101902@gcuf.edu.pk)
- **Muhammad Tahir ul Qamar**: [m.tahirulqamar@hotmail.com](mailto:m.tahirulqamar@hotmail.com)

