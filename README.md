Fungal species associated with apple tree roots
=======

This repository contains data and scripts used to analyze MiSeq data for PLP847 Advanced Mycology course, fall 2016 

####This workflow can be fully reproduced by running the following scripts provided all nessesary software is installed (see below) 
* First run the ITS1_fung_pipe.sh for sequence preprocessing
* Then run the Analysis.Rmd for R analysis

However as downloaded you can run just the R analysis as well 

###Objectives of study 
* Examine how fungal species diversity changes based on orchard managment type (conventional management vs. abandon orchard)
* See how fungal species diversity changes based on how far the sample was taken from the tree trunk 

###Hypotheses 
* Fungal species diversity will be lower in the conventionally managed orchard 
* Fungal species diversity will be different based on how far the sample was taken from the tree 

###Methods
Fungal ITS ampmlicons were amplified with primers ITS1f and ITS4 and subsiquently sequenced on an Illumina Miseq. Since these primers amplify the ITS1, 5.8S, and ITS2 reads 1 and 2 cannot be merged. Therefore, for time, we will only show the workflow for the forward reads. The reverse reads are provided and the workflow can be run on the reverse reads if desired. 

Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            		# raw and primary data, are not changed 
    |  |- raw/         		# raw data, will not be altered
    |  |- clean/       		# cleaned data, will not be altered once created
    |  |- input_files	#Input files for phyloseq
    |
    |- code/           # any programmatic code
    |  |- ITS_fung_pipe.sh # pipeline used for cleaning reads, clustering reads, and assigning taxonomy
    |  |- python_scripts/ # custom python scripts to perform various tasks
    |  |- Analysis.Rmd # R-phyloseq-analyses
    |
    |- reference_databases/ # databases used for analyses
    |
    |- Fungal ITS Apple Roots.Rproj # RStudio project for this study

Software
----------

###This workflow requires the following software

* USEARCHv8.1.1861 executable as usearch8 [install USEARCH8.1.1861](http://www.drive5.com/usearch/download.html)
* VSEARCHv2.3.2 executable as vsearch [install VSEARCHv2.3.2](https://github.com/torognes/vsearch)
* MacQIIMEv1.9.1 all scripts executable [install MacQIIMEv1.9.1](http://www.wernerlab.org/software/macqiime)
* cutadaptv1.11 executable as cutadapt [install cutadaptv1.11](http://cutadapt.readthedocs.io/en/stable/installation.html)
* FastQCv0.11.5 executable as fastqc [install FastQCv0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)

Datasets
----------

###Part I - Sequence preprocessing 
* [instructions on running the script](code/ITS1_fung_pipe.md)

###Part II - R analysis 
* [instructions on running the script](code/Analysis.md)

Miscelaneous
----------------
The initial file and directory structure of this project was developed by a group of participants in the Reproducible Science Curriculum Workshop, held at [NESCent] in December 2014. 

To access this repository template use [rr-init repository](https://github.com/Reproducible-Science-Curriculum/rr-init)

[NESCent]: http://nescent.org
[Rmarkdown]: http://rmarkdown.rstudio.com/
