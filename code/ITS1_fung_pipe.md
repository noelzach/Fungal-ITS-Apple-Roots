The following script provides code to do read preprosessing of fungal ITS sequences. This script was mostly written by Gian MN Benucci, Ph.D. (NICO) and was adapted and edited by Zachary Noel.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### This workflow requires the following software

-   USEARCHv8.1.1861 executable as usearch8 [install USEARCH8.1.1861](http://www.drive5.com/usearch/download.html)
-   VSEARCHv2.3.2 executable as vsearch [install VSEARCHv2.3.2](https://github.com/torognes/vsearch)
-   MacQIIMEv1.9.1 all scripts executable [install MacQIIMEv1.9.1](http://www.wernerlab.org/software/macqiime)
-   cutadaptv1.11 executable as cutadapt [install cutadaptv1.11](http://cutadapt.readthedocs.io/en/stable/installation.html)
-   FastQCv0.11.5 executable as fastqc [install FastQCv0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)

This is a small dataset and can be run on a local computer

The following is a set of instructions of what the ITS1\_fun\_pipe.sh is doing To run ITS1\_fung\_pipe.sh simply navigate to the code directory and execute the file.

-   CURRENT DIRECTORY: code
-   CHANGE DIRECTORY: data

``` bash
cd ../data
```

### Start with a sanity check:

#### DESCRIPTION: Make sure our mapping file is ready to go.

-   INPUT: mapping.txt
-   OUTPUT: mapping\_output

-   CURRENT DIRECTORY: data

``` bash
validate_mapping_file.py -m mapping.txt -o mapping_output
```

-   CURRENT DIRECTORY: data
-   CHANGE DIRECTORY: clean

``` bash
mkdir clean
cd clean
```

#### DESCRIPTION: Striping primers and adapters from the raw reads

-   INPUT: R1\_seqs.fastq
-   OUTPUT: R1\_stripped.fastq

-   CURRENT DIRECTORY: clean

``` bash
cutadapt -g CTTGGTCATTTAGAGGAAGTAA -e 0.01 --discard-untrimmed --match-read-wildcar ../raw/R1_seqs.fastq > R1_stripped.fastq
```

#### DESCRIPTION: Get some stats before trimming reads

-   INPUT: R1\_stripped.fastq
-   OUTPUT: R1\_stats\_result\_vsearch.txt

We have to make a decision on where to truncate the reads based on this output file for more information go to [this webpage to learn more](http://drive5.com/usearch/manual/cmd_fastq_stats.html)

-   CURRENT DIRECTORY: clean

``` bash
vsearch -fastq_stats R1_stripped.fastq -log R1_stats_result_vsearch.txt
```

#### DESCRIPTION: We will truncate to a read length of 180 bp and a expected error threshold of 0.1

-   INPUT: R1\_stripped.fastq
-   OUTPUT: R1\_filtered.fasta, R1\_filtered.fastq

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -fastq_filter R1_stripped.fastq -fastq_maxee 0.1 -fastq_trunclen 180 -fastq_maxns 0 -fastaout R1_filtered.fasta -fastqout R1_filtered.fastq
```

#### DESCRIPTION: Trimming 44 bp from the left of the reads to get rid of the conserved 5.8s region

-   INPUT: R1\_filtered.fastq
-   OUTPUT: R1\_trimmed.fasta, R1\_trimmed.fastq

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -fastq_filter R1_filtered.fastq -fastq_stripleft 44 -fastaout R1_trimmed.fasta -fastqout R1_trimmed.fastq
```

-   CURRENT DIRECTORY: clean
-   CHANGE DIRECTORY: fastqc

``` bash
mkdir fastqc
cd fastqc
```

Lets do another sanity check to make sure our sequences are good quality before moving on

#### DESCRIPTION: Do a quality check to make sure the reads going forward are good quality

-   INPUT: R1\_trimmed.fastq
-   OUTUT: R1\_trimmed\_fastqc.zip, R1\_trimmed\_fastqc.html

-   CURRENT DIRECTORY: fastqc

``` bash
fastqc R1_trimmed.fastq
```

-   CURRENT DIRECTORY: fastqc
-   CHANGE DIRECTORY: clean

``` bash
cd ..
```

#### DESCRIPTION: Dereplication, discard any duplicated reads

-   INPUT: R1\_trimmed.fasta
-   OUTPUT: R1\_dereplicated.fasta

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -derep_fulllength R1_trimmed.fasta -fastaout R1_dereplicated.fasta -sizeout -strand both
```

#### DESCRIPTION: Removal of singletons

-   INPUT: R1\_dereplicated.fasta
-   OUTPUT: R1\_no\_singletons.fasta

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -sortbysize R1_dereplicated.fasta -fastaout R1_no_singletons.fasta -minsize 2
```

#### DESCRIPTION: Clustering OTUs

-   INPUT: R1\_no\_singletons.fasta
-   OUTPUT: R1\_otus.fasta

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -cluster_otus R1_no_singletons.fasta -fulldp -otus R1_otus.fasta
```

#### DESCRIPTION: Chimera check

-   INPUT: R1\_otus.fasta
-   OUTPUT: R1\_otus\_no\_chimera\_ref.fasta

-   CURRENT DIRECTORY: clean

``` bash
vsearch -uchime_ref R1_otus.fasta -nonchimeras R1_otus_no_chimera_ref.fasta --db ../../reference_databases/uchime_sh_refs_dynamic_develop_985_01.01.2016.ITS1.fasta
```

#### DESCRIPTION: Renaming OTUs

-   INPUT: R1\_otus\_no\_chimera\_ref.fasta
-   OUTPUT: R1\_otus\_numbered.fasta

-   CURRENT DIRECTORY: clean

``` bash
python ../../code/python_scripts/fasta_number.py R1_otus_no_chimera_ref.fasta OTU > R1_otus_numbered.fasta
```

#### DESCRIPTION: Mapping the OTUs back to the raw reads

-   INPUT: R1\_otus\_numbered.fasta
-   OUTPUT: R1\_otu\_map.uc

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -usearch_global ../raw/R1_seqs.fna -db R1_otus_numbered.fasta -strand plus -id 0.97 -fulldp -uc R1_otu_map.uc -threads 2
```

#### DESCRIPTION: Creating an OTU table in .txt format

-   INPUT: R1\_otu\_map.uc
-   OUTPUT: R1\_otu\_table.txt

-   CURRENT DIRECTORY: clean

``` bash
python ../../code/python_scripts/uc2otutab_mod.py R1_otu_map.uc > R1_otu_table.txt
```

#### DESCRIPTION: Assigning taxonomy

-   INPUT: R1\_otus\_numbered.fasta
-   OUTPUT: R1\_taxonomy\_otutaxout.utax

-   CURRENT DIRECTORY: clean

``` bash
usearch8 -utax R1_otus_numbered.fasta -db ../../reference_databases/utax_db_22_08_16.db -strand both -utaxout R1_taxonomy_otutaxout.utax -utax_cutoff 0.8 -strand plus -threads 2
```

#### DESCRIPTION: Modify the taxonomy to the RDP taxonomy

-   INPUT: R1\_taxonomy\_otutaxout.utax
-   OUTPUT: R1\_taxonomy\_otutaxout.utax

-   CURRENT DIRECTORY: clean

``` bash
python ../../code/python_scripts/UTAX_to_RDP.py R1_taxonomy_otutaxout.utax
```

#### DESCRIPTION: Convert our OTU table to .biom format using QIIME

-   INPUT: R1\_otu\_table.txt
-   OUTPUT: R1\_otu\_table.biom

-   CURRENT DIRECTORY: clean

``` bash
biom convert -i R1_otu_table.txt -o R1_otu_table.biom --table-type="OTU table" --to-json
```

#### DESCRIPTION: Add taxonomy to the .biome file

-   INPUT:
-   .biom file: R1\_otu\_table.biom
-   tax file: R1\_taxonomy\_otutaxout\_\_RDP.txt
-   OUTPUT: R1\_otu\_table\_tax.biom

-   CURRENT DIRECTORY: clean

``` bash
biom add-metadata -i R1_otu_table.biom -o R1_otu_table_tax.biom --observation-metadata-fp R1_taxonomy_otutaxout__RDP.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence
```

-   CURRENT DIRECTORY: clean
-   CHANGE DIRECTORY: input\_files

``` bash
mkdir ../input_files
cd ../input_files
```

#### DESCRIPTION: Convert back to .json .biome format

-   INPUT: R1\_otu\_table\_tax.biom
-   OUTPUT: R1\_otu\_table\_tax\_json.biom

-   CURRENT DIRECTORY: input\_files

``` bash
biom convert -i ../clean/R1_otu_table_tax.biom -o R1_otu_table_tax_json.biom --table-type="OTU table" --to-json
```

#### DESCRIPTION: copy files into correct directories for R analysis

``` bash
cp ../clean/R1_taxonomy_otutaxout__RDP.txt . 
cp ../mapping.txt . 
```

Now you are ready to run R analysis on the files located in the input\_files directory
