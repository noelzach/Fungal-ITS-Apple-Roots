echo Script to do sequence preprocessing of fungal ITS 


#current directory code
cd ../data
#current directory data
validate_mapping_file.py -m mapping.txt -o mapping_output
mkdir clean
cd clean
#current directory clean
cutadapt -g CTTGGTCATTTAGAGGAAGTAA -e 0.01 --discard-untrimmed --match-read-wildcar ../raw/R1_seqs.fastq > R1_stripped.fastq
vsearch -fastq_stats R1_stripped.fastq -log R1_stats_result_vsearch.txt
usearch8 -fastq_filter R1_stripped.fastq -fastq_maxee 0.1 -fastq_trunclen 180 -fastq_maxns 0 -fastaout R1_filtered.fasta -fastqout R1_filtered.fastq
usearch8 -fastq_filter R1_filtered.fastq -fastq_stripleft 44 -fastaout R1_trimmed.fasta -fastqout R1_trimmed.fastq
mkdir fastqc
cd fastqc 
#current directory fastqc
fastqc ../R1_trimmed.fastq
cd .. 
#current directory clean
usearch8 -derep_fulllength R1_trimmed.fasta -fastaout R1_dereplicated.fasta -sizeout -strand both
usearch8 -sortbysize R1_dereplicated.fasta -fastaout R1_no_singletons.fasta -minsize 2
usearch8 -cluster_otus R1_no_singletons.fasta -fulldp -otus R1_otus.fasta
vsearch -uchime_ref R1_otus.fasta -nonchimeras R1_otus_no_chimera_ref.fasta --db ../../reference_databases/uchime_sh_refs_dynamic_develop_985_01.01.2016.ITS1.fasta
python ../../code/python_scripts/fasta_number.py R1_otus_no_chimera_ref.fasta OTU > R1_otus_numbered.fasta
usearch8 -usearch_global ../raw/R1_seqs.fna -db R1_otus_numbered.fasta -strand plus -id 0.97 -fulldp -uc R1_otu_map.uc -threads 2
python ../../code/python_scripts/uc2otutab_mod.py R1_otu_map.uc > R1_otu_table.txt
usearch8 -utax R1_otus_numbered.fasta -db ../../reference_databases/utax_db_22_08_16.db -strand both -utaxout R1_taxonomy_otutaxout.utax -utax_cutoff 0.8 -strand plus -threads 2
python ../../code/python_scripts/UTAX_to_RDP.py R1_taxonomy_otutaxout.utax
biom convert -i R1_otu_table.txt -o R1_otu_table.biom --table-type="OTU table" --to-json
biom add-metadata -i R1_otu_table.biom -o R1_otu_table_tax.biom --observation-metadata-fp R1_taxonomy_otutaxout__RDP.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence
mkdir ../input_files
cd ../input_files
#current directory input_files
biom convert -i ../clean/R1_otu_table_tax.biom -o R1_otu_table_tax_json.biom --table-type="OTU table" --to-json
cp ../clean/R1_taxonomy_otutaxout__RDP.txt . 
cp ../mapping.txt . 