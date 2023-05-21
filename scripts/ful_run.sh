#!/bin/bash

# Prompt user for keyword1_phenotype and keyword2_metab
echo "Enter keyword1_phenotype:"
read keyword1_phenotype
echo "Enter keyword2_metab:"
read keyword2_metab

# Get cluster value where row=keyword2_metab from metabolite_clusters.csv
CLUSTER=$(awk -F ',' -v k="$keyword2_metab" '$1==k{print $2}' "./$keyword1_phenotype/metabolite_clusters.csv")

# Get bacteria values where cluster=CLUSTER from clusters_int.csv
awk -F ',' -v c="$CLUSTER" '$2==c{print $1}' "./$keyword1_phenotype/clusters_int.csv" > metab_genomes_tmp.txt

# Combine metab_genomes_tmp.txt and golden_standart.txt into full_genome_set_tmp.txt
cat "./$keyword1_phenotype/golden_standart.txt" metab_genomes_tmp.txt > full_genome_set_tmp.txt

# Get idx values where name=values of full_genome_set_tmp.txt from cnode_name_to_idx.csv
awk -F ',' 'NR==FNR{a[$1]=$2;next}($1 in a){print a[$1]}' full_genome_set_tmp.txt "./$keyword1_phenotype/cnode_name_to_idx.csv" > full_genome_set_ids_tmp.txt

# Run cnode_predict.sh with model.b and full_genome_set_ids_tmp.txt as inputs
bash ./cnode_predict.sh "./$keyword1_phenotype/model.b" full_genome_set_ids_tmp.txt
