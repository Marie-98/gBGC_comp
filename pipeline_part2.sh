cd ~/Murinae/Hydromyini/

### PHYLOGENY ################################################################################

## create a file of concatenation of all exons

# concatenate exons sequences per specie (one fasta file per species that contains all exons sequences)
mkdir Aligned_concatenated_sequences
cd Aligned_concatenated_sequences

cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt | parallel -j20 ~/scripts/script_concatenate_exons.sh {} Hydromyini

# reformate fasta files
for sp in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt); do sed -i ':a;N;/>/!s/\n//;ta;P;D' ${sp}_all_seq_NT_aligned_final.fasta ; done

for sp in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt); do sed -i 's/!/N/g' ${sp}_all_seq_NT_aligned_final.fasta ; done

mkdir wt_chevrons
cd wt_chevrons

for sp in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt); do grep -v ">" ../${sp}_all_seq_NT_aligned_final.fasta > ${sp}_all_seq_NT_aligned_final.fasta  ; done

# concatenate sequences of all species in one fasta file
cd ../
mkdir All_seq
cd All_seq

for sp in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt) ; do echo ">${sp}" ; cat ~/Murinae/Hydromyini/Aligned_concatenated_sequences/wt_chevrons/${sp}_all_seq_NT_aligned_final.fasta ; done >> Hydromyini_all_seq_NT_aligned_final.fasta

## make the phylogeny on the concatenation of all exons
cd ~/Murinae/Hydromyini/
mkdir Phylogeny
cd Phylogeny

# you need to install IQTREE-2

iqtree2 -s ~/Murinae/Hydromyini/Aligned_concatenated_sequences/All_seq/Hydromyini_all_seq_NT_aligned_final.fasta  -nt 20 -m GTR -bb 1000 2>error_iqtree.txt >out_iqtree.txt

mv ~/Murinae/Hydromyini/Aligned_concatenated_sequences/All_seq/Hydromyini_all_seq_NT_aligned_final.fasta.treefile ~/Murinae/Hydromyini/Phylogeny/

## Root the phylogeny
Rscript ~/scripts/script_root_tree.R ~/Murinae/Hydromyini/Phylogeny/Hydromyini_all_seq_NT_aligned_final.fasta.treefile Haeromys_minahassae Hydromyini_rooted_tree.treefile

## make one phylogeny per gene by pruning the global phylogeny
	# keep only the common species between all exons of a same gene
cd ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/alignments_wt_ref_filtered/alignments_wt_ref_filtered_wt_stop/

find -maxdepth 1 -type f -name "*.fasta" | cut -f1,2,3 -d"_" | cut -f2 -d"/" > Genes_exons_pos_list.csv
sed -i 's/_/\t/g' Genes_exons_pos_list.csv

cut -f1 Genes_exons_pos_list.csv | sort -u > Genes_list.txt

cd ~/Murinae/Hydromyini/Aligned_Sequences/
mkdir FINAL
cd FINAL

cat ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/alignments_wt_ref_filtered/alignments_wt_ref_filtered_wt_stop/Genes_list.txt | parallel -j20 python3 ~/scripts/seq_com_exons_same_gene.py {} ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/alignments_wt_ref_filtered/alignments_wt_ref_filtered_wt_stop/Genes_exons_pos_list.csv ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/alignments_wt_ref_filtered/alignments_wt_ref_filtered_wt_stop/ 2>error_seq_com_genes.txt >out_seq_com_genes.txt

mkdir not_used
for align in *_FINAL_align_NT.fasta ; do 
if (($(grep -c ">" $align) < 3)) ; then mv $align not_used ; fi 
done

	# make one tree per gene
find -maxdepth 1 -type f -name "*.fasta" | cut -f1 -d'_' | cut -f2 -d'/' | sort -u > Final_genes_list.txt
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d'_' | sort -u > Final_exons_list.txt
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d'/' | cut -f1,2,3 -d'_' > Final_genes_exons_pos_list.txt
sed -i 's/_/\t/g' Final_genes_exons_pos_list.txt

cd ~/Murinae/Hydromyini/Phylogeny/
mkdir gene_tree
cd gene_tree

cat ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/Final_genes_list.txt | parallel -j20 ~/scripts/script_arbre_gene.sh Hydromyini {} 2> error_gene_tree.txt >out_gene_tree.txt


### SUPPLEMENTARY FILTER OF PARALOGS #########################################################

cd ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/

for file in *.fasta ; do len=$(seqkit fx2tab $file -n -l | head -1 | cut -f2) ; echo -e "${file}\t${len}" ; done > align_len.csv

cd ~/Murinae/Hydromyini/
mkdir Sup_filter_paralogs
cd Sup_filter_paralogs

   ## !! ## chek that you have the python module ete3 installed ##
cat ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/Final_exons_list.txt | parallel -j20 ~/scripts/script_filtre_sup_paralogues.sh {} Hydromyini > tab_lg_tree_for_sup_filter_paralogues.csv

~/scripts/script_list_paralog_exons.R tab_lg_tree_for_sup_filter_paralogues.csv


### SUBSTITUTIONS MAPPING ####################################################################

cd ~/Murinae/Hydromyini/
mkdir Substitutions_mapping

### make the list of exons with at least 10 species (to make the mapping only on exons with >=10 sp)

cd ~/Murinae/Hydromyini/Aligned_Sequences/FINAL
for exon in $(cat Final_exons_list.txt) ; do gene=$(grep $exon Final_genes_exons_pos_list.txt | cut -f1) ; pos=$(grep $exon Final_genes_exons_pos_list.txt | cut -f3) ; nb=$(grep -c ">" ${gene}_${exon}_${pos}_FINAL_align_NT.fasta) ; echo -e "${gene}\t${exon}\t${nb}" ; done > Genes_exons_nbsp.csv

LC_ALL=en_US awk -v min=10 '{if ($3>=min) {print $0}}' Genes_exons_nbsp.csv > Genes_exons_nbsp_min10.csv

cut -f2 Genes_exons_nbsp_min10.csv >  list_exons_min10sp.txt

### make the substitutions mapping

cd ~/Murinae/Hydromyini/Substitutions_mapping/
mkdir output

# you need to install bppml and mapnh (put them in your bin)

cat ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/list_exons_min10sp.txt | parallel -j20 ~/scripts/script_mapNH_bpp.sh Hydromyini {} 2> error_sortie_mapNH.txt > out_mapNH.txt

find -type f -name "*.dnd_1" | cut -f1 -d"_" | cut -f2 -d"/" > list_exons_mapping_ok.txt

cd output
mkdir correct_file_names
cat ../list_exons_mapping_ok.txt | parallel -j20 ~/scripts/transform_mapping_file_name.sh {}


### DETECTION OF gBGC EPISODES ####################################################################

#######
## 1 ## Calculate statistics ###
#######

### Distribution of synonymous substitutions

cd ~/Murinae/Hydromyini/
mkdir Detecting_gBGC
cd Detecting_gBGC
mkdir Distrib_subst_S
cd Distrib_subst_S

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 Rscript ~/scripts/script_nb_subst.R {} ~/Murinae/Hydromyini/Substitutions_mapping/output/correct_file_names/ {}_stat_nb_subst.csv 2> errors_script_nb_subst.txt > out_script_nb_subst.txt

### Moran's I

# list exons > 1000bp (to calculate Moran's I only on them)
cd ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/

for exon in $(cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt) ; do len=$(grep $exon align_len.csv | cut -f2); echo -e "${exon}\t${len}" ; done > exons_len.csv

LC_ALL=en_US awk -v min=1000 '{if ($2>=min) {print $1}}' exons_len.csv > list_exons_sup1000bp.txt 

cd ~/Murinae/Hydromyini/Detecting_gBGC
mkdir I_Moran
cd I_Moran

cat ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/list_exons_sup1000bp.txt | parallel -j20 python3 ~/scripts/script_I_Moran.py {} ~/Murinae/Hydromyini/Substitutions_mapping/output/correct_file_names/ I_Moran_{}.csv 2> errors_script_I_Moran.txt > out_script_I_Moran.txt


#######
## 2 ## Randomisations of statistics ###
#######

cd ~/Murinae/Hydromyini/
mkdir Randomisations
cd Randomisations
mkdir Br_lg_tables
mkdir Distrib_subst
mkdir I_Moran

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_simul_stat_V2.sh {} Hydromyini 1000 > out_simul.txt 2> error_simul.txt

mkdir used
for exon in $(cat ~/Murinae/Aligned_Sequences/Hydromyini/FINAL/list_exons_sup1000bp.txt) ; do
mv ${exon}_I_simul.csv used ;
done
rm *.csv
cd used
mv *.csv ..
cd ..
rm -r used


#######
## 3 ## Detection of gBGC episodes ###
#######

cd ~/Murinae/Hydromyini/Detecting_gBGC/
mkdir Detection
cd Detection
mkdir signif_stat
cd signif_stat

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_signif_stat_V2.sh {} Hydromyini 1000 >out_commande_signif_stat.txt 2> error_commande_signif_stat.txt

cd ~/Murinae/Hydromyini/Substitutions_mapping/
for ex in $(cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt) ; do len=$(grep ${ex} ~/Murinae/Hydromyini/Aligned_Sequences/FINAL/exons_len.csv | cut -f2) ; gene=$(grep ${ex} ~/Murinae/Reference_exons/Gene_exons_ref_list.txt | cut -f1) ; pos=$(grep -w ${ex} ~/Murinae/Reference_exons/Gene_exons_ref_list.txt | cut -f3) ; echo -e "${gene}\t${ex}\t${pos}\t${len}" ; done > gene_exon_lenseq_mapping_ok.txt

cd ~/Murinae/Hydromyini/Detecting_gBGC/Detection/
mkdir episodes
cd episodes

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_detection_episodes.sh {} Hydromyini 1000 995 > out_detection_episodes.txt 2> errors_detection_episodes.txt


### DISTRIBUTION & RANDOMISATION OF NON-SYNONYMOUS SUBSTITUTIONS ####################################

#######
## 1 ## Distribution ###
#######

cd ~/Murinae/Hydromyini/
mkdir Impact_episodes_NS
cd Impact_episodes_NS

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 Rscript ~/scripts/script_nb_subst_NS.R {} ~/Murinae/Hydromyini/Substitutions_mapping/output/correct_file_names/ {}_stat_nb_subst_NS.csv 2> errors_script_nb_subst_NS.txt > out_script_nb_subst_NS.txt


#######
## 2 ## Randomisation ###
#######

cd ~/Murinae/Hydromyini/Randomisations/
mkdir Distrib_NS
cd Distrib_NS

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_simul_distrib_subst_NS.sh {} Hydromyini 1000 > sortie_simul_NS.txt 2> error_simul_NS.txt


### CALCULATE GC CONTENT #########################################################################

cd ~/Murinae/Hydromyini/
mkdir GC_content
cd GC_content

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_GC_content.sh {} Hydromyini >out_GC_content.txt 2>errors_GC_content.txt


### BRANCH RELATIONS #############################################################################

cd ~/Murinae/Hydromyini/Substitutions_mapping/

for exon in $(cat list_exons_mapping_ok.txt) ; do
gene=$(grep $exon ~/Murinae/Reference_exons/Gene_exons_ref_list.txt | cut -f1) ; echo -e "${gene}\t${exon}" ; done > list_genes_exons_mapping_ok.txt

cut -f1 list_genes_exons_mapping_ok.txt | sort -u > list_genes_mapping_ok.txt

cd ~/Murinae/Hydromyini/Detecting_gBGC/Detection/
mkdir relations_branches
cd relations_branches

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_genes_mapping_ok.txt | parallel -j20 ~/scripts/script_num_nodes_tree.sh {} Hydromyini 2> error_rel_br.txt > out_rel_br.txt


### FINAL TABLE WITH INFORMATIONS FOR ALL EXONS FOR GENERAL RESULTS ###############################

cd ~/Murinae/Hydromyini/
mkdir Final_data
cd Final_data

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_genes_mapping_ok.txt | parallel -j20 ~/scripts/script_table_recap_V2.py {} ~/Murinae/Hydromyini/Substitutions_mapping/gene_exon_lenseq_mapping_ok.txt ~/Murinae/Hydromyini/Detecting_gBGC/Detection/relations_branches/{}_rel_br.csv ~/Murinae/Hydromyini/Randomisations/Br_lg_tables/{}_br_len.csv ~/Murinae/Hydromyini/Detecting_gBGC/Detection/episodes/{}_episodes_detection.csv ~/Murinae/Hydromyini/Detecting_gBGC/Distrib_subst_S/ ~/Murinae/Hydromyini/Detecting_gBGC/Detection/signif_stat/ ~/Murinae/Hydromyini/Impact_episodes_NS/ ~/Murinae/Hydromyini/GC_content/ {}_tab_recap.csv > sortie_tab_recap.txt 2> error_tab_recap.txt

for file in *_recap.csv ; do sed 1d $file ; done > tab_recap_all_genes.csv

ls *_recap.csv | cut -f1 -d'_' > list_genes.txt


### FINAL TABLES FOR SYNONYMOUS DATA ANALYSES (OVERFLOWING OF EPISODES) ############################


### Randomisation of synonymous substitutions by blocking episode branches ###

cd ~/Murinae/Hydromyini/Randomisations/
mkdir Distrib_S_blockbrep
cd Distrib_S_blockbrep

cat ~/Murinae/Hydromyini/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_simul_distrib_subst_S_blokbrep.sh {} Hydromyini 1000 > sortie_simul_S.txt 2> error_simul_S.txt


### Make tables ###

cd ~/Murinae/Hydromyini/
mkdir Final_data_analyse_Synonymous
cd Final_data_analyse_Synonymous

cat ~/Murinae/Hydromyini/Final_data/list_genes.txt | parallel -j20 python3 ~/scripts/analyse_random_S_V3.py {} ~/Murinae/Hydromyini/Final_Data/ ~/Murinae/Hydromyini/Randomisations/Distrib_S_blockbrep/ ~/Murinae/Hydromyini/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomS.txt 2>error_analyse_randomS.txt

~/scripts/script_concat_tab_analyse_randomS_V3.sh # concatenate the tables of all genes


### FINAL TABLES FOR COMPENSATION ANALYSES ########################################################

cd ~/Murinae/Hydromyini/
mkdir Final_data_analyse_NonSynonymous_Compensation
cd Final_data_analyse_NonSynonymous_Compensation

cat ~/Murinae/Hydromyini/Final_data/list_genes.txt | parallel -j20 python3 ~/scripts/analyse_random_NS_V3.py {} ~/Murinae/Hydromyini/Final_data/ ~/Murinae/Hydromyini/Randomisations/Distrib_NS/ ~/Murinae/Hydromyini/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

~/scripts/script_concat_tab_analyse_randomNS_V3.sh

# control upstream branches :
mkdir control_br_asc
cd control_br_asc/

cat ~/Murinae/Hydromyini/Donnees_finales/list_genes.txt | parallel -j20 python3 ~/scripts/analyse_random_NS_test_br_asc.py {} ~/Murinae/Hydromyini/Donnees_finales/ ~/Murinae/Hydromyini/Randomisations/distrib_NS/ ~/Murinae/Hydromyini/Sup_filter_paralogs/paralog_exons.txt >sortie_analyse_randomNS_test_br_asc.txt 2>error_analyse_randomNS_test_br_asc.txt

~/scripts/script_concat_tab_analyse_randomNS_test_br_asc.sh


### SORTING EPISODE TYPE (punctual, two-branches, or all-exons) ####################################

cd ~/Murinae/Hydromyini/
mkdir Sort_episodes_types
cd Sort_episodes_types

### make table for the sorting ###
cat ~/Murinae/Hydromyini/Final_data/list_genes.txt | parallel -j20 python3 ~/scripts/script_tab_for_sort_ep.py {} ~/Murinae/Hydromyini/Final_data/ ~/Murinae/Hydromyini/Sup_filter_paralogs/paralog_exons.txt 2>error.txt >out.txt

for file in *_tab_for_sort_ep.csv ; do sed 1d $file ; done > tab_for_sort_ep_all_genes_hydromyini.csv
## !!! ## AJOUTER TITRE COLONNES ###

### make the sorting ####
~/bin/sort_ep_type/MLsort tab_for_sort_ep_all_genes_hydromyini.csv optfile.txt posterior_optimum_hydromyini.csv



### DISTANCE BETWEEN DELETERIOUS MUTATIONS AND COMPENSATORY MUTATIONS ############################

cd ~/Murinae/Hydromyini/
mkdir Dist_delet_comp_mut
cd Dist_delet_comp_mut
mkdir Data
mkdir Randomisations

### make table of distances between deleterious and compensatory mutations ###

cd Data/

cat ~/Murinae/Hydromyini/Final_Data/list_genes.txt | parallel -j20 python3 ~/scripts/dist_comp.py {} ~/Murinae/Hydromyini/Final_data/ ~/Murinae/Hydromyini/Sup_filter_paralogs/paralog_exons.txt ~/Murinae/Hydromyini/Substitutions_mapping/output/correct_file_names/ >out_dist_comp.txt 2>error_dist_comp.txt

for file in *_tab_dist_comp.csv ; do sed 1d $file ; done > tab_dist_comp_all_genes.csv


### randomise the positions of compensatory mutations and therefore of the distances, and get the number of 1-1 cases with a distance less than X (3 in the article) ###

cd ../Randomisations/

python3 ~/scripts/random_dist_comp.py ~/Murinae/Hydromyini/Dist_delet_comp_mut/tab_dist_comp_all_genes.csv 3 1000 2>error_random_dist_comp.txt >out_random_dist_comp.txt

