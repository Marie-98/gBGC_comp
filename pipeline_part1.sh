cd ~/
mkdir Murinae
cd Murinae/
mkdir Hydromyini
cd Hydromyini/

### READS ####################################################################################

#######
## 1 ## Load reads with SRATools ###
#######

# need a .txt file with accession numbers (1 number per line) --> in the "in_files" folder

mkdir Data
cd Data # put the file of accession numbers here
mkdir reads
cd reads

for i in $(cat ~/Murinae/Hydromyini/Data/Hydromyini_accession_numbers.txt) ; do fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $i ; done

#######
## 2 ## Correct reads names ###
#######

cd ~/Murinae/Hydromyini/Data/
mkdir reads_correct_names
cd reads_correct_names

for ind_id in $(cat ~/Murinae/Hydromyini/Data/Hydromyini_accession_numbers.txt) ; do sed "s/\(^@\)\([0-9]*\/1\)$/@${ind_id}_\2/g" ~/Murinae/Hydromyini/Data/reads/${ind_id}_1.fastq > ${ind_id}_1.fastq ; done

for ind_id in $(cat ~/Murinae/Hydromyini/Data/Hydromyini_accession_numbers.txt) ; do sed "s/\(^@\)\([0-9]*\/2\)$/@${ind_id}_\2/g"  ~/Murinae/Hydromyini/Data/reads/${ind_id}_2.fastq > ${ind_id}_2.fastq ; done

## !! # now you should suppress the initial reads files (uncorrected names) to not take too much memory

#######
## 3 ## Clean reads ###
#######

# need to install Trimmomatic (put it in your bin)

cd ~/Murinae/Hydromyini/Data/
mkdir cleaned_reads
cd cleaned_reads

cat ~/Murinae/Hydromyini/Data/Hydromyini_accession_numbers.txt | parallel -j 20 java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 3 -trimlog {}.trimmomatic.log ~/Murinae/Hydromyini/Data/reads_correct_names/{}_1.fastq ~/Murinae/Hydromyini/Data/reads_correct_names/{}_2.fastq {}.R1.clean.paired.fastq.gz {}.R1.clean.unpaired.fastq.gz {}.R2.clean.paired.fastq.gz {}.R2.clean.unpaired.fastq.gz ILLUMINACLIP:/media/newvol/mriffis/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#######
## 4 ## Pool reads ###
#######

# need a .csv file with col 1 : species, col 2 : id of individus (accession numbers), with 1 line per individu --> in the "in_files" folder
# need a .txt file with list of species (1 specie per line) --> in the "in_files" folder

cd ~/Murinae/Hydromyini/Data/ # put the 2 files here
mkdir cleaned_pooled_reads
cd cleaned_pooled_reads

for specie in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt); do for indiv in $(grep $specie ~/Murinae/Hydromyini/Data/recap_sequencing_Hydromyini.csv | cut -f 2 -d','); do cat ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R1.clean.paired.fastq.gz ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R1.clean.unpaired.fastq.gz >> ${specie}_all_R1.clean.fastq.gz; done; done

for specie in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt); do for indiv in $(grep $specie ~/Murinae/Hydromyini/Data/recap_sequencing_Hydromyini.csv | cut -f 2 -d','); do cat ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R2.clean.paired.fastq.gz ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R2.clean.unpaired.fastq.gz >> ${specie}_all_R2.clean.fastq.gz; done; done


### ASSEMBLY #################################################################################

# need to install Trinity

cd ~/Murinae/Hydromyini/
mkdir Assembly
cd Assembly

# do assembly with Trinity
cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt | parallel -j20 ~/scripts/script_trinity_assemblage.sh {} ~/Murinae/Hydromyini/Data/cleaned_pooled_reads/ >results_assembly_Hydromyini.txt 2>errors_assembly_Hydromyini.txt

# !! # Be careful Trinity creates a lot of files in addition of the assembly files, that take non negligible memory ; you can suppress them and only keep the assembly files. 

# get some stats on the assembly
for sp in $(cat ~/Murinae/Hydromyini/Data/cleaned_pooled_reads/list_species_Hydromyini.txt) ; do
nb_contigs=$(grep -c ">" trinity_OUT_${sp}.Trinity.fasta)
echo -e "${sp}\t${nb_contigs}"
done >> nb_contigs.csv

for sp in $(cat ~/Murinae/Hydromyini/Data/cleaned_pooled_reads/list_species_Hydromyini.txt) ; do
echo ${sp}
seqkit stats trinity_OUT_${sp}.Trinity.fasta -a
done >> stats_contigs.txt


### SET OF CODING ORTHOLOGOUS EXONS ##########################################################

# need the fasta file with the informations and sequences of the mouse reference exons --> in the "in_files" folder
# need the fasta files with the informations and UTR5' and UTR3' sequences of the reference exons --> in the "in_files" folder

cd ~/Murinae/
mkdir Reference_exons
cd Reference_exons # put here the fasta files

#######
## 1 ## Retrieve the UTR sequences of the reference exons ###
#######

# need to install BLAST tool

# blast the UTR sequences against the exons sequences
~/scripts/script_blast_UTR_exons.sh ~/Murinae/Reference_exons/Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.txt ~/Murinae/Reference_exons/exons_5UTR.fasta ~/Murinae/Reference_exons/exons_3UTR.fasta 20

# list of the reference exons
grep ">" Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.txt | cut -f2 -d'|' > list_exons_ref_all.txt

# need to install the seqkit tool

# filter to get list of exons > 100 bp
seqkit fx2tab Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.txt -l -n | LC_ALL=en_US awk -v taillemin=99 '{if ($2>taillemin) {print $0}}' | cut -f2 -d'|' | sort -u > List_exons_ref_lgmin100_sorted.txt

# list of the reference genes
grep ">" Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.txt | cut -f2 -d'>' | cut -f1 -d'|' | sort -u > Mouse_gene_ref_list_all.txt
	
# table of exons positions in the gene
cat list_exons_ref_all.txt | parallel -j20 ~/scripts/get_exons_position_table.sh {}

# get the positions of exons at the extremities of the gene (start or end)
	# end
for gene in $(cat Mouse_gene_ref_list_all.txt) ; do
nb_exons=$(grep -w "${gene}" gene_exon_position_ref_all.txt | sort -r -k3 | head -1 | cut -f3) ;
echo -e "${gene}\t${nb_exons}" >> ref_gene_nb_exons_all.txt ;
done

for exon in $(cat list_exons_ref_all.txt) ; do
gene=$(grep -w "${exon}" gene_exon_position_ref_all.txt | cut -f1) ;
exon_position=$(grep -w "${exon}" gene_exon_position_ref_all.txt | cut -f3) ;
nb_exons=$(grep -w "${gene}" ref_gene_nb_exons_all.txt | cut -f2) ;
if (( "$exon_position" == "$nb_exons" )) ; then echo -e "${exon}\tEND" >> list_exons_end_gene_all.txt ; fi ;
done

	# start
for exon in $(cat list_exons_ref_all.txt) ; do
exon_position=$(grep -w "${exon}" gene_exon_position_ref_all.txt | cut -f3) ;
if (( $exon_position == 1 )) ; then echo -e "${exon}\tSTART" >> list_exons_start_gene_all.txt ; fi ;
done

	# all
cat list_exons_end_gene_all.txt list_exons_start_gene_all.txt > list_all_exons_extremities_genes_positions.txt

# list of exons at extremities of the genes
cut -f1 list_all_exons_extremities_genes_positions.txt | sort -u > list_all_exons_extremities_genes.txt

# get exons sequences without the UTR sequences
mkdir coding_exons
cd coding_exons

cat ../list_all_exons_extremities_genes.txt | parallel -j20 ~/scripts/script_get_exon_seq_sans_utr_V3.sh {} ../list_all_exons_extremities_genes_positions.txt ../blast_UTR_5_vs_exons.txt ../blast_UTR_3_vs_exons.txt ../Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.fasta 2 >> Mouse_exons_ref_extremities_genes_coding_seq.fasta 2>error.txt

# filter exons for which the length of blast is < of the length of the UTR sequence
awk '{if ($2==$3) {print $1}}' utr_blast_lengths.csv | sort -u > List_exons_ref_without_utr_filtered.txt

for exon in $(cat List_exons_ref_without_utr_filtered.txt); do
complete_name=$(grep $exon ../Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.fasta)
echo -e "${complete_name}" >> List_exons_ref_without_utr_filtered_complete_name.txt;
done

sed -i 's/>//g' List_exons_ref_without_utr_filtered_complete_name.txt

seqkit grep -n -f List_exons_ref_without_utr_filtered_complete_name.txt Mouse_exons_ref_extremities_genes_coding_seq.fasta -o Exons_ref_extremities_genes_coding_seq_filtered.fasta

# filter to keep only sequences > 100 bp
seqkit fx2tab Exons_ref_extremities_genes_coding_seq_filtered.fasta -l -n | LC_ALL=en_US awk -v taillemin=99 '{if ($2>taillemin) {print $0}}' | cut -f1 > List_exons_ref_without_utr_filtered_lgmin100.txt

seqkit grep -n -f List_exons_ref_without_utr_filtered_lgmin100.txt Exons_ref_extremities_genes_coding_seq_filtered.fasta -o Exons_ref_extremities_genes_coding_seq_filtered_lgmin100.fasta


# get the list and sequences of references exons not at the extremities of the genes
cd ../

sort -u list_all_exons_extremities_genes.txt > list_all_exons_extremities_genes_sorted.txt

comm -23 List_exons_ref_lgmin100_sorted.txt list_all_exons_extremities_genes_sorted.txt > List_exons_ref_lgmin100_not_extremities.txt

for exon in $(cat List_exons_ref_lgmin100_not_extremities.txt) ; do
complete_name=$(grep $exon Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.fasta)
echo -e "${complete_name}" >> List_exons_ref_lgmin100_not_extremities_complete_name.txt;
done

sed -i 's/>//g' List_exons_ref_lgmin100_not_extremities_complete_name.txt

seqkit grep -n -f List_exons_ref_lgmin100_not_extremities_complete_name.txt Mouse_Exons_ref_geneID_exonID_transcriptID_CDSlength_exonrank.fasta -o Mouse_Exons_ref_not_extremities_lgmin100.fasta

# make the final fasta file of reference exons : all exons ref, without the UTR sequences and longer than 100 bp
cat Mouse_Exons_ref_not_extremities_lgmin100.fasta coding_exons/Exons_ref_extremities_genes_coding_seq_filtered_lgmin100.fasta > Mouse_exons_ref_coding_seq_lgmin100.fasta

grep ">" Mouse_exons_ref_coding_seq_lgmin100.fasta | cut -f2 -d'|' > List_mouse_exons_ref_coding_seq_lgmin100.txt

for exon in $(cat List_mouse_exons_ref_coding_seq_lgmin100.txt) ; do
gene=$(grep $exon Mouse_exons_ref_coding_seq_lgmin100.fasta | head -1 | cut -f2 -d'>' | cut -f1 -d'|')
pos=$(grep $exon Mouse_exons_ref_coding_seq_lgmin100.fasta | head -1 | cut -f5 -d'|')
echo -e "${gene}\t${exon}\t${pos}" >> Gene_exons_ref_list.txt ;
done

#######
## 2 ## Blast of reference exons against assembled contigs ###
#######

cd ~/Murinae/Hydromyini/
mkdir Blast
cd Blast
mkdir blast_exon_vs_contigs
cd blast_exon_vs_contigs

# make blast #
cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt | parallel -j5 ~/scripts/script_blast_exons_vs_contigs.sh ~/Murinae/Reference_exons/Mouse_exons_ref_coding_seq_lgmin100.fasta Hydromyini {} 99 85 0.5 4 >results_blast_exons_vs_contigs_Hydromyini.txt 2>errors_blast_exons_vs_contigs_Hydromyini.txt

# suppress the database
cd ~/Murinae/Hydromyini/Assembly/
rm *.ndb
rm *.nhr
rm *.nin
rm *.not
rm *.nsq
rm *.ntf
rm *.nto

#######
## 3 ## Reciprocal blast ###
#######

cd ~/Murinae/Hydromyini/Blast/
mkdir reciprocal_blast
cd reciprocal_blast

makeblastdb -in ~/Murinae/Reference_exons/Mouse_exons_ref_coding_seq_lgmin100.fasta -dbtype nucl

for sp in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt) ; do 
mkdir ${sp}
cd ${sp}
cat ~/Murinae/Hydromyini/Blast/blast_exon_vs_contigs/${sp}/exons_list_blast_${sp}.txt | parallel -j20 ~/scripts/script_blast_reciproque_V2.sh {} ~/Murinae/Hydromyini/Blast/blast_exon_vs_contigs/${sp}/blast_exons_vs_contigs_filtre_${sp}.txt ~/Murinae/Hydromyini/Assembly/trinity_OUT_${sp}.Trinity.fasta ${sp} 1 >out_rec_blast_${sp}.txt 2>warnings_rec_blast_${sp}.txt
cd ../
done


### ALIGNMENT ################################################################################

# concatenate sequences of species in order to have one fasta file per exon that contain the sequences of all the species that have this exon
cd ~/Murinae/Hydromyini/
mkdir Sequences
cd Sequences

cat ~/Murinae/Reference_exons/List_mouse_exons_ref_coding_seq_lgmin100.txt | parallel -j20 ~/scripts/script_concatener_especes_exons.sh {} Hydromyini

# makes fasta files for reference exons
cd ~/Murinae/Reference_exons/
mkdir seq_ref_for_aln
cd seq_ref_for_aln

cat ~/Murinae/Reference_exons/List_mouse_exons_ref_coding_seq_lgmin100.txt | parallel -j20 ~/scripts/script_seq_ref.sh {}

for file in *_seqname.txt ; do rm $file ; done

# align sequences based on the reference sequences
cd ~/Murinae/Hydromyini/Sequences
	# make exons list
find -type f -name "*.fasta" | cut -f2 -d'_' > Exons_list.txt

cd ~/Murinae/Hydromyini/
mkdir Aligned_Sequences
cd Aligned_Sequences

# you need to install OMM-MACSE alignment tool (put it in your bin)

	# align sequences
cat ~/Murinae/Hydromyini/Sequences/Exons_list.txt | parallel -j20 ~/scripts/script_alignement_omm_macse.sh {} Hydromyini > out_alignment_Hydromyini.txt 2> errors_alignment_Hydromyini.txt
	# move all alignments in the same directory
mkdir alignments

cat ~/Murinae/Hydromyini/Sequences/Exons_list.txt | parallel -j20 ~/scripts/script_mv_align_omm.sh {} Hydromyini
	# make list of aligned exons
cd alignments

find -type f -name "*.aln" | cut -f2 -d'_' > list_exons_align.txt
	# retrieve the reference sequence from the alignments
cd ~/Murinae/Reference_exons
grep ">" Mouse_exons_ref_coding_seq_lgmin100.fasta > Complete_names_Mouse_exons_ref_coding_seq_lgmin100.txt
sed -i 's/>//g' Complete_names_Mouse_exons_ref_coding_seq_lgmin100.txt

cd ~/Murinae/Hydromyini/Aligned_Sequences/alignments/
mkdir alignments_wt_ref
cd alignments_wt_ref

cat ~/Murinae/Hydromyini/Aligned_Sequences/alignments/list_exons_align.txt | parallel -j20 ~/scripts/script_remove_ref_from_align.sh {} Hydromyini


### FILTER ALIGNMENTS ########################################################################

#######
## 1 ## Mapping reads ###
#######

# mapping reads on contigs in order to calcul mean coverage depth and heterozygosity

# for species with several individus, make one reads file per individu
cd ~/Murinae/Hydromyini/Data/cleaned_pooled_reads/

cut -f2 -d',' ~/Murinae/Hydromyini/Data/recap_sequencing_Hydromyini.csv > list_indiv_Hydromyini.txt

for indiv in $(cat list_indiv_Hydromyini.txt) ; do cat ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R1.clean.paired.fastq.gz ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R1.clean.unpaired.fastq.gz >> ${indiv}_all_R1.clean.fastq.gz ; done

for indiv in $(cat list_indiv_Hydromyini.txt) ; do cat ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R2.clean.paired.fastq.gz ~/Murinae/Hydromyini/Data/cleaned_reads/${indiv}.R2.clean.unpaired.fastq.gz >> ${indiv}_all_R2.clean.fastq.gz ; done

# mapping reads
cd ~/Murinae/Hydromyini/
mkdir Mapping
cd Mapping

# you need to install BWA tool

cat ~/Murinae/Hydromyini/Data/cleaned_pooled_reads/list_indiv_Hydromyini.txt | parallel -j4 ~/scripts/script_mapping.sh 5 {} Hydromyini >out_mapping_Hydromyini.txt 2>errors_mapping_Hydromyini.txt

# suppress the index files
cd ~/Murinae/Hydromyini/Assembly/
rm *.amb
rm *.ann
rm *.bwt
rm *.pac
rm *.sa
rm *.fai

## !! # Be careful the mapping of reads creates many intermediate files that take a lot of memory (around 2T for these data), you can keep only the "*_sorted_bam.bam" files and suppress all the others

#######
## 2 ## Heterozygosity and mean coverage depth ###
#######

## heterozygosity
cd ~/Murinae/Hydromyini/
mkdir Heterozygosity
cd Heterozygosity

# list contigs for each specie
mkdir contigs_lists
cd contigs_lists

cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt | parallel -j 20 ~/scripts/script_contig_list.sh {} Hydromyini

# calcul heterozygosity

# you need to install Angsd (put it in your bin)

cd ../

cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt | parallel -j 20 ~/scripts/script_angsd.sh {} Hydromyini 2>error_het.txt >out_het.txt


# make table that contains heterozygosity and mean coverage depth
cd ~/Murinae/Hydromyini/
mkdir Exons_infos
cd Exons_infos

# you need to install Samtools

cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt | parallel -j 20 ~/scripts/script_table_info_exons_V2.sh {} Hydromyini 2>error.txt >out.txt

for SPECIE in $(cat ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt) ; do
LC_ALL=en_US awk '{if ($3!=0) {print $0}}' ${SPECIE}_Exons_tab_infos.csv > ${SPECIE}_Exons_tab_infos_final.csv ;
done

#######
## 3 ## Filter alignments based on heterozygosity and mean coverage depth ###
#######

cd ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/

python3 ~/scripts/Filtre_exons.py ~/Murinae/Hydromyini/Data/list_species_Hydromyini.txt ~/Murinae/Hydromyini/Aligned_Sequences/alignments/list_exons_align.txt ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/ ~/Murinae/Hydromyini/Exons_infos/ 4 5 3 >out_filter_exons.txt 2>error_filter_exons.txt 

mkdir alignments_wt_ref_filtered
for file in *_filtered.fasta ; do mv $file alignments_wt_ref_filtered ; done

cd alignments_wt_ref_filtered
mkdir not_used

for align in *_filtered.fasta ; do 
if (($(grep -c ">" $align) < 3)) ; then mv $align not_used ; fi 
done

find -maxdepth 1 -type f -name "*.fasta"  | cut -f2 -d"_"  > Exons_filtered_final_list.txt 

#######
## 4 ## Remove sequences with stop codons ### 
#######

mkdir alignments_wt_ref_filtered_wt_stop
cd alignments_wt_ref_filtered_wt_stop

cat ~/Murinae/Hydromyini/Aligned_Sequences/alignments/alignments_wt_ref/alignments_wt_ref_filtered/Exons_filtered_final_list.txt | parallel -j20 ~/scripts/script_remove_seq_with_stop_codons_V2.sh {} Hydromyini >out_command_stop.txt 2>error_command_stop.txt

mkdir not_used

for align in *_wt_stop.fasta ; do 
if (($(grep -c ">" $align) < 3)) ; then mv $align not_used ; fi 
done

find -maxdepth 1 -type f -name "*.fasta"  | cut -f2 -d"_"  > Exons_filtered_final_list.txt 

