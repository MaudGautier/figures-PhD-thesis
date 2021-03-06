mkmount

# Figure 5.2
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/03_General_Statistics/04_Statistics_on_recombinants/Recombinants_per_hotspot.txt_BIS ./data/Recombinants_per_hotspot.txt

# Figure 5.3
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/03_General_Statistics/02_Bias_checks/3_Capture/ALL_columns_with_proportions_and_DMC1_and_dBGC.txt ./data/ALL_columns_with_proportions_and_DMC1_and_dBGC.txt

# Figure 6.1
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/08_Compare_with_Paigen//02_Results/PositiveControls_with_all_peaks_compared_to_nb_events_detected_with_informative_fragments_Summed_nb_events_on_interval.txt ./data/PositiveControls_with_all_peaks_compared_to_nb_events_detected_with_informative_fragments_Summed_nb_events_on_interval.txt

# Figure 6.2 and 6.3
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/03_General_Statistics/05_Distribution_of_recombinants/Intensity_and_counts_good_classes_and_DMC1_with_nb_frags_per_hotspot.txt ./data/Intensity_and_counts_good_classes_and_DMC1_with_nb_frags_per_hotspot.txt

# Figure 6.5
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/04_dBGC_Analysis/ALL.Fragments.givers_with_DMC1_SUM_59_60.txt ./data/ALL.Fragments.givers_with_DMC1_SUM_59_60.txt

# Figure 6.6 (and all hotspots with recombinants)
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/06_Drawings/Recombinants_dataset_SNP_table.txt_BIS ./data/Recombinants_dataset_SNP_table.txt_BIS
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/06_Drawings/Per_base_coverages/SUM.recal_reads.coverage_per_base_final ./data/SUM.recal_reads.coverage_per_base_final
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/03_General_Statistics/00_Hotspot_information/1_Polymorphism/Variants_x_Hotspots.txt ./data/Variants_x_Hotspots.txt
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/03_General_Statistics/00_Hotspot_information/2_Peak_detectability/Recombination_rates.txt ./data/Recombination_rates.txt 

# scp gautier@pbil-gates.univ-lyon1.fr:/beegfs/data/gautier/0_Genomes/Mus_musculus/domesticus/02_Baits/baits_mm10.bed ./data/
cut -f4 ./data/baits_mm10.bed > ./data/list_hotspots.txt

# Figure 6.7
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/02_Recombinants_and_False_Positives_Dataset/Recombination_events_per_hotspot_with_informative_fragments.txt ./data/Recombination_events_per_hotspot_with_informative_fragments.txt

# Figures 7.2
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/1_Hotspots/2_P9_target_motif_analysis_DEFINITIVE_VERSION/5_Processed_findings/P9_Dom.all_selected_hotspots.Best_score.joined_haplotypes.txt ./data/P9_Dom.all_selected_hotspots.Best_score.joined_haplotypes.txt
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/1_Hotspots/2_P9_target_motif_analysis_DEFINITIVE_VERSION/5_Processed_findings/P9_Cast.all_selected_hotspots.Best_score.joined_haplotypes.txt ./data/P9_Cast.all_selected_hotspots.Best_score.joined_haplotypes.txt

# Figure 7.3
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/04_dBGC_Analysis/ALL.Fragments.givers.txt ./data/ALL.Fragments.givers.txt

# Figure 7.4
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/6_ThirdSequencing/05_Analyses_of_Recombinants/03_General_Statistics/00_Hotspot_information/1_Polymorphism/Variants_x_Hotspots_with_categories_and_polym.txt ./data/Variants_x_Hotspots_with_categories_and_polym.txt

# Figure 7.5
cp /Users/maudgautier/Documents/These/02_pandata_sshfs/2_dBGC/8_Simulations_ABC/Results/Simul_relationship_SNP_density_prop_NCO1_real_mutational_process/Summary/summary.txt ./data/NCO1_summary.txt

# Figure 8.2
cp /Users/maudgautier/Documents/These/_RANGE/HFM1_for_figures/List_genotypes_hotspots_refined_with_chr_size.txt ./data/HFM1_List_genotypes_hotspots_refined_with_chr_size.txt

# Figure 8.3
cp /Users/maudgautier/Documents/These/_RANGE/HFM1_for_figures/All_HFM1_events_target_names.txt ./data/HFM1_All_recombinants_target_names.txt
cp /Users/maudgautier/Documents/These/_RANGE/HFM1_for_figures/Tot_fragments.txt ./data/HFM1_all_fragments.txt
