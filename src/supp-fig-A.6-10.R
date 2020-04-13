source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "supp-figure-A.6"
width_png <- 1200*0.75
height_png <- 1000*0.75



# Get data ----------------------------------------------------------------

## Identical as for Figure 8.2
table_HFM1_bakgrounds_per_hotspot_with_chr_size <- read.table("./data/HFM1_List_genotypes_hotspots_refined_with_chr_size.txt", header=T, comment="")

# Reorder chromosomes
chrOrderbis<-c(paste("chr",1:19,sep=""),"chrX","chrY","chrM")
table_HFM1_bakgrounds_per_hotspot_with_chr_size$Chr <- factor(table_HFM1_bakgrounds_per_hotspot_with_chr_size$Chr, levels=chrOrderbis)

# Replace NOCOV with NA
table_HFM1_bakgrounds_per_hotspot_with_chr_size[table_HFM1_bakgrounds_per_hotspot_with_chr_size=="NOCOV"] <- NA


## Specific to Figure 8.3 and A.6
# Recombinants in HFM1 analysis
tab_HFM1_events <- read.table("./data/HFM1_All_recombinants_target_names.txt", 
                              header=T, comment="")

# All fragments for HFM1 analysis
tab_HFM1_fragments_sequenced <- read.table("./data/HFM1_all_fragments.txt", 
                                           header=T)



# Plot --------------------------------------------------------------------


for (pair_samples in list(c("S28353", "S28355"),
                          c("S28353", "S28367"),
                          c("S28353", "S28371"),
                          c("S28367", "S28355"),
                          c("S28367", "S28371"),
                          c("S28355", "S28371"))) { 
    
    # Get sample names
    sample1 <- pair_samples[1]
    sample2 <- pair_samples[2]
    
    # Open output file
    png(paste("./output/", fig_name, "-", sample1, "-vs-", sample2, ".png", sep = ""), 
        width = width_png, height = height_png)
    
    # Plot corrrelation
    plot_correl <- get_pairwise_correl_Rec1S(table_HFM1_bakgrounds_per_hotspot_with_chr_size, 
                                             tab_HFM1_fragments_sequenced, 
                                             sample1, sample2)
    print(plot_correl)
    
    # Close output file
    dev.off()
}



