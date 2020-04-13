source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.6" # One Example of recombinants
width_png <- 1254
height_png <- 974
## TODO: make loop and function to create all the recombinants


# Get data ----------------------------------------------------------------

tab_frag <- read.table('./data/Recombinants_dataset_SNP_table.txt_BIS', header = TRUE)

# Order levels
tab_frag$Read_ID2 <- factor(tab_frag$Read_ID, levels = unique(tab_frag$Read_ID[order(tab_frag$DSB, tab_frag$middle_first_NA_zone)]))

selected_hotspot <- 'P9peak.chr11_10175985_10177504' # 'P9peak.chr15_93973100_93974147'
tab_frag$Genotype <- tab_frag$Zone

myTypes <- c(22,25,8)
names(myTypes) <- levels(tab_frag$Var_type)
myLabels <- c("DEL", "INS", "SNP")
names(myLabels) <- levels(tab_frag$Var_type)

# table coverage per base
table_cov <- read.table('./data/SUM.recal_reads.coverage_per_base_final', header=T)

# Table SNPs
table_SNPs <- read.table('./data/Variants_x_Hotspots.txt', header = T)

# Tabel recombination rates
table_recomb_rates <- read.table('./data/Recombination_rates.txt', header = T)

# Plot --------------------------------------------------------------------




# ggsave(paste0("/Users/maudgautier/Documents/These/R_projects/04_Analysis_Final_Data/Images/4_Recombinants_RED_YELLOW_CORRECTIONS_LAURENT/", selected_hotspot, ".png"))
# 
png(paste("./output/", fig_name, ".png", sep = ""), 
 width = width_png, height = height_png)
print(draw_hotspot_recombinants(tab_frag, table_cov, table_SNPs, 
                                table_recomb_rates, selected_hotspot))
dev.off()





