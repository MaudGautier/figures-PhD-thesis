source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.6" # One Example of recombinants
width_png <- 1254
height_png <- 974
## TODO: make loop and function to create all the recombinants


# Get data ----------------------------------------------------------------

### Dataset 1: Recombinants

# Get table
tab_frag <- read.table('./data/Recombinants_dataset_SNP_table.txt_BIS', header = TRUE)

# Sort reads
tab_frag$Read_ID2 <- factor(tab_frag$Read_ID, levels = unique(tab_frag$Read_ID[order(tab_frag$DSB, tab_frag$middle_first_NA_zone)]))
tab_frag$Genotype <- tab_frag$Zone


## INTULE ??
# myTypes <- c(22,25,8)
# names(myTypes) <- levels(tab_frag$Var_type)
# myLabels <- c("DEL", "INS", "SNP")
# names(myLabels) <- levels(tab_frag$Var_type)

### Dataset 2: Per-base coverage
table_cov <- read.table('./data/SUM.recal_reads.coverage_per_base_final', header=T)

### Dataset 3: Variants
table_SNPs <- read.table('./data/Variants_x_Hotspots.txt', header = T)

### Dataset 4: Swith intervals distribution
table_recomb_rates <- read.table('./data/Recombination_rates.txt', header = T)



# Plot --------------------------------------------------------------------

for (selected_hotspot in readLines("./data/list_hotspots.txt")) { 
    
    # selected_hotspot <- 'P9peak.chr11_10175985_10177504'

    # Open output file
    png(paste("./output/draw_recombinants/", selected_hotspot, ".png", sep = ""), 
        width = width_png, height = height_png)
    
    # Plot recombinants
    print(draw_hotspot_recombinants(tab_frag, table_cov, table_SNPs, 
                                table_recomb_rates, selected_hotspot))

    # Close output file
    dev.off()
}

