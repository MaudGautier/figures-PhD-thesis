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






########## REDO DRAW RECOMBINNNT WITH RED-YELLOW COLORS AND COMMENTS FROM LAURENT #####


categ <- label_categ(tab_frag[which(tab_frag$Hotspot_ID == selected_hotspot),]$Hotspot_category[1])
nb_recomb <- length(unique(tab_frag[which(tab_frag$Hotspot_ID == selected_hotspot),]$Read_ID))

plot1 <- ggplot(tab_frag[which(tab_frag$Hotspot_ID == selected_hotspot & tab_frag$Status != "Complex"),], aes(colour=Genotype)) + 
    theme_classic() +
    geom_segment(mapping=aes(x=Start-DSB, y=Read_ID2, xend=Stop-DSB, yend=Read_ID2, size=0.7)) +
    geom_point(aes(x=SNP_pos-DSB, y = Read_ID2, colour=SNP_genot, shape=SNP_genot, size=2)) +
    #geom_segment(mapping=aes(x=Tract_start-DSB, y=Read_ID_bis, xend=Tract_end-DSB, yend=Read_ID_bis, size=4), colour = "lightgrey") +
    scale_size_identity() +
    #coord_cartesian(xlim = c(-1000, 1000)) +
    #scale_x_continuous(limits=c(-1000, 1000)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position=c(0.94,0.5)) + 
    labs(title = paste0(selected_hotspot, "     -     ", categ), x = "Coordinates (relative to the hotspot center)", y = '') +
    facet_grid(Status ~ ., scales = 'free', switch = "y", labeller=status_labeller) +
    geom_vline(aes(xintercept=0), linetype=2, color ="darkgrey") +
    geom_vline(aes(xintercept=500), linetype=3, color ="darkgrey") +
    geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
    scale_shape_manual('Marker', values = c(16,16), labels = c("B6", "CAST"), drop = FALSE, na.translate = FALSE, 
                       guide = guide_legend(override.aes = list(
                           linetype = c(rep("blank", 2)),
                           shape = c(rep(16, 2)),
                           colour = c("#B95353", "#EFBC00")))) +
    scale_colour_manual('Haplotype', values=c("#B95353", "#EFBC00"), labels = c("B6", "CAST"), 
                        na.value = "grey50", 
                        guide = guide_legend(override.aes = list(
                            linetype = c(rep("solid", 2), "solid"),
                            shape = c(rep(NA, 2), NA)))) +
    theme(text=element_text(size = 15), 
          plot.title = element_text(hjust = 0.5))


sub_table_cov <- table_cov[which(table_cov$Target_ID == selected_hotspot), ]
plot2 <- ggplot(sub_table_cov, aes(x = Position - DSB, y = Coverage)) +
    geom_smooth(method = "loess", se = FALSE, span = 1/10, color = "darkgrey") + 
    #coord_cartesian(xlim = c(-1000, 1000)) +
    scale_x_continuous(limits=c(-1000, 1000)) +
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_classic() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_vline(aes(xintercept=0), linetype=2, color ="darkgrey") +
    geom_vline(aes(xintercept=500), linetype=3, color ="darkgrey") +
    geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
    labs(x = '', y = "Coverage") +
    theme(text=element_text(size = 13))

plot3 <- ggplot(table_SNPs[which(table_SNPs$PEAK_NAME == selected_hotspot),]) + 
    geom_point(aes(x=POS-DSB, y = PEAK_NAME, shape = SNP_TYPE, size = 2, colour="val1")) + 
    scale_size_identity() +
    theme_classic() +
    scale_x_continuous(limits=c(-1000, 1000)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.title.x=element_blank(), legend.position="none") + 
    labs(x = '', y = "Markers") + 
    scale_shape_manual(values = c(16,16,16), labels = myLabels, drop = FALSE) +
    geom_vline(aes(xintercept=0), linetype=2, color ="darkgrey") +
    geom_vline(aes(xintercept=500), linetype=3, color ="darkgrey") +
    geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
    scale_color_manual(values = c("val1" = "darkgrey")) +
    theme(text=element_text(size = 13))

plot5 <- ggplot(table_recomb_rates[which(table_recomb_rates$Peak_name == selected_hotspot),]) +
    geom_rect(mapping=aes(xmin = Start - DSB, xmax = Stop - DSB, ymin = 0, ymax=Potential_recombinants, fill = "grey68"), color = "grey51") +
    theme_classic() + 
    geom_rect(mapping=aes(xmin = Start -DSB, xmax = Stop -DSB, ymin = 0, ymax=Potential_COs, fill = "lightgrey"), color = "grey51") +
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
    geom_vline(aes(xintercept=0), linetype=2, color ="darkgrey") +
    geom_vline(aes(xintercept=500), linetype=3, color ="darkgrey") +
    geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
    scale_x_continuous(limits=c(-1000, 1000)) +
    labs(x = '', y = "# Informative\nfragments") +
    scale_fill_manual('',
                      values = c('grey68', 'lightgrey'), 
                      labels = c('Rec-2S', 'Rec-1S'),
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.position=c(0.94,0.5)) +
    theme(text=element_text(size = 13))


tracks(plot5, plot2, plot3, plot1,
       heights=c(1.5,1.1,1.1,10), 
       main = paste0(selected_hotspot, " (", categ, "): ", nb_recomb, " recombination events"), 
       xlab = "Position relative to the inferred DSB site (bp)",
       xlim=c(-1000,1000))

# ggsave(paste0("/Users/maudgautier/Documents/These/R_projects/04_Analysis_Final_Data/Images/4_Recombinants_RED_YELLOW_CORRECTIONS_LAURENT/", selected_hotspot, ".png"))
# 
# png(paste("./output/", fig_name, ".png", sep = ""), 
#     width = width_png, height = height_png)
# 
# 
# dev.off()





