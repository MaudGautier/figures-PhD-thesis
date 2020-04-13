# Functions used for several scripts --------------------------------------

# Data formatting
formatter1000000 <- function(x) { x/1000000 }
formatter100 <- function(x) { x/100 }

# Standard error function 
# (standard deviation divided by the square root of the sample size)
std <- function(x) sd(x)/sqrt(length(x))

# Rename categories of hotspots
relabel_hotspot_category <- function(x) {
    if (x=="tC.hB+C") {return("tC.sym")}
    if (x=="tC.hB")   {return("tC.chC")}
    if (x=="tC.hC")   {return("tC.chB")}
    if (x=="tB.hB")   {return("tB.chC")}
    if (x=="tB.hB+C") {return("tB.sym")}
    if (x=="tB.hC")   {return("tB.chB")}
    if (x=="NOV.tB.hC") {return("NOV.tB.chB")}
    if (x=="NOV.tC.hB") {return("NOV.tC.chC")}
    if (x=="NOV.NA")  {return("NOV.sym")}
    if (x=="INDEP.hB") {return("INDEP.chC")}
    if (x=="INDEP.hB+C") {return("INDEP.sym")}
    if (x=="INDEP.hC") {return("INDEP.chB")}
}

# Regroup hotspots by PRDM9 target
relabel_hotspots_by_P9_target <- function(x) {
    if (is.na(x)) {return(NA)}
    if (x=="tC.hB+C") {return("CAST")}
    if (x=="tC.hB")   {return("CAST")}
    if (x=="tC.hC")   {return("CAST")}
    if (x=="tB.hB")   {return("B6")}
    if (x=="tB.hB+C") {return("B6")}
    if (x=="tB.hC")   {return("B6")}
    if (x=="NOV.tB.hC") {return("Novel_B6")}
    if (x=="NOV.tC.hB") {return("Novel_CAST")}
    if (x=="NOV.NA")  {return(NA)}
    if (x=="INDEP.hB") {return(NA)}
    if (x=="INDEP.hB+C") {return(NA)}
    if (x=="INDEP.hC") {return(NA)}
}

# Relabel types of recombination events
rename_events <- function(variable,value){
    status_names <- list(
        'CO'="Rec-1S",
        'Complex'="Rec-MS",
        'NCO'="Rec-2S"
    )
    return(status_names[value])
}



# Figure 6.4 --------------------------------------------------------------

def_sizes <- function(x) { ifelse(x == "SNP", 2, 2) }



# Figure 6.6 --------------------------------------------------------------

# Draw recombinants of hotspot
draw_hotspot_recombinants <- function(tab_frag, table_cov, table_SNPs, 
                                      table_recomb_rates, selected_hotspot) {
    
    # Info for the title (Hotspot category and number of recombinants)
    categ <- relabel_hotspot_category(tab_frag[which(tab_frag$Hotspot_ID == selected_hotspot),]$Hotspot_category[1])
    nb_recomb <- length(unique(tab_frag[which(tab_frag$Hotspot_ID == selected_hotspot),]$Read_ID))
    
    # Plot recombinants
    plot1 <- ggplot(tab_frag[which(tab_frag$Hotspot_ID == selected_hotspot & tab_frag$Status != "Complex"),], aes(colour=Genotype)) + 
        theme_classic() +
        geom_segment(mapping = aes(x = Start - DSB, y = Read_ID2, 
                                   xend = Stop - DSB, yend = Read_ID2, 
                                   size = 0.7)) +
        geom_point(aes(x = SNP_pos - DSB, y = Read_ID2, 
                       colour = SNP_genot, shape = SNP_genot, size = 2)) +
        scale_size_identity() +
        theme(axis.text.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              legend.position=c(0.94,0.5)) + 
        labs(title = paste0(selected_hotspot, "     -     ", categ), 
             x = "Coordinates (relative to the hotspot center)", y = '') +
        facet_grid(Status ~ ., 
                   scales = 'free', switch = "y", 
                   labeller = rename_events) +
        geom_vline(aes(xintercept = 0),    linetype=2, color ="darkgrey") +
        geom_vline(aes(xintercept = 500),  linetype=3, color ="darkgrey") +
        geom_vline(aes(xintercept = -500), linetype=3, color ="darkgrey") +
        scale_shape_manual('Marker', 
                           values = c(16,16), 
                           labels = c("B6", "CAST"), drop = FALSE, 
                           na.translate = FALSE, 
                           guide = guide_legend(override.aes = list(
                               linetype = c(rep("blank", 2)),
                               shape = c(rep(16, 2)),
                               colour = c(color_B6_2, color_CA_2)))) +
        scale_colour_manual('Haplotype', 
                            values=c(color_B6_2, color_CA_2), 
                            labels = c("B6", "CAST"), 
                            na.value = "grey50", 
                            guide = guide_legend(override.aes = list(
                                linetype = c(rep("solid", 2), "solid"),
                                shape = c(rep(NA, 2), NA)))) +
        theme(text=element_text(size = 15), 
              plot.title = element_text(hjust = 0.5))
    
    # Plot coverage
    sub_table_cov <- table_cov[which(table_cov$Target_ID == selected_hotspot), ]
    plot2 <- ggplot(sub_table_cov, aes(x = Position - DSB, y = Coverage)) +
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = "darkgrey") + 
        scale_x_continuous(limits = c(-1000, 1000)) +
        scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
        theme_classic() + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        geom_vline(aes(xintercept = 0),    linetype=2, color ="darkgrey") +
        geom_vline(aes(xintercept = 500),  linetype=3, color ="darkgrey") +
        geom_vline(aes(xintercept = -500), linetype=3, color ="darkgrey") +
        labs(x = '', y = "Coverage") +
        theme(text = element_text(size = 13))
    
    # Plot variants
    plot3 <- ggplot(table_SNPs[which(table_SNPs$PEAK_NAME == selected_hotspot),]) + 
        geom_point(aes(x = POS-DSB, y = PEAK_NAME, shape = SNP_TYPE, size = 2, colour="val1")) + 
        scale_size_identity() +
        theme_classic() +
        scale_x_continuous(limits=c(-1000, 1000)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), 
              axis.title.x = element_blank(), legend.position = "none") + 
        labs(x = '', y = "Markers") + 
        scale_shape_manual(values = c(16,16,16), labels = myLabels, drop = FALSE) +
        geom_vline(aes(xintercept=0),    linetype=2, color ="darkgrey") +
        geom_vline(aes(xintercept=500),  linetype=3, color ="darkgrey") +
        geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
        scale_color_manual(values = c("val1" = "darkgrey")) +
        theme(text=element_text(size = 13))
    
    # Plot interval distribution
    plot5 <- ggplot(table_recomb_rates[which(table_recomb_rates$Peak_name == selected_hotspot),]) +
        geom_rect(mapping = aes(xmin = Start - DSB, xmax = Stop - DSB, ymin = 0, ymax = Potential_recombinants, fill = "grey68"), color = "grey51") +
        theme_classic() + 
        geom_rect(mapping=aes(xmin = Start -DSB, xmax = Stop -DSB, ymin = 0, ymax=Potential_COs, fill = "lightgrey"), color = "grey51") +
        scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
        geom_vline(aes(xintercept=0),    linetype=2, color ="darkgrey") +
        geom_vline(aes(xintercept=500),  linetype=3, color ="darkgrey") +
        geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
        scale_x_continuous(limits=c(-1000, 1000)) +
        labs(x = '', y = "# Informative\nfragments") +
        scale_fill_manual('',
                          values = c('grey68', 'lightgrey'), 
                          labels = c('Rec-2S', 'Rec-1S'),
                          guide = guide_legend(override.aes = list(alpha = 1))) +
        theme(legend.position=c(0.94,0.5)) +
        theme(text=element_text(size = 13))
    
    # Combine plots
    combination <- ggbio::tracks(plot5, plot2, plot3, plot1,
                                 heights=c(1.5,1.1,1.1,10), 
                                 main = paste0(selected_hotspot, " (", categ, "): ", nb_recomb, " recombination events"), 
                                 xlab = "Position relative to the inferred DSB site (bp)",
                                 xlim=c(-1000,1000))
    
    return(combination)
}



# Figure 8.3 --------------------------------------------------------------

## Pairwise correlations for all events
get_pairwise_correl_events <- function(table_HFM1_bakgrounds_per_hotspot_with_chr_size, 
                                       tab_HFM1_fragments_sequenced,
                                       tab_HFM1_events,
                                       sample1, sample2) {
    
    table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak <- table_HFM1_bakgrounds_per_hotspot_with_chr_size[grep("P9peak",table_HFM1_bakgrounds_per_hotspot_with_chr_size$X.Hotspot_ID),]
    
    # 1.2. Get col number
    col1 <- which(colnames(table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak)==sample1)
    col2 <- which(colnames(table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak)==sample2)
    
    # 1.3. Get list shared
    list_shared <- table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak[which(table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak[col1]=="HET" & table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak[col2]=="HET"),]$X.Hotspot_ID   
    
    
    ### 2. Subset hotspots shared in table
    subset_events_of_samples_in_hotspots_shared <- tab_HFM1_events[which(tab_HFM1_events$TARGET_NO_EDGES %in% list_shared & tab_HFM1_events$INDIVIDUAL_NAME %in% c(sample1, sample2)),]
    
    
    ### 3. Melt table into nb of events per peak name per sample
    
    # 3.1. RECUPERE TOUS LES HOTSPOTS AVEC AU MOINS UN DIFFERENT DE 0
    count_events_of_samples_in_hotspots_shared <- dcast(subset_events_of_samples_in_hotspots_shared, TARGET_NO_EDGES ~ INDIVIDUAL_NAME, fun.aggregate = length, value.var="INDIVIDUAL_NAME")
    #NB: penser a le faire en ne prenant que les CO
    
    # 3.2. Ajouter la liste de tous les hotspots shared et du nombre de reads sur chaque hotspot pour ces individus
    TARGET_NO_EDGES <- list_shared
    count_events_of_samples_in_hotspots_shared_ALL <- merge(count_events_of_samples_in_hotspots_shared, as.data.frame(TARGET_NO_EDGES), sort=T, by="TARGET_NO_EDGES", all=T)
    count_events_of_samples_in_hotspots_shared_ALL[is.na(count_events_of_samples_in_hotspots_shared_ALL)] <- 0
    
    # 3.3. Penser a ajouter apres le nombre de reads par hotspot
    count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags <- merge(count_events_of_samples_in_hotspots_shared,tab_HFM1_fragments_sequenced)
    lab_frags_sample1 <- paste0("Tot_nb_fragments_", sample1)
    lab_div_sample1 <- paste0("Ratio_", sample1)
    count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_div_sample1]] <- count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[sample1]]/count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_frags_sample1]]
    
    lab_frags_sample2 <- paste0("Tot_nb_fragments_", sample2)
    lab_div_sample2 <- paste0("Ratio_", sample2)
    count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_div_sample2]] <- count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[sample2]]/count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_frags_sample2]]
    
    ## 4. Print correlation
    plot_correl_events <- ggplot(count_events_of_samples_in_hotspots_shared_ALL) +
        geom_point(aes_string(x = sample1, y = sample2), size=0.8) +
        labs(x = paste("Number of recombination events in mouse", sample1),
             y = paste("Number of recombination events in mouse", sample2))+
        theme_classic() +
        geom_abline(intercept = lm(count_events_of_samples_in_hotspots_shared_ALL[[sample2]] ~ count_events_of_samples_in_hotspots_shared_ALL[[sample1]])$coeff[[1]], 
                    slope = lm(count_events_of_samples_in_hotspots_shared_ALL[[sample2]] ~ count_events_of_samples_in_hotspots_shared_ALL[[sample1]])$coeff[[2]], 
                    linetype=3) +
        geom_abline(intercept = 0, slope=1, linetype = 5) +
        xlim(0,max(count_events_of_samples_in_hotspots_shared_ALL[[sample1]],count_events_of_samples_in_hotspots_shared_ALL[[sample2]])) +
        ylim(0,max(count_events_of_samples_in_hotspots_shared_ALL[[sample1]],count_events_of_samples_in_hotspots_shared_ALL[[sample2]])) +
        theme(text=element_text(family="LM Roman 10", size=30))

    return(plot_correl_events)
    
    
}

## Pairwise correlations for Rec-1S events only
get_pairwise_correl_Rec1S <- function(table_HFM1_bakgrounds_per_hotspot_with_chr_size, 
                                      tab_HFM1_fragments_sequenced,
                                      tab_HFM1_events,
                                      sample1, sample2) {
    
    table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak <- table_HFM1_bakgrounds_per_hotspot_with_chr_size[grep("P9peak",table_HFM1_bakgrounds_per_hotspot_with_chr_size$X.Hotspot_ID),]
    
    # 1.2. Get col number
    col1 <- which(colnames(table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak)==sample1)
    col2 <- which(colnames(table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak)==sample2)
    
    # 1.3. Get list shared
    list_shared <- table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak[which(table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak[col1]=="HET" & table_HFM1_bakgrounds_per_hotspot_with_chr_size_onlyP9peak[col2]=="HET"),]$X.Hotspot_ID   
    
    
    ### 2. Subset hotspots shared in table
    subset_events_of_samples_in_hotspots_shared <- tab_HFM1_events[which(tab_HFM1_events$TARGET_NO_EDGES %in% list_shared & tab_HFM1_events$INDIVIDUAL_NAME %in% c(sample1, sample2) & tab_HFM1_events$PRODUCT=="Rec-1S"),]
    
    
    ### 3. Melt table into nb of events per peak name per sample
    
    # 3.1. RECUPERE TOUS LES HOTSPOTS AVEC AU MOINS UN DIFFERENT DE 0
    count_events_of_samples_in_hotspots_shared <- dcast(subset_events_of_samples_in_hotspots_shared, TARGET_NO_EDGES ~ INDIVIDUAL_NAME, fun.aggregate = length, value.var="INDIVIDUAL_NAME")
    #NB: penser a le faire en ne prenant que les CO
    
    # 3.2. Ajouter la liste de tous les hotspots shared et du nombre de reads sur chaque hotspot pour ces individus
    TARGET_NO_EDGES <- list_shared
    count_events_of_samples_in_hotspots_shared_ALL <- merge(count_events_of_samples_in_hotspots_shared, as.data.frame(TARGET_NO_EDGES), sort=T, by="TARGET_NO_EDGES", all=T)
    count_events_of_samples_in_hotspots_shared_ALL[is.na(count_events_of_samples_in_hotspots_shared_ALL)] <- 0
    
    # 3.3. Penser a ajouter apres le nombre de reads par hotspot
    count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags <- merge(count_events_of_samples_in_hotspots_shared,tab_HFM1_fragments_sequenced)
    lab_frags_sample1 <- paste0("Tot_nb_fragments_", sample1)
    lab_div_sample1 <- paste0("Ratio_", sample1)
    count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_div_sample1]] <- count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[sample1]]/count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_frags_sample1]]
    
    lab_frags_sample2 <- paste0("Tot_nb_fragments_", sample2)
    lab_div_sample2 <- paste0("Ratio_", sample2)
    count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_div_sample2]] <- count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[sample2]]/count_events_of_samples_in_hotspots_shared_ALL_with_tot_frags[[lab_frags_sample2]]
    
    ## 4. Print correlation
    plot_correl_CO_events <- ggplot(count_events_of_samples_in_hotspots_shared_ALL) +
        geom_point(aes_string(x = sample1, y = sample2), size=0.8) +
        labs(x = paste("Number of Rec-1S events in mouse", sample1),
             y = paste("Number of Rec-1S events in mouse", sample2))+
        theme_classic() +
        geom_abline(intercept = lm(count_events_of_samples_in_hotspots_shared_ALL[[sample2]] ~ count_events_of_samples_in_hotspots_shared_ALL[[sample1]])$coeff[[1]], 
                    slope = lm(count_events_of_samples_in_hotspots_shared_ALL[[sample2]] ~ count_events_of_samples_in_hotspots_shared_ALL[[sample1]])$coeff[[2]], 
                    linetype=3) +
        geom_abline(intercept = 0, slope=1, linetype = 5) +
        xlim(0,max(count_events_of_samples_in_hotspots_shared_ALL[[sample1]],count_events_of_samples_in_hotspots_shared_ALL[[sample2]])) +
        ylim(0,max(count_events_of_samples_in_hotspots_shared_ALL[[sample1]],count_events_of_samples_in_hotspots_shared_ALL[[sample2]])) +
        theme(text=element_text(family="LM Roman 10", size=30))
    
    return(plot_correl_CO_events)
    
}


