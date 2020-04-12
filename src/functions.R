
formatter1000000 <- function(x){ 
    x/1000000 
}
formatter100 <- function(x){ 
    x/100 
}


std <- function(x) sd(x)/sqrt(length(x))
# sd_upper=function(x) { mean(x)+sd(x) }
# sd_lower=function(x) { mean(x)-sd(x) }
# std_upper=function(x) { mean(x)+std(x) }
# std_lower=function(x) { mean(x)-std(x) }


def_sizes <- function(x){ifelse(x == "SNP", 2,2)}


# Get right colours
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


# Rename categories of hotspots
label_categ <- function(x) {
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
status_labeller <- function(variable,value){
    status_names <- list(
        'CO'="Rec-1S",
        'Complex'="Rec-MS",
        'NCO'="Rec-2S"
    )
    return(status_names[value])
}



label_hot_categ_NOV_sep <- function(x) {
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


plot.multi.dens <- function(s, title = "")
{
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s)) {
        junk.x = c(junk.x, density(s[[i]])$x)
        junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = title)
    for(i in 1:length(s)) {
        lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
    }
}



# Automatise figure 8.3 ---------------------------------------------------


## Pairwise correlations for all events
get_pairwise_correl_events <- function(table_HFM1_bakgrounds_per_hotspot_with_chr_size, sample1, sample2) {
    
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
get_pairwise_correl_Rec1S <- function(sample1, sample2) {
    
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




