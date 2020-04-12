source("src/libraries.R")



######## Au clair les graphes proportions de GC Sommés ########

# Get data
tab_GC_profile <- read.table('./data/Variants_x_Hotspots_with_categories_and_polym.txt', header=T)
# Ajout variables
#tab_GC_profile$rel_pos <- tab_GC_profile$POS-tab_GC_profile$DSB
tab_GC_profile$fact_rel_pos <- as.factor(tab_GC_profile$POS-tab_GC_profile$DSB)
tab_GC_profile$fact_rel_pos_LEVELS_ORDERED <- factor(tab_GC_profile$fact_rel_pos, levels = as.factor(c(seq(-1500,1500))))


# Category
list_categ_CAST <- c("tC.hB", "NOV.tC.hB", "tC.hB+C") ; title_categ <- expression(paste("PRDM9"^"CAST", "-targeted hotspots")) 
# Un poil mieux avec les symétriques inclus

subtab_CAST <- tab_GC_profile[which(tab_GC_profile$CATEGORY %in% list_categ_CAST & tab_GC_profile$POLYM_TYPE =="WS"),]
(nb_hotspots <- length(unique(subtab_CAST$PEAK_NAME)))


vec_GC_counts_B6_for_CAST <- (table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"G"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"C"])
vec_GC_counts_CAST_for_CAST <- (table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"G"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"C"])

# Crée vecteurs de comptage des A C T G par base
vec_GC_props_B6_for_CAST <- (table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"G"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"C"])/(table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"G"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"C"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"A"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$B6))[,"T"])
vec_GC_props_CAST_for_CAST <- (table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"G"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"C"])/((table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"G"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"C"])+(table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"A"]+table(subtab_CAST$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_CAST$CAST))[,"T"]))


wid_choice <- 300
require(zoo)
data_plot_no_melt_CAST <- data.frame(Pos = seq(as.numeric(names(rollapply(vec_GC_counts_B6_for_CAST, width = wid_choice, by=1, FUN=mean, align="center")[1])),as.numeric(names(tail(rollapply(vec_GC_counts_B6_for_CAST, width = wid_choice, by=1, FUN=mean, align="center"),1)))), 
                                     mean_prop_GC_B6 = rollapply(vec_GC_props_B6_for_CAST, width = wid_choice, by=1, FUN=mean, align="center", na.rm=T),# / nb_hotspots, 
                                     mean_prop_GC_CAST = rollapply(vec_GC_props_CAST_for_CAST, width = wid_choice, by=1, FUN=mean, align="center", na.rm=T))# / nb_hotspots)
require(reshape2)
data_plot_CAST <- melt(data_plot_no_melt_CAST, id = "Pos") 
data_plot_CAST$variable <- factor(data_plot_CAST$variable, levels = c("mean_prop_GC_B6", "mean_prop_GC_CAST"), labels = c("B6", "CAST"))


# # subset CAST
# ggplot(data = droplevels(subset(data_plot_CAST,data_plot_CAST$variable %in% c("CAST")))) + 
#     geom_line(aes(x = Pos, y = value, colour = variable)) +
#     ylim(c(0,1)) +
#     geom_hline(aes(yintercept=0.5), linetype=2, color ="lightgrey") +
#     xlim(c(-1500,1500)) +
#     theme_classic() +
#     theme(text=element_text(family="LM Roman 10", size = 35), 
#           plot.title = element_text(hjust = 0.5)) +
#     labs(x="Position relative to the inferred DSB site (bp)", y="Proportion of S alleles in WS polymorphic sites\noriginating from the CAST lineage") +
#     scale_colour_manual("Haplotype", values=c("#FFCC00"), labels = c("CAST")) + 
#     theme(legend.position="none") +
#     geom_vline(aes(xintercept=0), linetype=2, color ="lightgrey") +
#     geom_vline(aes(xintercept=-250), linetype=3, color ="lightgrey") +
#     geom_vline(aes(xintercept=250), linetype=3, color ="lightgrey")
# 







list_categ_B6 <- c("tB.hC", "NOV.tB.hC", "tB.hB+C") ; title_categ <- expression(paste("PRDM9"^"B6", "-targeted hotspots"))


subtab_B6 <- tab_GC_profile[which(tab_GC_profile$CATEGORY %in% list_categ_B6 & tab_GC_profile$POLYM_TYPE =="WS"),]
(nb_hotspots <- length(unique(subtab_B6$PEAK_NAME)))


vec_GC_counts_B6_for_B6 <- (table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"G"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"C"])
vec_GC_counts_CAST_for_B6 <- (table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"G"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"C"])

# Crée vecteurs de comptage des A C T G par base
vec_GC_props_B6_for_B6 <- (table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"G"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"C"])/(table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"G"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"C"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"A"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$B6))[,"T"])
vec_GC_props_CAST_for_B6 <- (table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"G"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"C"])/((table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"G"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"C"])+(table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"A"]+table(subtab_B6$fact_rel_pos_LEVELS_ORDERED, droplevels(subtab_B6$CAST))[,"T"]))


wid_choice <- 300
require(zoo)
data_plot_no_melt_B6 <- data.frame(Pos = seq(as.numeric(names(rollapply(vec_GC_counts_B6_for_B6, width = wid_choice, by=1, FUN=mean, align="center")[1])),as.numeric(names(tail(rollapply(vec_GC_counts_B6_for_B6, width = wid_choice, by=1, FUN=mean, align="center"),1)))), 
                                   mean_prop_GC_B6 = rollapply(vec_GC_props_B6_for_B6, width = wid_choice, by=1, FUN=mean, align="center", na.rm=T),# / nb_hotspots, 
                                   mean_prop_GC_CAST = rollapply(vec_GC_props_CAST_for_B6, width = wid_choice, by=1, FUN=mean, align="center", na.rm=T))# / nb_hotspots)
require(reshape2)
data_plot_B6 <- melt(data_plot_no_melt_B6, id = "Pos") 
data_plot_B6$variable <- factor(data_plot_B6$variable, levels = c("mean_prop_GC_B6", "mean_prop_GC_CAST"), labels = c("B6", "CAST"))
# 
# 
# 
# # subset B6
# ggplot(data = droplevels(subset(data_plot_B6,data_plot_B6$variable %in% c("B6")))) + 
#     geom_line(aes(x = Pos, y = value, colour = variable)) +
#     ylim(c(0,1)) +
#     geom_hline(aes(yintercept=0.5), linetype=2, color ="lightgrey") +
#     xlim(c(-1500,1500)) +
#     theme_classic() +
#     theme(text=element_text(family="LM Roman 10", size = 35), 
#           plot.title = element_text(hjust = 0.5)) +
#     labs(x="Position relative to the inferred DSB site (bp)", y="Proportion of S alleles in WS polymorphic sites\noriginating from the B6 lineage") +
#     scale_colour_manual("Haplotype", values=c("#993333"), labels = c("B6")) + 
#     theme(legend.position="none") +
#     geom_vline(aes(xintercept=0), linetype=2, color ="lightgrey") +
#     geom_vline(aes(xintercept=-250), linetype=3, color ="lightgrey") +
#     geom_vline(aes(xintercept=250), linetype=3, color ="lightgrey")
# 




# Combine the two

ggplot() + 
    geom_line(data = droplevels(subset(data_plot_CAST,data_plot_CAST$variable %in% c("CAST"))), aes(x = Pos, y = value, colour = variable), size=1.1) +
    geom_line(data = droplevels(subset(data_plot_B6,data_plot_B6$variable %in% c("B6"))), aes(x = Pos, y = value, colour = variable), size=1) +
    ylim(c(0,1)) +
    geom_hline(aes(yintercept=0.5), linetype=2, color ="lightgrey") +
    xlim(c(-1500,1500)) +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size = 35), 
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(margin = margin(l = 10, r = 15))) +
    labs(x="Position relative to the PRDM9 ChIP-seq peak summit (bp)", y="Proportion of S alleles originating from the self lineage") +
    #scale_colour_manual("Haplotype", values=c("#993333", "#FFCC00"), labels = c("B6", "CAST")) + 
    scale_colour_manual("Haplotype", values=c("#B95353", "#EFBC00"), labels = c("B6", "CAST")) + 
    theme(legend.position="none") +
    geom_vline(aes(xintercept=0), linetype=2, color ="lightgrey") +
    geom_vline(aes(xintercept=-250), linetype=3, color ="lightgrey") +
    geom_vline(aes(xintercept=250), linetype=3, color ="lightgrey") 


