source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-7.3"
width_png <- 1200*0.75
height_png <- 1100*0.75



# Get data ----------------------------------------------------------------

tab_givers <- read.table('./data/ALL.Fragments.givers.txt', header=T)

tab_givers$tot_B6 <- tab_givers$Nb_B_giver_in_complex + tab_givers$Nb_B_giver_in_simple
tab_givers$tot_CAST <- tab_givers$Nb_C_giver_in_complex + tab_givers$Nb_C_giver_in_simple
tab_givers$tot <- tab_givers$tot_CAST + tab_givers$tot_B6
tab_givers$prop_CAST <- tab_givers$tot_CAST/tab_givers$tot
tab_givers$Hot_categ_bis_nov_sep <- sapply(tab_givers$Hotspot_category, relabel_hotspots_by_P9_target)
tab_givers$dBGC_coeff <- 2*tab_givers$prop_CAST-1

alpha_val=0.5


# Plot --------------------------------------------------------------------

# Open output file
png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

# Plot
ggplot() + 
    geom_density(aes(tab_givers[which(tab_givers$Hot_categ_bis_nov_sep=="Novel_B6"),]$dBGC_coeff, y=(..count..)/sum(..count..), fill="darkgrey", colour= "darkgrey"), alpha=alpha_val) +
    geom_density(aes(tab_givers[which(tab_givers$Hot_categ_bis_nov_sep=="Novel_CAST"),]$dBGC_coeff, y=(..count..)/sum(..count..), fill="lightgrey", colour= "lightgrey"), alpha=alpha_val) +
    geom_density(aes(tab_givers[which(tab_givers$Hot_categ_bis_nov_sep=="B6"),]$dBGC_coeff, y=(..count..)/sum(..count..), fill="#993333", colour= "#993333"), alpha=alpha_val) +
    geom_density(aes(tab_givers[which(tab_givers$Hot_categ_bis_nov_sep=="CAST"),]$dBGC_coeff, y=(..count..)/sum(..count..), fill="#FFCC00", colour= "#FFCC00"), alpha=alpha_val) +
    scale_y_continuous(labels= scales::percent) +
    scale_fill_manual('PRDM9 target',
                      values = c("#993333", "#FFCC00", "darkgrey", "lightgrey"),
                      labels = c(bquote(PRDM9^Dom2), bquote(PRDM9^Cst), bquote(PRDM9^Dom2~(novel)), bquote(PRDM9^Cst~(novel))),
                      #labels = c(bquote(PRDM9^Dom2), bquote(PRDM9^Cst), expression(paste(PRDM9^Dom2, ' (novel)')), expression(paste(PRDM9^Cst, ' (novel)'))),
                      guide = guide_legend(override.aes = list(alpha = alpha_val))) +
    scale_colour_manual('PRDM9 target',
                        values = c("#993333", "#FFCC00", "darkgrey", "lightgrey"),
                        labels = c(bquote(PRDM9^Dom2), bquote(PRDM9^Cst), bquote(PRDM9^Dom2~(novel)), bquote(PRDM9^Cst~(novel))),
                        #labels = c(bquote(PRDM9^Dom2), bquote(PRDM9^Cst), expression(paste(PRDM9^Dom2, ' (novel)')),expression(paste(PRDM9^Cst, ' (novel)'))),
                        guide = guide_legend(override.aes = list(alpha = alpha_val))) +
    theme_classic() +
    labs(y = "Proportion of hotspots of the category (%)",
         x = "dBGC coefficient") +
    theme(text=element_text(family="LM Roman 10", size=35), 
          legend.position=c(0.5,0.88),
          legend.key.size = unit(2, 'lines'))

# Close output file
dev.off()

