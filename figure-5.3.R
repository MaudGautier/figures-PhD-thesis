source("libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-5.3"
width_png <- 1677*0.75
height_png <- 1227*0.75



# Get data ----------------------------------------------------------------

tab_capture_eff <- read.table('./data/ALL_columns_with_proportions_and_DMC1_and_dBGC.txt', header=T)



# Plot --------------------------------------------------------------------

png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

ggplot(tab_capture_eff, aes(x=tot_B6_frags/(tot_B6_frags + tot_CAST_frags + tot_Recomb_frags), ..count../sum(..count..))) + 
    geom_histogram(fill = "grey", binwidth = 0.01) + 
    theme_classic() +
    coord_cartesian(xlim = c(0,1)) + 
    geom_vline(xintercept=quantile(tab_capture_eff$tot_B6_frags/(tab_capture_eff$tot_B6_frags+tab_capture_eff$tot_CAST_frags+ tab_capture_eff$tot_Recomb_frags), c(0.025), na.rm=T)[[1]], linetype=3) + 
    geom_vline(xintercept=quantile(tab_capture_eff$tot_B6_frags/(tab_capture_eff$tot_B6_frags+tab_capture_eff$tot_CAST_frags+ tab_capture_eff$tot_Recomb_frags), c(0.975), na.rm=T)[[1]], linetype=3) + 
    geom_vline(xintercept=quantile(tab_capture_eff$tot_B6_frags/(tab_capture_eff$tot_B6_frags+tab_capture_eff$tot_CAST_frags+ tab_capture_eff$tot_Recomb_frags), c(0.5), na.rm=T)[[1]], linetype=2) + 
    labs(y = "Proportion of targets", x = "Proportion of non-recombinant fragments coming from the B6 haplotype") +
    scale_y_continuous(labels= scales::percent) +
    theme(text=element_text(family="LM Roman 10", size = 35))

dev.off()

