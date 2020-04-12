source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.1"
width_png <- 1400*0.75
height_png <- 1184*0.75



# Get data ----------------------------------------------------------------

tab_paigen <- read.table("./data/PositiveControls_with_all_peaks_compared_to_nb_events_detected_with_informative_fragments_Summed_nb_events_on_interval", header=T)
tab_paigen$CO_nb_corrected_for_paigen <- tab_paigen$Tot_recombinants/tab_paigen$LgInterval_kb



# Plot --------------------------------------------------------------------


png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

ggplot(tab_paigen) +
    geom_point(aes(x = log(1+RR_Male_cM_Mb), y = log(1+CO_nb_corrected_for_paigen)), size=2) +
    labs(x = expression(paste('CO rate measured by Paigen et al. (2008) (cM/Mb) (log-transformed)')),
         y = expression(paste('Number of events per kb (log-transformed)'))) +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=30)) +
    geom_abline(intercept = lm(log(1+tab_paigen$CO_nb_corrected_for_paigen) ~ log(1+tab_paigen$RR_Male_cM_Mb))$coeff[[1]], 
                slope = lm(log(1+tab_paigen$CO_nb_corrected_for_paigen) ~ log(1+tab_paigen$RR_Male_cM_Mb))$coeff[[2]], 
                linetype=3)

dev.off()

