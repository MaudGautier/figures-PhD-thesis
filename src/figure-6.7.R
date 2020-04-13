source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.7"
width_png <- 1550*0.75
height_png <- 1227*0.75



# Get data ----------------------------------------------------------------

tab_correl_extrapolation_CO_rate <- read.table("./data/Recombination_events_per_hotspot_with_informative_fragments.txt", header=T, comment="")

tab_correl_extrapolation_CO_rate$CO_rate_ABC <- tab_correl_extrapolation_CO_rate$Tot_recombinants/(0.0668807*tab_correl_extrapolation_CO_rate$Tot_fragments)*0.119 * 10**5*2
tab_correl_extrapolation_CO_rate$CO_rate_inf_fragments <- tab_correl_extrapolation_CO_rate$Nb_COs/(tab_correl_extrapolation_CO_rate$Nb_inf_fragments*500)*10**2*10**6*2



# Plot --------------------------------------------------------------------

# Open output file
png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

# Plot
ggplot(tab_correl_extrapolation_CO_rate) +
    geom_point(aes(x = CO_rate_ABC, y = CO_rate_inf_fragments), size=0.8) +
    labs(x = "CO rate (cM/Mb) extrapolated from the ABC",
         y = "CO rate (cM/Mb) extrapolated from informative fragments")+
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=30)) +
    geom_abline(intercept = lm(tab_correl_extrapolation_CO_rate$CO_rate_ABC ~ tab_correl_extrapolation_CO_rate$CO_rate_inf_fragments)$coeff[[1]], 
                slope = lm(tab_correl_extrapolation_CO_rate$CO_rate_ABC ~ tab_correl_extrapolation_CO_rate$CO_rate_inf_fragments)$coeff[[2]], 
                linetype=3) #+
#geom_abline(intercept = 0, slope =1, linetype=1)

# Close output file
dev.off()

