source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.7"
width_png <- 1550*0.75
height_png <- 1227*0.75



# Get data ----------------------------------------------------------------

## Get table
tab_correl_extrapolation_CO_rate <- read.table("./data/Recombination_events_per_hotspot_with_informative_fragments.txt", 
                                               header=T, comment="")

## Extrapolate CO rate

## Method 1 (from the ABC)
# CO rate = nb_events / (detectability * nb_fragments) * perc_CO * 1000 * 2 * 100 
# NB: 1000 CO rate measured on 1-kb hotspots, 2 because bivalents to get cM, 100 because cM.
tab_correl_extrapolation_CO_rate$CO_rate_ABC <- tab_correl_extrapolation_CO_rate$Tot_recombinants/(0.0668807 * tab_correl_extrapolation_CO_rate$Tot_fragments) * 0.119 * 10**5 * 2

## Method 2 (from the number of 'informative' fragments, i.e. fragments having at least 4 markers, i.e. fragments that can be detected as COs)
# CO rate = nb_COs / (nb_informative_fragments * length_sequenced_per_fragment) * 100 * 1000000 * 2
# NB: 100 because cM, 1000000 because Mb, 2 because bivalent to get cM
tab_correl_extrapolation_CO_rate$CO_rate_inf_fragments <- tab_correl_extrapolation_CO_rate$Nb_COs/(tab_correl_extrapolation_CO_rate$Nb_inf_fragments*500) * 10**2 * 10**6 * 2



# Plot --------------------------------------------------------------------

# Open output file
png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

# Get correlation
lr <- lm(tab_correl_extrapolation_CO_rate$CO_rate_ABC ~ tab_correl_extrapolation_CO_rate$CO_rate_inf_fragments)

# Plot
ggplot(tab_correl_extrapolation_CO_rate) +
    geom_point(aes(x = CO_rate_ABC, y = CO_rate_inf_fragments), size = 0.8) +
    labs(x = "CO rate (cM/Mb) extrapolated from the ABC",
         y = "CO rate (cM/Mb) extrapolated from informative fragments")+
    theme_classic() +
    theme(text=element_text(family = "LM Roman 10", size = 30)) +
    geom_abline(intercept = lr$coeff[[1]], 
                slope     = lr$coeff[[2]], 
                linetype = 3) 

# Close output file
dev.off()

