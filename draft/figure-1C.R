source("libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-1C"
width_png <- 1400*0.75
height_png <- 1184*0.75



# Get data ----------------------------------------------------------------

tab_paigen <- read.table("./data/PositiveControls_with_all_peaks_compared_to_nb_events_detected_with_informative_fragments_Summed_nb_events_on_interval.txt", header=T)
tab_paigen$CO_nb_corrected_for_paigen <- tab_paigen$Tot_recombinants/tab_paigen$LgInterval_kb

tab_paigen$CO_corrected_log10  <- log10(tab_paigen$RR_Male_cM_Mb)
tab_paigen$RR_Male_cM_Mb_log10 <- log10(tab_paigen$RR_Male_cM_Mb)

# Plot --------------------------------------------------------------------

png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

lr <- lm(log10(tab_paigen$CO_nb_corrected_for_paigen) ~ log10(tab_paigen$RR_Male_cM_Mb))


ggplot(tab_paigen) +
    geom_point(aes(x = RR_Male_cM_Mb, y = CO_nb_corrected_for_paigen), size=2) +
    labs(x = expression(paste('CO rate measured by Paigen et al. (2008) (cM/Mb) (log-transformed)')),
         y = expression(paste('Number of events per kb (log-transformed)'))) +
    theme_classic() +
    scale_x_log10(breaks=trans_breaks("log10", function(x) 10^x), 
                  labels = function(x) format(x, scientific = FALSE)) +
    scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x)) +
    theme(text=element_text(family="LM Roman 10", size=30)) +
    geom_abline(intercept = (lr$coeff[[1]]), 
                slope = (lr$coeff[[2]]), 
                linetype=3)

dev.off()



# Notes:
# Régression linéaire log-log: log(y) = f(log(x))

# Droite loin des points => Vérification à la main du calcul des coeff regression linéaire
## Slope
# x1 ~ 1   ; y1 ~ 0.1
# x2 ~ 100 ; y2 ~ 10
# => slope = (y2-y1)/(x2-x1) ~ (10-0.1)/(100-1) ~ 0.1

## Intercept
# y1 = slope * x1 + intercept
# => intercept = y1 - slope * x1
# => intercept ~ 0.1 - 0.1 * 1 
# => intercept ~ 0

## Conclusion:
# Résultats differents de la regression lineaire calculee par lm
# CAR: beaucoup de points sont tout à gauche (valeur RR_cM_Mb très faibles)
# (Visible quand on met l'option limits = c(0.00000000000000000001, 1000) dans scale_x_log10)
# Donc il faudrait une relation NON-lineaire si on veut plotter la relation

# Tentative de correction en utilisant log1p (i.e. log10(1+x)) au lieu de log10
# Dans ce cas, la transformation des axes devient : 10**x - 1 (et plus 10**x)
# Démonstration:
# log10(1+x) = z
# => 10**log10(1+x) = 10**z
# => 1 + x = 10^z
# => x = 10^z - 1

# Mais ne résout pas le problème



## TODO:
# Trouver quelle relation (log-log) à plotter
# OU choisir de ne pas plotter la régression et seulement afficher les points


