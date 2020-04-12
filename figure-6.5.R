source("libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.5"
width_png <- 1254
height_png <- 974



# Get data ----------------------------------------------------------------


table_donor <- read.table("./data/ALL.Fragments.givers_with_DMC1_SUM_59_60.txt", header = TRUE)
val <- 5
table_donor_nb_recomb_ok <- table_donor[which(table_donor$Nb_B_giver_in_simple + table_donor$Nb_C_giver_in_simple + table_donor$Nb_B_giver_in_complex + table_donor$Nb_C_giver_in_complex  >= val),]


# Plot --------------------------------------------------------------------

png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

ggplot(table_donor_nb_recomb_ok) +
    geom_point(aes(x = NA., y = ((1/2)*Nb_C_giver_in_simple + Nb_C_giver_in_complex)/(Nb_C_giver_in_complex + Nb_B_giver_in_complex + (1/2)*(Nb_C_giver_in_simple + Nb_B_giver_in_simple)))) +
    labs(x = "Expected proportion of CAST-donor fragments (based on DMC1 SSDS)",
         y = "Observed proportion of CAST-donor fragments")+
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=35)) +
    geom_abline(intercept = 0, slope = 1, linetype=2)


dev.off()



## POUR INFO - figure coloree en supplement ?


# table_donor_nb_recomb_ok_with_colors <- table_donor_nb_recomb_ok
# 
# P9_CAST <- c("tC.hB", "NOV.tC.hB", "tC.hB+C")
# P9_B6 <- c("tB.hC", "NOV.tB.hC", "tB.hB+C")
# 
# table_donor_nb_recomb_ok_with_colors <- within(table_donor_nb_recomb_ok_with_colors, {
#     PRDM9_target = ifelse(Hotspot_category %in% P9_B6, "B6", ifelse(Hotspot_category %in% P9_CAST, "CAST", "NA"))
# })    

# Idem sans les unknown (=PLOT FINAL FINAL COLORE) - avec les bon label pour y ET SUR L'ORDI du labo
# ggplot(table_donor_nb_recomb_ok_with_colors[which(table_donor_nb_recomb_ok_with_colors$PRDM9_target != "NA"),]) +
#     geom_point(aes(x = NA., y = ((1/2)*Nb_C_giver_in_simple + Nb_C_giver_in_complex)/(Nb_C_giver_in_complex + Nb_B_giver_in_complex + (1/2)*(Nb_C_giver_in_simple + Nb_B_giver_in_simple)), color = PRDM9_target, shape = PRDM9_target), size = 1.5) +
#     labs(x = "Expected proportion of CAST-donor fragments (based on DMC1 SSDS)",
#          y = "Observed proportion of CAST-donor fragments") +
#     scale_colour_manual("PRDM9 target", values=c("#993333", "#FFCC00"), labels = c("B6", "CAST")) +
#     scale_shape_manual("PRDM9 target", values=c(8, 19)) +
#     scale_x_continuous(labels = scales::percent) +
#     scale_y_continuous(labels = scales::percent) +
#     theme_classic() +
#     theme(text=element_text(family="LM Roman 10", size=25), 
#           legend.position=c(0.14,0.92), 
#           legend.key.size = unit(1.5, 'lines')) + 
#     geom_abline(intercept = 0, slope = 1, linetype=2)

