source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "supp-figure-A.2"
width_png <- 1700*0.75
height_png <- 1320*0.75



# Get data ----------------------------------------------------------------

# Get table
table_donor <- read.table("./data/ALL.Fragments.givers_with_DMC1_SUM_59_60.txt", 
                          header = TRUE)

# Proportion of CAST donor. NB: half weight for Rec1S because, for each recombination event, 2 Rec1S can be detected VS only 1 Rec2S.
table_donor$prop_CAST_donor <- ((1/2) * table_donor$Nb_C_giver_in_simple + table_donor$Nb_C_giver_in_complex) / (table_donor$Nb_C_giver_in_complex + table_donor$Nb_B_giver_in_complex + (1/2) * (table_donor$Nb_C_giver_in_simple + table_donor$Nb_B_giver_in_simple))

# Select those that have a minimum of 5 recombinants with donor identified
val <- 5
table_donor_nb_recomb_ok <- table_donor[which(table_donor$Nb_B_giver_in_simple + table_donor$Nb_C_giver_in_simple + table_donor$Nb_B_giver_in_complex + table_donor$Nb_C_giver_in_complex  >= val),]

# Add colors based on PRDM9 target
P9_CAST <- c("tC.hB", "NOV.tC.hB", "tC.hB+C")
P9_B6 <- c("tB.hC", "NOV.tB.hC", "tB.hB+C")
table_donor_nb_recomb_ok <- within(table_donor_nb_recomb_ok, {
    PRDM9_target = ifelse(Hotspot_category %in% P9_B6, "B6", 
                          ifelse(Hotspot_category %in% P9_CAST, "CAST", 
                                 "NA"))
    })


# Plot --------------------------------------------------------------------

# Open output file
png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

# Plot
ggplot(table_donor_nb_recomb_ok[which(table_donor_nb_recomb_ok$PRDM9_target != "NA"),]) +
    geom_point(aes(x = NA., y = prop_CAST_donor, color = PRDM9_target, shape = PRDM9_target), size = 3) +
    labs(x = "Expected proportion of CAST-donor fragments (based on DMC1 SSDS)",
         y = "Observed proportion of CAST-donor fragments") +
    scale_colour_manual("PRDM9 target",
                        values=c(color_B6_2, color_CA_2),
                        labels = c("B6", "CAST")) +
    scale_shape_manual("PRDM9 target", values=c(8, 19)) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=35),
          legend.position=c(0.14,0.92),
          legend.key.size = unit(1.5, 'lines')) +
    geom_abline(intercept = 0, slope = 1, linetype=2)

# Close output file
dev.off()

