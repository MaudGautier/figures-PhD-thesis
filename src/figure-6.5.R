source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.5"
width_png <- 1254
height_png <- 974



# Get data ----------------------------------------------------------------

# Get table
table_donor <- read.table("./data/ALL.Fragments.givers_with_DMC1_SUM_59_60.txt", 
                          header = TRUE)

# Proportion of CAST donor. NB: half weight for Rec1S because, for each recombination event, 2 Rec1S can be detected VS only 1 Rec2S.
table_donor$prop_CAST_donor <- ((1/2) * table_donor$Nb_C_giver_in_simple + table_donor$Nb_C_giver_in_complex) / (table_donor$Nb_C_giver_in_complex + table_donor$Nb_B_giver_in_complex + (1/2) * (table_donor$Nb_C_giver_in_simple + table_donor$Nb_B_giver_in_simple))

# Select those that have a minimum of 5 recombinants with donor identified
val <- 5
table_donor_nb_recomb_ok <- table_donor[which(table_donor$Nb_B_giver_in_simple + table_donor$Nb_C_giver_in_simple + table_donor$Nb_B_giver_in_complex + table_donor$Nb_C_giver_in_complex  >= val),]


# Plot --------------------------------------------------------------------

# Open output file
png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

# Plot
ggplot(table_donor_nb_recomb_ok) +
    geom_point(aes(x = NA., y = prop_CAST_donor)) +
    labs(x = "Expected proportion of CAST-donor fragments (based on DMC1 SSDS)",
         y = "Observed proportion of CAST-donor fragments")+
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(text = element_text(family = "LM Roman 10", size=35)) +
    geom_abline(intercept = 0, slope = 1, linetype=2)

# Close output file
dev.off()
