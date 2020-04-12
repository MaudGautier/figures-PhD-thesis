source("libraries.R")


# Parameters --------------------------------------------------------------

fig_name_2 <- "figure-6.2-DMC1"
fig_name_3 <- "figure-6.3-DMC1"
width_png <- 1200*0.75
height_png <- 949*0.75



# Get data ----------------------------------------------------------------


table_intensity <- read.table("./data/Intensity_and_counts_good_classes_and_DMC1_with_nb_frags_per_hotspot.txt", header = TRUE)


# 3.1. Create 10 groups and add column to the table
nb_groups <- 10
vec_groups <- NULL
for (i in seq(1, nb_groups)) {
    
    group_vals <- rep(i, ceiling(nrow(table_intensity)/nb_groups))
    vec_groups <- c(vec_groups, group_vals)
    
}
vec_groups <- vec_groups[1:nrow(table_intensity)]
table_intensity_sorted_intensity <- table_intensity[order(table_intensity$DMC1_B6+table_intensity$DMC1_CAST),]

# 3.2. Get mean values for each group
table_groups <- cbind(table_intensity_sorted_intensity, vec_groups)
mean_val <- list()
for (i in unique(vec_groups)) {
    mean_val[i] <- median(table_groups[which(table_groups$vec_groups == i),]$DMC1_B6+table_groups[which(table_groups$vec_groups == i),]$DMC1_CAST, na.rm=T)
}

# 3.3. Create vector of mean values and add it to each group
vec_group_means <- c()
for (i in seq(1, nb_groups)) {
    vec_group_means <- c(vec_group_means,
                         rep(mean_val[[i]], ceiling(nrow(table_intensity)/nb_groups)))
}
vec_group_means <- vec_group_means[1:nrow(table_intensity)]
table_groups <- cbind(table_intensity_sorted_intensity, vec_groups, vec_group_means)






#table_groups$COrate <- table_groups$X1/(0.0668807*table_groups$Tot_nb_fragments)*100000 * 0.119
# CORRECTION LAURENT X 2 (CO rate = nbrecombinants / (detectability * nb_frags) * 100000 * propCO * 2) # 100.000 = 100 * 1000 (100 pour cM, 1000 pour hotspots de 1 kb mesures)
# voir aussi: calcul du taux de recombinaison moyen across hotspots: 18821/(0.0668807*228984512)*1000*0.119*2*100 = 29.24912 cm/Mb
table_groups$COrate <- table_groups$X1/(0.0668807*table_groups$Tot_nb_fragments)*100000 * 0.119*2

sumup_table_groups_std <- aggregate(.~vec_group_means, data = table_groups, std)
# Then get means
sumup_table_groups_mean <- aggregate(.~vec_group_means, data = table_groups, mean)


# Plot --------------------------------------------------------------------


# png(paste("./output/", fig_name, ".png", sep = ""), 
#     width = width_png, height = height_png)

ggplot(sumup_table_groups_mean) +
    geom_errorbar(aes(as.numeric(vec_group_means), ymax=sumup_table_groups_mean$COrate +sumup_table_groups_std$COrate , ymin=sumup_table_groups_mean$COrate -sumup_table_groups_std$COrate)) +
    geom_point(aes(as.numeric(vec_group_means), y = sumup_table_groups_mean$COrate)) +
    geom_abline(intercept = lm(sumup_table_groups_mean$COrate ~ sumup_table_groups_mean$vec_group_means)$coeff[[1]],
                slope = lm(sumup_table_groups_mean$COrate ~ sumup_table_groups_mean$vec_group_means)$coeff[[2]], linetype=3) +
    scale_x_continuous() +
    scale_y_continuous(limits = quantile(table_groups$X1/(0.0668807*table_groups$Tot_nb_fragments)*100000 * 0.119*2, c(.00001, 0.9), na.rm=T), 
                       sec.axis = sec_axis(~.*0.02297+0.21222, name = "# Recombination events / sequenced Mb"), position="right") +
    labs(y = "CO Rate (cM/Mb)",
         x = "Median DMC1 signal") +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=28),
          axis.title.y.right = element_text(margin = margin(l = 10, r = 10)),
          axis.title.y = element_text(margin = margin(l = 10, r = 10)), 
          axis.title.x = element_text(margin = margin(t = 10, b = 10)))


# dev.off()





# Asymmetric --------------------------------------------------------------


col1 <- brewer.pal(n=4, name='Accent')[1]
col2 <- brewer.pal(n=4, name='Accent')[3]
displayecement = 14

# ATTENTION: DANS CE CAS IL FAUDRAIT AUSSI REFAIRE LES DROITES
model_1bis <- lm(X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119*2 ~ vec_group_means, data = table_groups[which(table_groups$class %in% c("tC.hB+C", "tB.hB+C", "tC.hB", "tB.hC")),])
model_2bis <- lm(X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119*2 ~ vec_group_means, data = table_groups[which(table_groups$class %in% c("NOV.tC.hB", "NOV.tB.hC")),])
model_allbis <- lm(X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119 ~ vec_group_means, data = table_groups[which(table_groups$class %in% c("NOV.tC.hB", "NOV.tB.hC", "tC.hB+C", "tB.hB+C", "tC.hB", "tB.hC")),])



### IDEM POUR REDACTION MAIS CHANGE JUSTE LE LABEL X POUR APPPENDIX MANUSCRIT - Corrected CO X 2
ggplot() +
    geom_abline(slope = model_1bis$coefficients[[2]], intercept = model_1bis$coefficients[[1]], col = col1, lty = 2, alpha = 0.7) + 
    geom_abline(slope = model_2bis$coefficients[[2]], intercept = model_2bis$coefficients[[1]], col = col2, lty = 2, alpha = 0.8) + 
    geom_boxplot(aes(x= as.numeric(vec_group_means) - displayecement,
                     y =X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119*2,
                     group = as.numeric(vec_group_means),
                     fill = col2),
                 lwd = 0.75, width = displayecement*2, alpha = 0.8, outlier.color = col2, outlier.size = 2, 
                 data = table_groups[which(table_groups$class %in% c("NOV.tC.hB", "NOV.tB.hC")),]) +
    geom_boxplot(aes(x= as.numeric(vec_group_means) + displayecement,
                     y =X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119*2,
                     group = as.numeric(vec_group_means),
                     fill = col1),
                 lwd = 0.75,width = displayecement*2, alpha = 0.7, outlier.color = col1, outlier.size = 2, 
                 data = table_groups[which(table_groups$class %in% c("tC.hB+C", "tB.hB+C", "tC.hB", "tB.hC")),]) +
    scale_x_continuous() +
    scale_y_continuous(sec.axis = sec_axis(~.*0.02297+0.21222, 
                                           breaks = c(0,2,4,6,8),
                                           name = "# Recombination events / sequenced Mb"), position="right") +
    labs(y = "CO Rate (cM/Mb)",
         x = "Median DMC1 signal") +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=28),
          axis.title.y.right = element_text(margin = margin(l = 10, r = 10)),
          axis.title.y = element_text(margin = margin(l = 10, r = 10)), 
          axis.title.x = element_text(margin = margin(t = 10, b = 10)),
          legend.position = c(0.51, 0.90), 
          legend.key.size = unit(1.8, 'lines')) +
    scale_fill_identity(guide = "legend", name="Hotspot type", breaks = c(col1, col2), labels = c("Symmetric", "Asymmetric"))









