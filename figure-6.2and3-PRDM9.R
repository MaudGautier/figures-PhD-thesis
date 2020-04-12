## TODO: Rewrite script to simplify (function for DMC1 and DMC2 ?)

source("libraries.R")


# Parameters --------------------------------------------------------------

fig_name_2 <- "figure-6.2-PRDM9"
fig_name_3 <- "figure-6.3-PRDM9"
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
table_intensity_sorted_intensity <- table_intensity[order(table_intensity$intensity),]

# 3.2. Get mean values for each group
table_groups <- cbind(table_intensity_sorted_intensity, vec_groups)
mean_val <- list()
for (i in unique(vec_groups)) {
    mean_val[i] <- median(table_groups[which(table_groups$vec_groups == i),]$intensity, na.rm=T)
}

# 3.3. Create vector of mean values and add it to each group
vec_group_means <- c()
for (i in seq(1, nb_groups)) {
    vec_group_means <- c(vec_group_means,
                               rep(mean_val[[i]], ceiling(nrow(table_intensity)/nb_groups)))
}
vec_group_means <- vec_group_means[1:nrow(table_intensity)]
table_groups <- cbind(table_intensity_sorted_intensity, vec_groups, vec_group_means)


table_groups$ratio <- table_groups$X1 / table_groups$Detectability


table_groups$COrate <- table_groups$X1/(0.0668807*table_groups$Tot_nb_fragments)*100*10**6/500 * 0.119 *2
sumup_table_groups_std <- aggregate(.~vec_group_means, data = table_groups, std)
# Then get means
sumup_table_groups_mean <- aggregate(.~vec_group_means, data = table_groups, mean)




# Plot --------------------------------------------------------------------


png(paste("./output/", fig_name_2, ".png", sep = ""), 
    width = width_png, height = height_png)

ggplot(sumup_table_groups_mean) +
    geom_errorbar(aes(as.numeric(vec_group_means), ymax=sumup_table_groups_mean$COrate +sumup_table_groups_std$COrate , ymin=sumup_table_groups_mean$COrate -sumup_table_groups_std$COrate)) +
    geom_point(aes(as.numeric(vec_group_means), y = sumup_table_groups_mean$COrate)) +
    scale_x_continuous() +
    scale_y_continuous(limits = quantile(table_groups$X1/(0.0668807*table_groups$Tot_nb_fragments)*10**5 * 0.119*2, c(.00001, 0.99), na.rm=T), 
                       sec.axis = sec_axis(~.*0.02297+0.21222, name = "# Recombination events / sequenced Mb"), position="right") +
    labs(y = "CO Rate (cM/Mb)",
         x = "Median PRDM9 signal") +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=28),
          axis.title.y.right = element_text(margin = margin(l = 10, r = 10)),
          axis.title.y = element_text(margin = margin(l = 10, r = 10)), 
          axis.title.x = element_text(margin = margin(t = 10, b = 10)))


dev.off()






# Asymmetric --------------------------------------------------------------



# model_1bis_PRDM9 <- lm(X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119 *2 ~ vec_group_means_PRDM9, data = table_groups_PRDM9[which(table_groups_PRDM9$class %in% c("tC.hB+C", "tB.hB+C", "tC.hB", "tB.hC")),])
# model_2bis_PRDM9 <- lm(X1/(0.0668807*Tot_nb_fragments)*100000 * 0.119 *2 ~ vec_group_means_PRDM9, data = table_groups_PRDM9[which(table_groups_PRDM9$class %in% c("NOV.tC.hB", "NOV.tB.hC")),])

model_1 <- lm(X1/(Detectability/1000000) ~ vec_group_means, data = table_groups[which(table_groups$class %in% c("tC.hB+C", "tB.hB+C", "tC.hB", "tB.hC")),])
model_2 <- lm(X1/(Detectability/1000000) ~ vec_group_means, data = table_groups[which(table_groups$class %in% c("NOV.tC.hB", "NOV.tB.hC")),])


## GRAPHE PRDM9 # OK MAIS SANS LEGENDE ET OUTLIERS
displayecement=0.0007
width_bar=displayecement*2
col1 <- brewer.pal(n=4, name='Accent')[1]
col2 <- brewer.pal(n=4, name='Accent')[3]



png(paste("./output/", fig_name_3, ".png", sep = ""), 
    width = width_png, height = height_png)


# PRDM9 With legend - avec les outliers
ggplot() +
    geom_abline(slope = model_1$coefficients[[2]], intercept = model_1$coefficients[[1]], col = col1, lty = 2, alpha = 0.7) + 
    geom_abline(slope = model_2$coefficients[[2]], intercept = model_2$coefficients[[1]], col = col2, lty = 2, alpha = 0.8) + 
    geom_boxplot(aes(x= as.numeric(vec_group_means) - displayecement,
                     y =X1/(Detectability/1000000),
                     group = as.numeric(vec_group_means),
                     fill = col2),
                 lwd = 0.75, width = displayecement*2, alpha = 0.8, outlier.color = col2, outlier.size = 2, 
                 data = table_groups[which(table_groups$class %in% c("NOV.tC.hB", "NOV.tB.hC")),]) +
    geom_boxplot(aes(x= as.numeric(vec_group_means) + displayecement,
                     y =X1/(Detectability/1000000),
                     group = as.numeric(vec_group_means),
                     fill = col1),
                 lwd = 0.75,width = displayecement*2, alpha = 0.7, outlier.color = col1, outlier.size = 2, 
                 data = table_groups[which(table_groups$class %in% c("tC.hB+C", "tB.hB+C", "tC.hB", "tB.hC")),]) +
    scale_x_continuous() +
    scale_y_continuous(sec.axis = sec_axis(~.*(sum(as.numeric(table_groups$Detectability), na.rm=T)/1000000) / (0.0668807*sum(as.numeric(table_groups$Tot_nb_fragments), na.rm=T)) *100000 * 0.119*2, name = "CO Rate (cM/Mb)")) +
    labs(y = "# Recombination events / sequenced Mb",
         x = paste0("Median PRDM9 signal")) +
    theme_classic() +
    theme(text=element_text(family="LM Roman 10", size=28),
          axis.title.y.right = element_text(margin = margin(l = 10, r = 10)),
          axis.title.y = element_text(margin = margin(l = 10, r = 10)), 
          axis.title.x = element_text(margin = margin(t = 10, b = 10)),
          legend.position = c(0.5, 0.90), 
          legend.key.size = unit(2, 'lines')) +
    scale_fill_identity(guide = "legend", name="Hotspot type", breaks = c(col1, col2), labels = c("Symmetric", "Asymmetric"))


dev.off()
