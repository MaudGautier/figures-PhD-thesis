

file_detectability_across_hotspots <- '/Users/maudgautier/Documents/These/_RANGE/Copy_cluster_tmp/2_Data_mouse/2_dBGC/09_Statistics_recombinants_REDONE/Recombination_rates/FINAL_DETECTABILITY_ACROSS_HOTSPOTS'
table_detectability_across_hotspots <- read.table(file_detectability_across_hotspots, header = T)
head(table_detectability_across_hotspots)

ggplot(table_detectability_across_hotspots) +
    geom_rect(mapping=aes(xmin = Position, xmax = Position, ymin = 0, ymax=Fragments), color = "lightgrey") +
    geom_rect(mapping=aes(xmin = Position, xmax = Position, ymin = 0, ymax=Recombinants*5000), color = "grey51") +
    theme_classic() + 
    geom_vline(aes(xintercept=0), linetype=2, color ="darkgrey") +
    geom_vline(aes(xintercept=500), linetype=3, color ="darkgrey") +
    geom_vline(aes(xintercept=-500), linetype=3, color ="darkgrey") +
    scale_x_continuous("Position (Relative to the peak centre)", limits=c(-1000, 1000)) +
    scale_y_continuous("Number of Fragments (light grey)", sec.axis = sec_axis(~ . /5000, name = "Number of Recombinants (dark grey)"), labels = scales::comma) +
    theme(text=element_text(family="LM Roman 10", size = 15)) +
    labs(title = "Distribution of Recombinants and fragments along all hotspots")







