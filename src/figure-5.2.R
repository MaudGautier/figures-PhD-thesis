source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-5.2"
width_png <- 2165*0.75
height_png <- 1240*0.75



# Get data ----------------------------------------------------------------

file_recomb_per_hotspot_with_chr_size <- './data/Recombinants_per_hotspot.txt'
table_recomb_per_hotspot_with_chr_size <- read.table(file_recomb_per_hotspot_with_chr_size, header=T)

chrOrder<-c(paste("chr",1:19,sep=""),"chrX","chrY","chrM")
table_recomb_per_hotspot_with_chr_size$Hotspot_chr <- factor(table_recomb_per_hotspot_with_chr_size$Hotspot_chr, levels=chrOrder)



# Plot --------------------------------------------------------------------

# Open output file
png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

# Plot
ggplot() +
    geom_segment(data = table_recomb_per_hotspot_with_chr_size,
                 aes(x = 5, xend = 5, y = 0, yend = Chrom_size),
                 lineend = "round", color = "lightgrey", size = 5) +
    facet_grid(. ~ table_recomb_per_hotspot_with_chr_size$Hotspot_chr) +
    geom_segment(data = table_recomb_per_hotspot_with_chr_size,
                 aes(x = 0, xend = nb_recomb,
                     y = Hotspot_start+(Hotspot_start - Hotspot_end)/2, yend = Hotspot_start+(Hotspot_start - Hotspot_end)/2),
                 size = 1,
                 inherit.aes=FALSE) +
    theme_bw() +
    theme(text=element_text(family="LM Roman 10", size = 35), 
          axis.title.y = element_text(margin = margin(l = 10, r = 10)), 
          axis.title.x = element_text(margin = margin(t = 10, b = 10))) +
    labs(y="Position on chromosome (Mb)", x="Number of recombination events (x 100)") +
    scale_y_continuous(labels=formatter1000000) +
    scale_x_continuous(labels=formatter100)

# Close output file
dev.off()

