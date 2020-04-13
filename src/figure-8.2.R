source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-8.2"
width_png <- 2165*0.77
height_png <- 1240*0.77



# Get data ----------------------------------------------------------------

table_HFM1_bakgrounds_per_hotspot_with_chr_size <- read.table("./data/HFM1_List_genotypes_hotspots_refined_with_chr_size.txt", header=T, comment="")

# Reorder chromosomes
chrOrderbis<-c(paste("chr",1:19,sep=""),"chrX","chrY","chrM")
table_HFM1_bakgrounds_per_hotspot_with_chr_size$Chr <- factor(table_HFM1_bakgrounds_per_hotspot_with_chr_size$Chr, levels=chrOrderbis)

# Replace NOCOV with NA
table_HFM1_bakgrounds_per_hotspot_with_chr_size[table_HFM1_bakgrounds_per_hotspot_with_chr_size=="NOCOV"] <- NA



# Plot --------------------------------------------------------------------

for (sample in c("S28353", "S28355", "S28367", "S28371")) {
    
    # Open output file
    png(paste("./output/", fig_name, "-", sample, ".png", sep = ""), 
        width = width_png, height = height_png)
    
    # Plot
    p <- ggplot() +
        geom_segment(data = table_HFM1_bakgrounds_per_hotspot_with_chr_size,
                     aes(x = 5, xend = 5, y = 0, yend = Chrom_size),
                     lineend = "round", color = "grey88", size = 10) +
        facet_grid(. ~ table_HFM1_bakgrounds_per_hotspot_with_chr_size$Chr) +
        geom_segment(data = table_HFM1_bakgrounds_per_hotspot_with_chr_size,
                     aes(x = 3.1, xend = 6.9,
                         y = Start+(Start - Stop)/2, yend = Start+(Start - Stop)/2,
                         colour = get(sample)), # column name == sample
                     size = 1,
                     inherit.aes=FALSE) +
        theme_minimal() +
        theme(text = element_text(family="LM Roman 10", size = 35), 
              axis.title.y = element_text(margin = margin(l = 10, r = 10)), 
              axis.title.x = element_text(margin = margin(t = 10, b = 10)),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(), 
              legend.key.size = unit(2, 'lines')) +
        labs(y = "Position on chromosome (Mb)", x = "") +
        scale_y_continuous(labels = formatter1000000) +
        scale_x_continuous(limits = c(0, 10)) +
        scale_colour_manual('Background', 
                            values = c(color_DOM, color_HET), 
                            labels = c("DOM/DOM", "DOM/CAST"), 
                            na.value = "grey50")
    print(p)
    
    # Close output file
    dev.off()
}

