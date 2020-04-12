source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-7.5"
width_png <- 1616*0.9
height_png <- 1165*0.9



# Get data ----------------------------------------------------------------

tab_rel_density_prop_NCO1 <- read.table("./data/NCO1_summary.txt", header = T)



# Plot --------------------------------------------------------------------

png(paste("./output/", fig_name, ".png", sep = ""), 
    width = width_png, height = height_png)

ggplot(tab_rel_density_prop_NCO1) + 
    geom_point(aes(y=PropNCO1, x=Divergence), size=3) + 
    theme_classic() +
    labs(y = "Proportion of NCO CT markers involved in NCO-1 events",
         x = "SNP density") +
    scale_y_continuous(labels = scales::percent) +
    theme(text=element_text(family="LM Roman 10", size = 35),
          axis.title.y = element_text(margin = margin(r = 20, l=20)), 
          axis.title.x = element_text(margin = margin(t = 20, b=10)))

dev.off()

