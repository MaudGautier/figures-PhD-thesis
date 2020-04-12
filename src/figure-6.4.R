## TODO: Rewrite script to simplify (function for DMC1 and DMC2 ?)

source("src/libraries.R")


# Parameters --------------------------------------------------------------

fig_name <- "figure-6.4"
width_png <- 941*0.75
height_png <- 676*0.75
### WARNING: TO DO ON OTHER MAC: CHECK SIZES


# Get data ----------------------------------------------------------------




Nsnps <- 10
N1 <- 7
N2 <- 7
N3 <- 6
N4 <- 5
N5 <- 0
tab_test <- structure(list(Hotspot = rep("Hotspot", N1+N2++N3+N4+N5),
                           DSB = rep(1500, N1+N2++N3+N4+N5),
                           Read_ID = structure(c(rep(1L, N1), rep(2L, N2), rep(3L, N3), rep(4L, N4)#, rep(5L, N5)
                           ),
                           .Label = c("R1","R2", "R3", "R4"), class = "factor"),
                           Status = c(rep("CO", N1), rep("CO", N2), rep("NCO", N3), rep("non", N4)#, rep("non", N5)
                           ),
                           Start = c(1200L, 1326L, 1350L, 1362L, 1378L, 1512L, 1520L, #1570L, 1800L, 1950L, 
                                     1368L, 1378L, 1512L, 1520L, 1570L, 1800L, 1950L, 
                                     1508L, 1512L, 1520L, 1570L, 1800L, 1950L, 
                                     1250L, 1326L, 1350L, 1362L, 1378L#, 
                                     #1250L, 1326L, 1350L, 1362L, 1378L, 1512L, 1520L, 1570L, 1800L, 1950L
                           ),
                           Stop =  c(1326L, 1350L, 1362L, 1378L, 1512L, 1520L, 1587L, #1800L, 1950L, 2003L,
                                     1378L, 1512L, 1520L, 1570L, 1800L, 1950L, 1972L,
                                     1512L, 1520L, 1570L, 1800L, 1950L, 2003L,
                                     1326L, 1350L, 1362L, 1378L, 1480L #,
                                     #1326L, 1350L, 1362L, 1378L, 1512L, 1520L, 1570L, 1800L, 1950L, 2003L
                           ),
                           Tract_start = c(1433.5L, NA, NA, NA, NA, NA, NA, #NA, NA, NA,
                                           1500L, NA, NA, NA, NA, NA, NA,
                                           1545L, NA, NA, NA, NA, NA,
                                           NA, NA, NA, NA, NA#,
                                           #1433.5L, NA, NA, NA, NA, NA, NA, NA, NA, NA
                           ),
                           Tract_end = c(1500L, NA, NA, NA, NA, NA, NA, #NA, NA, NA,
                                         1685L, NA, NA, NA, NA, NA, NA,
                                         1875L, NA, NA, NA, NA, NA,
                                         NA, NA, NA, NA, NA#,
                                         #1500L, NA, NA, NA, NA, NA, NA, NA, NA, NA
                           ),
                           Genotype = structure(c(1L, 1L, 1L, 1L, NA, 2L, 2L,# 2L, 2L, 2L,
                                                  1L, 1L, 1L, 1L, NA, 2L, 2L,
                                                  2L, 2L, NA, 1L, NA, 2L,
                                                  2L, 2L, 2L, 2L, 2L#, 
                                                  #1L, 1L, 1L, 1L, NA, NA, 2L, 2L, 2L, 2L
                           ), 
                           .Label = c("B6","CAST"), class = "factor"),
                           SNP_pos = c(NA, 1326L, 1350L, 1362L, 1378L, 1512L, 1520L, #1570L, 1800L, 1950L, 
                                       NA, 1378L, 1512L, 1520L, 1570L, 1800L, 1950L,
                                       NA, 1512L, 1520L, 1570L, 1800L, 1950L,
                                       NA, 1326L, 1350L, 1362L, 1378L#,
                                       #NA, 1326L, 1350L, 1362L, 1378L, 1512L, 1520L, 1570L, 1800L, 1950L
                           ),
                           SNP_genot = structure(c(NA, 1L, 1L, 1L, 1L, 2L, 2L, #2L, 2L, 2L,
                                                   NA, 1L, 1L, 1L, 1L, 2L, 2L,
                                                   NA, 2L, 2L, 1L, 1L, 2L,
                                                   NA, 2L, 2L, 2L, 2L#, 
                                                   #NA, 1L, 1L, 1L, NA, NA, 2L, 2L, 2L, 2L
                           ), 
                           .Label = c("B6","CAST"), class = "factor"),
                           SNP_type = structure(c(3L, 3L, 3L, 3L, 3L, 3L, 3L, #3L, 3L, 3L,
                                                  3L, 3L, 3L, 3L, 3L, 3L, 3L,
                                                  3L, 3L, 3L, 3L, 3L, 3L,
                                                  3L, 3L, 3L, 3L, 3L#, 
                                                  #3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L
                           ), 
                           .Label = c("DEL","INS","SNP"), class = "factor"),
                           start_first_NA_zone = c(1378L, 1378L, 1378L, 1378L, 1378L, 1378L, 1378L, #1378L, 1378L, 1378L, 
                                                   1378L, 1378L, 1378L, 1378L, 1378L, 1378L, 1378L,
                                                   1378L, 1378L, 1378L, 1378L, 1378L, 1378L,
                                                   1378L, 1378L, 1378L, 1378L, 1378L, 1378L#,
                                                   #1378L, 1378L, 1378L, 1378L, 1378L, 1378L, 1378L, 1378L, 1378L, 1378L
                           ),
                           middle_first_NA_zone = c(1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, #1433.5L, 1433.5L, 1433.5L,
                                                    1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L,
                                                    1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L,
                                                    1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L#,
                                                    #1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L, 1433.5L)
                           )),
                      .Names = c("Hotspot", "DSB", "Read_ID", "Status", "Start", "Stop", "Tract_start", "Tract_end", "Genotype", "SNP_pos", "SNP_genot", "SNP_type", "start_first_NA_zone", "middle_first_NA_zone"),
                      class = "data.frame", row.names = seq(1,N1+N2+N3+N4+N5))





df <- data.frame(
    ChIP_P9=round(c(rnorm(200, mean=1500, sd=110)))
)
# df1 <- df

densP9 <- density(df$ChIP_P9)
p1 <- ggplot(df, aes(x=-densP9$x[which.max(densP9$y)]+ChIP_P9, colour="grey50")) +
    geom_density() +
    theme_classic() +
    geom_vline(aes(xintercept=0), linetype=3, color ="darkgrey") +
    scale_x_continuous(limits=c(1200-1500,2003-1500)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          line = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = '', y = "") +
    scale_colour_manual(values = c("grey50" = "grey50")) 

tab_SNP_test <- structure(list(CHROM = rep("chrom", Nsnps),
                               POS = c(1326L, 1350L, 1362L, 1378L, 1489L, 1512L, 1520L, 1570L, 1800L, 1950L, 
                                       1326L, 1350L, 1362L, 1378L, 1489L, 1512L, 1520L, 1570L, 1800L, 1950L),
                               REF = structure(c(1L, 1L, 4L, 4L, 1L, 1L, 4L, 4L, 1L, 4L, 
                                                 1L, 1L, 4L, 4L, 1L, 1L, 4L, 4L, 1L, 4L), 
                                               .Label = c("A", "C", "G", "T"), class = "factor"),
                               ALT = structure(c(2L, 3L, 3L, 2L, 2L, 3L, 3L, 2L, 3L, 2L, 
                                                 2L, 3L, 3L, 2L, 2L, 3L, 3L, 2L, 3L, 2L), 
                                               .Label = c("A", "C", "G", "T"), class = "factor"),
                               PEAK_NAME = structure(c(rep(1L, 10), rep(2L, 10)), 
                                                     .Label = c("CAST", "B6"), class = "factor"),
                               DSB = rep(1500L, 20),
                               SNP_TYPE = structure(rep(3L, 20), .Label = c("DEL","INS","SNP"), class = "factor")),
                          .Names = c("CHROM", "POS", "REF", "ALT", "PEAK_NAME", "DSB", "SNP_TYPE"),
                          class = "data.frame", row.names = seq(1,2*Nsnps))
tab_SNP_test$PEAK_NAME2 <- factor(levels(tab_SNP_test$PEAK_NAME), c(rep("CAST",10), rep("B6",10)))






# Refait avec les bonnes couleurs (jaune/rouge)
p2_bis <- ggplot(tab_SNP_test) +
    geom_point(aes(x=POS-DSB, y = PEAK_NAME, shape = SNP_TYPE, size = def_sizes(SNP_TYPE), colour=PEAK_NAME)) + 
    scale_size_identity() +
    theme_classic() +
    scale_x_continuous(limits=c(1200-1500,2003-1500)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), line = element_blank(),
          axis.title.x=element_blank(), legend.position="none") + 
    labs(x = '', y = "") + 
    scale_colour_manual(values = c("B6" = "#A94343", "CAST" = "#EFBC00")) +
    annotate("text", x=-200, y=2.05, label= "B6 reference", colour = "#993333", hjust=1) + 
    annotate("text", x=-200, y=1.055, label= "CAST reference", colour = "#DFAC00", hjust=1) + 
    geom_vline(aes(xintercept=0), linetype=3, color ="darkgrey") 










p3_bis <- ggplot(tab_test, aes(colour=Genotype)) + 
    theme_classic() +
    geom_segment(mapping=aes(x=Tract_start-DSB, y=Read_ID, xend=Tract_end-DSB, yend=Read_ID, size=4), colour = "lightgrey") +
    geom_segment(mapping=aes(x=Start-DSB, y=Read_ID, xend=Stop-DSB, yend=Read_ID, size=0.7)) +
    geom_point(aes(x=SNP_pos-DSB, y = Read_ID, colour=SNP_genot, shape=SNP_type, size = def_sizes(SNP_type))) +
    scale_size_identity() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), line = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          #legend.position=c(0.94,0.1)
          legend.position="none") + 
    labs(title = "TERMINOLOGY", x = "", y = '') +
    annotate("text", x=(1685-1500)/2, y=1.80, label= "CO: B6 donor", colour = "grey50") + 
    annotate("text", x=(1433.5-1500)/2, y=0.80, label= "CO: Undefined donor", colour = "grey50") + 
    annotate("text", x=1545-1500+(1875-1545)/2, y=2.80, label= "NCO: B6 donor", colour = "grey50") + 
    annotate("text", x=(1480-1250)/2+(1250-1500), y=3.80, label= "Non-recombination event", colour = "grey50") + 
    #facet_grid(Status ~ ., scales = 'free', switch = "y", labeller=status_labeller_BIS) +
    #scale_shape_manual('Variant type', values = myTypes, labels = myLabels, drop = FALSE) +
    geom_vline(aes(xintercept=0), linetype=3, color ="darkgrey") + 
    
    annotate("segment", x = 1875-1500, y = 3.35, xend = 1875-1500, yend = 3.05, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    annotate("text", x=1875-1500, y=3.45, label= "Switch point = Switch interval midpoint", colour = "grey50") + 
    
    annotate("segment", x = 1875-1500, y = 2.875, xend = 1950-1500, yend = 2.875, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    annotate("segment", x = 1875-1500, y = 2.875, xend = 1800-1500, yend = 2.875, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    annotate("text", x=1875-1500, y=2.775, label= "Switch", colour = "grey50") + 
    annotate("text", x=1875-1500, y=2.625, label= "interval", colour = "grey50") + 
    
    # annotate("segment", x = 0, y = 2.15, xend = 0, yend = 2.025, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    # annotate("text", x=0, y=2.2, label= "Switch point = Inferred DSB position", colour = "grey50") + 
    annotate("segment", x = 0, y = 1.35, xend = 0, yend = 1.05, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    annotate("text", x=0, y=1.45, label= "Switch point = Inferred DSB position", colour = "grey50") + 
    
    # annotate("segment", x = 1466.75-1500, y = 1.05, xend = 1433.5-1500, yend = 1.05, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    # annotate("segment", x = 1466.75-1500, y = 1.05, xend = 1500-1500, yend = 1.05, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    # annotate("text", x=1466.75-1500, y=1.10, label= "Inferred conversion tract", colour = "grey50") + 
    annotate("segment", x = 1592.5-1500, y = 2.125, xend = 1685-1500, yend = 2.125, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    annotate("segment", x = 1592.5-1500, y = 2.125, xend = 1500-1500, yend = 2.125, colour = "grey50", arrow=arrow(length=unit(0.3,"cm"))) +
    annotate("text", x=1592.5-1500, y=2.25, label= "Inferred CT", colour = "grey50") + 
    
    scale_shape_discrete(guide=FALSE) +
    xlim(c(1200-1500,2003-1500)) + scale_colour_manual(values=c("#993333", "#FFCC00"), labels = c("B6", "CAST"), na.value = "grey50")


tracks("PRDM9 ChIP-seq" = p1,"Markers" = p2_bis, "Fragments" = p3_bis,
       heights=c(2,1,4))










# Plot --------------------------------------------------------------------


png(paste("./output/", fig_name_2, ".png", sep = ""), 
    width = width_png, height = height_png)


dev.off()




