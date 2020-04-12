source("src/libraries.R")



col_names <- c("ID", 
               "pos_dom", "start_dom", "end_dom", "strand_dom", "score_dom", "pval_dom", "qval_dom", 
               "pos_cas", "start_cas", "end_cas", "strand_cas", "score_cas", "pval_cas", "qval_cas")
table_dom <- read.table("./data/P9_Dom.all_selected_hotspots.Best_score.joined_haplotypes.txt", 
                        col.names = col_names)
table_cas <- read.table("./data/P9_Cast.all_selected_hotspots.Best_score.joined_haplotypes.txt", 
                        col.names = col_names)
table_dom_tB <- table_dom[grep("tB", table_dom$ID),]
table_dom_tC <- table_dom[grep("tC", table_dom$ID),]
table_cas_tB <- table_cas[grep("tB", table_cas$ID),]
table_cas_tC <- table_cas[grep("tC", table_cas$ID),]

# 3. Create the vectors without missing values
d_dom <- table_dom$start_dom + (table_dom$end_dom - table_dom$start_dom)/2 - 500
d_dom <- d_dom[!is.na(d_dom)]
d_cas <- table_dom$start_cas + (table_dom$end_cas - table_dom$start_cas)/2 - 500
d_cas <- d_cas[!is.na(d_cas)]


d_dom <- table_cas$start_dom + (table_cas$end_dom - table_cas$start_dom)/2 - 500
d_dom <- d_dom[!is.na(d_dom)]
d_cas <- table_cas$start_cas + (table_cas$end_cas - table_cas$start_cas)/2 - 500
d_cas <- d_cas[!is.na(d_cas)]




### Separate

# table_dom_tB

d_dom <- table_dom_tB$start_dom + (table_dom_tB$end_dom - table_dom_tB$start_dom)/2 - 500
d_dom <- d_dom[!is.na(d_dom)]
d_cas <- table_dom_tB$start_cas + (table_dom_tB$end_cas - table_dom_tB$start_cas)/2 - 500
d_cas <- d_cas[!is.na(d_cas)]
dat <- data.frame(dens = c(d_dom, d_cas), 
                  Haplotype = c(rep(c("B6"), length(d_dom)), rep(c("CAST"), length(d_cas))))

ggplot(dat, aes(x = dens, col = Haplotype)) + 
    theme_classic() + 
    geom_density(aes(y = ..count..), size=1) +
    labs(x = "Position relative to the PRDM9 peak summit (bp)",
         y = expression(paste("Number of ", italic(Prdm9^Dom2), " motifs"))) +
    theme(text=element_text(family="LM Roman 10", size = 35),
          legend.position=c(0.88,0.90),
          legend.key.size = unit(2, 'lines')) + 
    ylim(c(0,12.5)) +
    scale_color_manual(name = "Haplotype", breaks = c("B6", "CAST"), values = c("#B95353", "#EFBC00"))


# table_cas_tB #155
d_dom <- table_cas_tB$start_dom + (table_cas_tB$end_dom - table_cas_tB$start_dom)/2 - 500
d_dom <- d_dom[!is.na(d_dom)]
d_cas <- table_cas_tB$start_cas + (table_cas_tB$end_cas - table_cas_tB$start_cas)/2 - 500
d_cas <- d_cas[!is.na(d_cas)]
dat <- data.frame(dens = c(d_dom, d_cas), 
                  Haplotype = c(rep(c("B6"), length(d_dom)), rep(c("CAST"), length(d_cas))))

ggplot(dat, aes(x = dens, col = Haplotype)) + 
    theme_classic() + 
    geom_density(aes(y = ..count..), size=1) +
    labs(x = "Position relative to the PRDM9 peak summit (bp)",
         y = expression(paste("Number of ", italic(Prdm9^Cst), " motifs"))) +
    theme(text=element_text(family="LM Roman 10", size = 35),
          legend.position=c(0.88,0.90),
          legend.key.size = unit(2, 'lines')) + 
    ylim(c(0,12.5)) +
    scale_color_manual(name = "Haplotype", breaks = c("B6", "CAST"), values = c("#B95353", "#EFBC00"))



# table_dom_tC #217
d_dom <- table_dom_tC$start_dom + (table_dom_tC$end_dom - table_dom_tC$start_dom)/2 - 500
d_dom <- d_dom[!is.na(d_dom)]
d_cas <- table_dom_tC$start_cas + (table_dom_tC$end_cas - table_dom_tC$start_cas)/2 - 500
d_cas <- d_cas[!is.na(d_cas)]
dat <- data.frame(dens = c(d_dom, d_cas), 
                  Haplotype = c(rep(c("B6"), length(d_dom)), rep(c("CAST"), length(d_cas))))

ggplot(dat, aes(x = dens, col = Haplotype)) + 
    theme_classic() + 
    geom_density(aes(y = ..count..), size=1) +
    labs(x = "Position relative to the PRDM9 peak summit (bp)",
         y = expression(paste("Number of ", italic(Prdm9^Dom2), " motifs"))) +
    theme(text=element_text(family="LM Roman 10", size = 35),
          legend.position=c(0.88,0.90),
          legend.key.size = unit(2, 'lines')) + 
    ylim(c(0,12.5)) +
    scale_color_manual(name = "Haplotype", breaks = c("B6", "CAST"), values = c("#B95353", "#EFBC00"))



# table_cas_tC #237
d_dom <- table_cas_tC$start_dom + (table_cas_tC$end_dom - table_cas_tC$start_dom)/2 - 500
d_dom <- d_dom[!is.na(d_dom)]
d_cas <- table_cas_tC$start_cas + (table_cas_tC$end_cas - table_cas_tC$start_cas)/2 - 500
d_cas <- d_cas[!is.na(d_cas)]
dat <- data.frame(dens = c(d_dom, d_cas), 
                  Haplotype = c(rep(c("B6"), length(d_dom)), rep(c("CAST"), length(d_cas))))

ggplot(dat, aes(x = dens, col = Haplotype)) + 
    theme_classic() + 
    geom_density(aes(y = ..count..), size=1) +
    labs(x = "Position relative to the PRDM9 peak summit (bp)",
         y = expression(paste("Number of ", italic(Prdm9^Cst), " motifs"))) +
    theme(text=element_text(family="LM Roman 10", size = 35),
          legend.position=c(0.88,0.90),
          legend.key.size = unit(2, 'lines')) + 
    ylim(c(0,12.5)) +
    scale_color_manual(name = "Haplotype", breaks = c("B6", "CAST"), values = c("#B95353", "#EFBC00"))



