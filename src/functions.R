
formatter1000000 <- function(x){ 
    x/1000000 
}
formatter100 <- function(x){ 
    x/100 
}


std <- function(x) sd(x)/sqrt(length(x))
# sd_upper=function(x) { mean(x)+sd(x) }
# sd_lower=function(x) { mean(x)-sd(x) }
# std_upper=function(x) { mean(x)+std(x) }
# std_lower=function(x) { mean(x)-std(x) }


def_sizes <- function(x){ifelse(x == "SNP", 2,2)}


# Get right colours
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


# Rename categories of hotspots
label_categ <- function(x) {
    if (x=="tC.hB+C") {return("tC.sym")}
    if (x=="tC.hB")   {return("tC.chC")}
    if (x=="tC.hC")   {return("tC.chB")}
    if (x=="tB.hB")   {return("tB.chC")}
    if (x=="tB.hB+C") {return("tB.sym")}
    if (x=="tB.hC")   {return("tB.chB")}
    if (x=="NOV.tB.hC") {return("NOV.tB.chB")}
    if (x=="NOV.tC.hB") {return("NOV.tC.chC")}
    if (x=="NOV.NA")  {return("NOV.sym")}
    if (x=="INDEP.hB") {return("INDEP.chC")}
    if (x=="INDEP.hB+C") {return("INDEP.sym")}
    if (x=="INDEP.hC") {return("INDEP.chB")}
}
status_labeller <- function(variable,value){
    status_names <- list(
        'CO'="Rec-1S",
        'Complex'="Rec-MS",
        'NCO'="Rec-2S"
    )
    return(status_names[value])
}



label_hot_categ_NOV_sep <- function(x) {
    if (is.na(x)) {return(NA)}
    if (x=="tC.hB+C") {return("CAST")}
    if (x=="tC.hB")   {return("CAST")}
    if (x=="tC.hC")   {return("CAST")}
    if (x=="tB.hB")   {return("B6")}
    if (x=="tB.hB+C") {return("B6")}
    if (x=="tB.hC")   {return("B6")}
    if (x=="NOV.tB.hC") {return("Novel_B6")}
    if (x=="NOV.tC.hB") {return("Novel_CAST")}
    if (x=="NOV.NA")  {return(NA)}
    if (x=="INDEP.hB") {return(NA)}
    if (x=="INDEP.hB+C") {return(NA)}
    if (x=="INDEP.hC") {return(NA)}
    
}


plot.multi.dens <- function(s, title = "")
{
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s)) {
        junk.x = c(junk.x, density(s[[i]])$x)
        junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = title)
    for(i in 1:length(s)) {
        lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
    }
}

