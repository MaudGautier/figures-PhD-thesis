
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

