#!/usr/bin/env Rscript

package_list<-c("ggplot2", "reshape2", "curl", "httr", "plotly")
# , "grid", "plyr", "knitr", "VennDiagram", "gridExtra", "datasets", "digest", "Hmisc", "xtable", "reshape2", "data.table", "scales", "corrplot", "RColorBrewer", "lattice", "gplots", "MASS", "stringr", "flsa", "genlasso", "optparse", "pastecs", "plotrix", "zoo", "reshape", "chron","UpSetR", "plotly") 

for(p in package_list){
        # check if the package is already installed
        if(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)){
            write(paste0("Package already installed: ", p), stderr())
        }
        
        # check if package can't be loaded
        if(!require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)){
            write(paste0("Attempting to install package: ",p), stderr())
            # try to install & load the packages, give a message upon failure
            tryCatch(install.packages(p,repos="http://cran.rstudio.com/"),
                     warning = function(e){write(paste0("Failed to install pacakge: ", p), stderr())},
                     error = function(e){write(paste0("Failed to install pacakge: ", p), stderr())})
            tryCatch(library(p,character.only=TRUE,verbose=FALSE),
                     warning = function(e){write(paste0("Package not installed: ", p), stderr())},
                     error = function(e){write(paste0("Package not installed: ", p), stderr())})
            
            # try to install & load the packages, skip to next loop iteration upon failure
            tryCatch(install.packages(p,repos="http://cran.rstudio.com/"),warning = next)
            tryCatch(library(p,character.only=TRUE,verbose=FALSE),warning = next)
        }
    }
