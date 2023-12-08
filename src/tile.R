#!/usr/bin/env Rscript

### Tile.R
###
### Copyright(C)Chao Tang
### Contact: Chao Tang <tangchao@ihcams.ac.cn>
###
### Estimate the methylation level of samples in the resolution of tiles.
###
### Usage:
###   tile.R -i <input> -t <tile size> -n <number of CpG sites>
###
### Options:
###   -i    Specify the input folder containing the CG.txt files.
###   -t    Tile size to divide the genome.
###   -s    Number of CpG sites in a tile. 
###   -c    Number of cells covered the same tile.
###   -h    Show this message.


is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("optparse")){
    warning("Detect package \"optparse\" is not installed in your R enviroment.")
    warning("Trying to install the \"optparse\" package.")
    warning("If failed, please try to install it mannually.")
    install.packages("optparse")
}

## libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
                    make_option(c("-i", "--infile"), dest = "infile", default = "",
                                help = "input folder containing .CG.txt"),
                    make_option(c("-t", "--tile"), dest = "tile", default = 300, 
                                help = "tile size to divide genome"), 
                    make_option(c("-s", "--site"), dest = "site", default = 3,
                                help = "The minimal number of CpG sites in a tile"), 
                    make_option(c("-c", "--cell"), dest = "cell", default = 5, 
                                help = "The minimal number of cells covered a tile")
)


parser <- OptionParser(usage = "Tile.R [options]",
                       option_list=option_list, description = "Description: Estimate the methylation level of samples in the resolution of tiles.\
                       Contact: Chao Tang <tangchao@ihcams.ac.cn>.\
"
)

## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
if(infile == ""){ # default, STDIN
  infile <- setwd("./")
} else { # user specified
  if( file.access(infile) == -1 ){ # file not exists
    print_help(parser)
  }
}

tile <- arguments$tile
if(tile == ""){
    tile <- 300
}

site <- arguments$site
if(site == ""){
    site <- 3
}

cell <- arguments$cell
if(cell == ""){
    cell <- 5
}


f <- list.files(infile)
lst <- vector("list", length(f))
for(i in 1:length(lst)){
    w <- fread(file=paste0(infile, "/", f[i]))
    w$seg <- cut(w$V3, breaks = seq(0, to=((max(w$V3/tile)+1)*tile), by=tile))
    w <- unite(w, "Feature", sep="-", V1, seg)
    w2 <- w %>% group_by(Feature) %>% summarise(total_reads = sum(V8), met_reads = sum(V7), n=n())
    w2$Cell = str_sub(f[i], 1, -9)
    w2 <- w2[which(w2$n>=site), ]
    w2 <- w2[, c(1, 5, 2, 3)]
    lst[[i]] <- w2
}

## filter features by minimal covered cells
lst_cp = lst 
for(i in 1:length(lst_cp)){
        lst_cp[[i]] = lst_cp[[i]][, 1]

}
lst_cp = dplyr::bind_rows(lst_cp)
lst_cp = lst_cp %>% group_by(Feature) %>% summarise(n=n())
lst_cp = lst_cp[which(lst_cp$n>=cell), ]
feature = lst_cp$Feature

feature = data.frame(Feature = feature)
lst = dplyr::bind_rows(lst)
lst_final = left_join(feature, lst)

write.table(lst_final, file=paste0(infile, "/tile.txt"), sep="\t", quote=F, row.names=F)


