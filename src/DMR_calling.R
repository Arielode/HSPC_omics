#!/usr/bin/env Rscript

### DMR_calling.R
###
### Copyright(C)Chao Tang
### Contact: Chao Tang <tangchao@ihcams.ac.cn>
###
### Identify DMRs between two specified groups.
###
### Usage:
###   DMR_calling.R -f1 <input1> -f2 -t <input2> -d <methylation differences cutoff> -q <adjusted p value cutoff>
###
### Options:
###   -f1    group1 methylation data split by tiles.(can obtained by Tile.R)
###   -f2    group2 methylation data split by tiles.
###   -d    Methylation differences cutoff between compared groups.
###   -q    Adjusted p value cutoff.
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
                    make_option(c("-f1", "--file1"), dest = "file1", default = "",
                                help = "methylation profile of group1"),
                    make_option(c("-f2", "--file2"), dest = "file2", default = "", 
                                help = "methylation profile of group2"), 
                    make_option(c("-d", "--diff"), dest = "diff", default = 0.2,
                                help = "The minimal number of CpG sites in a tile"), 
                    make_option(c("-q", "--qvalue"), dest = "qvalue", default = 0.05, 
                                help = "The adjusted p value cutoff")
)


parser <- OptionParser(usage = "DMR_calling.R [options]",
                       option_list=option_list, description = "Description: Identify DMRs between two specified groups.\
                       Contact: Chao Tang <tangchao@ihcams.ac.cn>.\
"
)

## check arguments
arguments <- parse_args(parser)
file1 <- arguments$file1
if(file1 == ""){ # default, STDIN
  print("please specify the input file")
} else { # user specified
  if( file.access(file1) == -1 ){ # file not exists
    print_help(parser)
  }
}

file2 <- arguments$file2
if(file2 == ""){ # default, STDIN
  print("please specify the input file")
} else { # user specified
  if( file.access(file2) == -1 ){ # file not exists
    print_help(parser)
  }
}

diff <- arguments$diff
if(diff == ""){
    diff <- 0.2
}

qvalue <- arguments$qvalue
if(qvalue == ""){
    qvalue <- 0.05
}

txt1 = fread(file1, header=T)
txt2 = fread(file2, header=T)

txt1 = txt1 %>% group_by(Feature) %>% summarise(total_reads=sum(total_reads), met_reads=sum(met_reads))
txt2 = txt2 %>% group_by(Feature) %>% summarise(total_reads=sum(total_reads), met_reads=sum(met_reads))
inter = data.frame(Feature=intersect(txt1$Feature, txt2$Feature))

comp1_test = left_join(inter, txt1)
comp2_test = left_join(inter, txt2)
comp1_test$unmet_reads = comp1_test$total_reads - comp1_test$met_reads
comp2_test$unmet_reads = comp2_test$total_reads - comp2_test$met_reads

df=data.frame(pos=character(nrow(comp1_test)), pvalue=numeric(nrow(comp1_test)), met1=numeric(nrow(comp1_test)), met2=numeric(nrow(comp1_test)))

for(i in 1:nrow(comp1_test)){
    df$pos[i]=comp1_test$Feature[i]
    mat = matrix(c(comp1_test$met_reads[i], comp2_test$met_reads[i], comp1_test$unmet_reads[i], comp2_test$unmet_reads[i]), nrow=2, dimnames=list(comp=c("comp1", "comp2"), statu=c("met", "unmet")))
    df$pvalue[i]=tryCatch({ fisher.test(mat, alternative = "two.sided")$p.value }, error=function(e) NA)
	df$met1[i] = comp1_test$met_reads[i] / comp1_test$total_reads[i]
	df$met2[i] = comp2_test$met_reads[i] / comp2_test$total_reads[i]
}
q=p.adjust(df$pvalue, method="BH")
df$qvalue=q

df$methydiff <- df$met2 - df$met1
df <- df[df$qvalue <- qvalue & abs(df$methydiff) > diff, ]


write.table(df, file=paste0(file1, "_vs_", file2, "_DMR.txt"), sep="\t", quote=F, row.names=F)

