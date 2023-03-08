#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(dplyr)
library(tidyverse)


tmp <- read.delim(paste(args[1],".csv",sep=""),sep = ',')
tmp$site <- args[1]
tmp <- tmp[1:1000*floor(nrow(tmp)/1000),]
tmp[1,'Ucounts'] <- 0

tmp <- tmp %>%
  group_by(gr=gl(n()/1000,1000)) %>%
  mutate(Ucounts_1kb=sum(Ucounts),Bcounts_1kb=sum(Bcounts)) %>%
  slice(1000)


