#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ranger"))
suppressPackageStartupMessages(library("tidyverse"))


message("-building rg obj")
traindata = read.table(gzfile("ranger.test_data.gz"), header=T, stringsAsFactors=F, sep="\t")
rg = ranger(leiden ~ ., data=traindata)

message("-saving rg obj")
saveRDS(rg, file="ranger.rg_obj.rds")

message("-done")
quit(save = "no", status = 0, runLast = FALSE)

