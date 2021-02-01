#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ranger"))


parser = ArgumentParser()
parser$add_argument("--fusions", help="input fusions data file with all attributes assigned", required=TRUE, nargs=1)
parser$add_argument("--output", help="output filename", required=TRUE, nargs=1)
parser$add_argument("--ranger", help="ranger predictor rds obj", required=TRUE, nargs=1)

args = parser$parse_args()

dat_filename = args$fusions
out_filename = args$output
rg_rds_file = args$ranger

message("-parsing ", dat_filename)
data = read.table(dat_filename, header=T, sep="\t", stringsAsFactors = F)

orig_data = data


data = data %>% select(annot_splice, consensus_splice, FFPM, left_counter_ffpm, right_counter_ffpm, FAR_left, FAR_right, microh_brkpt_dist, num_microh)

pseudocount = 1

data$adj_FFPM=log2(data$FFPM + pseudocount)
data$adj_left_counter_ffpm = log2(data$left_counter_ffpm + pseudocount)
data$adj_right_counter_ffpm = log2(data$right_counter_ffpm + pseudocount)
data$adj_FAR_left = log2(data$FAR_left + pseudocount)
data$adj_FAR_right = log2(data$FAR_right + pseudocount)
data$adj_microh_brkpt_dist = log2(data$microh_brkpt_dist + pseudocount)
data$adj_num_microh = log2(data$num_microh + pseudocount)
data$adj_annot_splice = data$annot_splice
data$adj_consensus_splice = data$consensus_splice


## center and scale

##> data.colmeans
##      annot_splice   consensus_splice               FFPM  left_counter_ffpm right_counter_ffpm           FAR_left
##         0.8144065          0.9149324          0.1815195          1.3516459          1.3662279          1.2690494
##         FAR_right  microh_brkpt_dist         num_microh
##         1.3892456          9.8476273          0.2846550


data.scaled = data
data.scaled = data.scaled %>% mutate(adj_annot_splice = adj_annot_splice - 0.8144065,
                                     adj_consensus_splice = adj_consensus_splice - 0.9149324,
                                     adj_FFPM = adj_FFPM - 0.1815195,
                                     adj_left_counter_ffpm = adj_left_counter_ffpm - 1.3516459,
                                     adj_right_counter_ffpm = adj_right_counter_ffpm - 1.3662279,
                                     adj_FAR_left = adj_FAR_left - 1.2690494,
                                     adj_FAR_right = adj_FAR_right - 1.3892456,
                                     adj_microh_brkpt_dist = adj_microh_brkpt_dist - 9.8476273,
                                     adj_num_microh = adj_num_microh - 0.2846550)


##> data.sd
##      annot_splice   consensus_splice               FFPM  left_counter_ffpm right_counter_ffpm           FAR_left
 ##        0.3887820          0.2789849          0.3744377          2.4123070          2.6002049          1.2838213
##         FAR_right  microh_brkpt_dist         num_microh
##         1.2570087          2.4569856          0.9846331


data.scaled = data.scaled %>% mutate(adj_annot_splice = adj_annot_splice / 0.3887820,
                                     adj_consensus_splice = adj_consensus_splice / 0.2789849,
                                     adj_FFPM = adj_FFPM  / 0.3744377,
                                     adj_left_counter_ffpm = adj_left_counter_ffpm / 2.4123070,
                                     adj_right_counter_ffpm = adj_right_counter_ffpm / 2.6002049,
                                     adj_FAR_left = adj_FAR_left / 1.2838213,
                                     adj_FAR_right = adj_FAR_right / 1.2570087,
                                     adj_microh_brkpt_dist = adj_microh_brkpt_dist / 2.4569856,
                                     adj_num_microh = adj_num_microh / 0.9846331)


## trim outliers, set to range [-2,2]
data.scaled[data.scaled < -2] = -2
data.scaled[data.scaled > 2] = 2

## rescale so all values are between -2,2

## apply(data.scaled, 2, range)
##     annot_splice consensus_splice       FFPM left_counter_ffpm right_counter_ffpm   FAR_left FAR_right
##[1,]   -2.0000000       -2.0000000 -0.4747742        -0.5603125         -0.5254309 -0.9884938   -1.1052
##[2,]    0.4773718        0.3049183  2.0000000         2.0000000          2.0000000  2.0000000    2.0000
##     microh_brkpt_dist num_microh
##[1,]         -2.000000 -0.2890975
##[2,]          1.400183  2.0000000


scale_range_min = -2
scale_range_max = 2
scale_range_size = scale_range_max - scale_range_min

data.scaled = data.scaled %>% mutate(adj_annot_splice = (adj_annot_splice - -2) / 2.4773718 * scale_range_size + scale_range_min,
                                     adj_consensus_splice = (adj_consensus_splice - -2) / 2.3049183 * scale_range_size + scale_range_min,
                                     adj_FFPM = (adj_FFPM - -0.4747742) / 2.4747742 * scale_range_size + scale_range_min,
                                     adj_left_counter_ffpm = (adj_left_counter_ffpm - -0.5603125) / 2.5603125 * scale_range_size + scale_range_min,
                                     adj_right_counter_ffpm = (adj_right_counter_ffpm - -0.5254309) / 2.5254309 * scale_range_size + scale_range_min,
                                     adj_FAR_left = (adj_FAR_left - -0.9884938) / 3.9884938 * scale_range_size + scale_range_min,
                                     adj_FAR_right = (adj_FAR_right - -1.1052) / 3.1052 * scale_range_size + scale_range_min,
                                     adj_microh_brkpt_dist = (adj_microh_brkpt_dist - -2) / 3.400183 * scale_range_size + scale_range_min,
                                     adj_num_microh = (adj_num_microh - -0.2890975) / 2.2890975 * scale_range_size + scale_range_min)


## build predictor:
message("-building predictor")
#umap.merged.layout = read.table(train_data_filename, header=T, stringsAsFactors = F, sep="\t")


#rangerdata = umap.merged.layout %>% select(leiden, FFPM, annot_splice, consensus_splice,
#                       left_counter_ffpm, right_counter_ffpm,  FAR_left, FAR_right,
#                       microh_brkpt_dist, num_microh)
#
#rangerdata$leiden = factor(rangerdata$leiden)
# rg = ranger(leiden ~ ., data=rangerdata)


rg = readRDS(rg_rds_file)


message("-predicting fusion clusters")
pred = predict(rg, data=data.scaled)

orig_data$pred_cluster = pred$predictions

## annotate clusters according to attribute types.
orig_data = orig_data %>%
    mutate(fusion_cluster_att = ifelse(leiden == 4, "cosmic-like", "NA")) %>%
    mutate(fusion_cluster_att = ifelse(leiden %in% c(47,20,41,15), "expr_microH_RT_artifact?", fusion_cluster_att)) %>%
    mutate(fusion_cluster_att = ifelse(leiden %in% c(57,56,60), "high_FAR_microH_bioinf_artifact?", fusion_cluster_att)) %>%
    mutate(fusion_cluster_att = ifelse(leiden %in% c(49,51), "high_counter_evidence", fusion_cluster_att))


write.table(orig_data, file=out_filename, quote=F, sep="\t", row.names=F)

