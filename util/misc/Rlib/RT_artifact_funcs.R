#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse"))

## functions:

# stingray plot

local_lib_dir=dirname(sys.frame(1)$ofile) #https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script

stingray_plot = function(data, include_background=FALSE) {

    stringray_boilerplate_plot = ggplot()

    if (include_background) {
        ## build stingray plot template.


        background_rds_file = file.path(local_lib_dir, "../data/background_ffpm_far.rds")

        background_ffpm_far = readRDS(background_rds_file) # based on GTEx
        background_ffpm_far = background_ffpm_far %>% mutate(log2FFPM = log2(FFPM + 1),
                                                             log2_FAR_left = log2(FAR_left+1),
                                                             log2_FAR_right = -1*log2(FAR_right+1))  %>%
            gather(key=FAR_type, value=log2_FAR, log2_FAR_left, log2_FAR_right)


        stringray_boilerplate_plot = stringray_boilerplate_plot  + geom_point(data=background_ffpm_far, aes(y=log2_FAR, x=log2FFPM), color='lightgray', alpha=0.2)
    }

    ## now plot the given data on top of it.
    plotdata = data %>% mutate(log2FFPM = log2(FFPM + 1),
                             log2_FAR_left = log2(FAR_left+1),
                             log2_FAR_right = -1*log2(FAR_right+1)) %>%
                      gather(key=FAR_type, value=log2_FAR, log2_FAR_left, log2_FAR_right)


    ymax = max(abs(plotdata$log2_FAR), abs(background_ffpm_far$log2_FAR))
    p = stringray_boilerplate_plot + geom_point(data=plotdata, aes(y=log2_FAR, x=log2FFPM, color=FAR_type)) + ylim(-1*ymax, ymax)

    return(p)
}




breakpoint_plot = function(fusion_data, microhomology_data, title, fusion_brkpt_size_by=c("FFPM", "num_samples")) {

    #message(fusion_data)
    
    fusion_brkpt_size_by = match.arg(fusion_brkpt_size_by)

    splice_type_colors = c("ONLY_REF_SPLICE" = "purple",
                           "INCL_NON_REF_SPLICE" = "orange",
                           "NO_JUNCTION_READS_IDENTIFIED" = "red")

    splice_dinuc_shapes = c("Consensus" = 16,
                            "Non" = 17)
    
    
    if (fusion_brkpt_size_by == "FFPM") {

        p = fusion_data %>% ggplot() +
            geom_point(aes(x=LeftLocalBreakpoint, y=RightLocalBreakpoint, size=FFPM,
                           color=SpliceType, shape=SpliceDinuc), alpha=0.5) +
                           scale_color_manual(values=splice_type_colors) +
                           scale_shape_manual(values=splice_dinuc_shapes)

    } else if (fusion_brkpt_size_by == "num_samples") {

        p = fusion_data %>% mutate(brkpt = paste(LeftBreakpoint, RightBreakpoint, sep="::")) %>%
                            group_by(brkpt) %>% mutate(num_samples = n()) %>% ungroup() %>%
                            ggplot() +
                            geom_point(aes(x=LeftLocalBreakpoint, y=RightLocalBreakpoint, size=num_samples,
                                           color=SpliceType, shape=SpliceDinuc), alpha=0.3) +
                            scale_color_manual(values=splice_type_colors) +
                            scale_shape_manual(values=splice_dinuc_shapes)

    } else {
        ## shouldn't happen.
        stop("Error, not recognizing fusion_brkpt_size_by_param: ", fusion_brkpt_size_by)
    }


    p = p + ggtitle(title)

    ## plot gene structures.
    ## geneA on x-axis, geneB on y-axis

    geneA_info = microhomology_data %>% filter(feature_type == "GeneA")
    geneB_info = microhomology_data %>% filter(feature_type == "GeneB")

    ## plot gene structures.
    padding = 1000

    geneA_lend = min(geneA_info$lend)
    geneA_rend = max(geneA_info$rend)
    geneA_length = geneA_rend - geneA_lend + 2*padding

    geneB_lend = min(geneB_info$lend)
    geneB_rend = max(geneB_info$rend)
    geneB_length = geneB_rend - geneB_lend + 2*padding

    ## set up scales for view
    p = p +
        xlim(geneA_lend - padding, geneA_rend + padding) +
        ylim(geneB_lend - padding, geneB_rend + padding)


    ## draw geneA
    geneA_minY = geneB_lend - 0.95*padding
    geneA_maxY = geneA_minY + 0.02*geneB_length

    p = p + geom_rect(data=geneA_info, aes(xmin=lend, xmax=rend, ymin=geneA_minY, ymax=geneA_maxY), fill=NA, color='black', alpha=0.5)

    geneA_midY = mean(c(geneA_minY, geneA_maxY))

    p = p + geom_segment(aes(x=geneA_lend, xend=geneA_rend, y=geneA_midY, yend=geneA_midY), alpha=0.5) # geneB center line

    ## draw geneB

    geneB_minX = geneA_lend - 0.95*padding
    geneB_maxX = geneB_minX + 0.03*geneA_length

    p = p + geom_rect(data=geneB_info, aes(ymin=lend, ymax=rend, xmin=geneB_minX, xmax=geneB_maxX), fill=NA, color='black', alpha=0.5)

    geneB_midX = mean(c(geneB_minX, geneB_maxX))
    p = p + geom_segment(aes(x=geneB_midX, xend=geneB_midX, y=geneB_lend, yend=geneB_rend), alpha=0.5) # geneB center line

    microhomology_info = microhomology_data %>% filter(feature_type == "MicroH")

    if (nrow(microhomology_data) > 0 ) {
        p = p + geom_point(data=microhomology_info, aes(x=lend, y=rend), size=0.5)
    }

    return(p)
}


