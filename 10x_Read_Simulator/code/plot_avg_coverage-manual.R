#!/usr/bin/env Rscript
# R 3.3

# USAGE: plot_avg_coverage.R data/depth.averages.txt /path/to/outdir

# DESCRIPTION:
# This script will read in a text file generated with the 'average_coverages.py' script
# and output a plot showing the average coverage per chromosome per genome (column)

library("reshape2")
library("ggplot2")
library("plotly")

# get commands passed to the script
# get commands passed to the script
args <- commandArgs(TRUE)
cov_file <- args[1]
outdir <- args[2]
# cov_file <- "/home/devsci4/Structural_Variants_CSHL/10x_Read_Simulator/test_output/depth13Genome.depth.averages.txt"
# outdir <- "/home/devsci4/Structural_Variants_CSHL/10x_Read_Simulator/test_output/"
# > load("/home/devsci4/Structural_Variants_CSHL/10x_Read_Simulator/test_output/plot_avg-manual.Rdata")

# read in the file
coverage_df <- read.table(cov_file)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
# colnames(coverage_df)[-1] <- paste("genome_", seq_along(colnames(coverage_df)[-1]), sep="")

sample_names <- c("HG00512"
, "HG00513"
, "HG00514"
, "HG00731"
, "HG00732"
, "HG00733"
, "NA19238"
, "NA19239"
, "NA19240"
, "NA24143"
, "NA24149"
, "NA24385"
, "NA12878")
colnames(coverage_df)[-1] <- sample_names

# chrom_vec <- dput(scan("hg38_chrom_list.txt",what = ""))
# chrom_vec <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
# "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
# "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
# "chrX", "chrY")
# chrom_vec <- factor(chrom_vec, levels = chrom_vec)
# print(chrom_vec)
chrom_vec <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
multi_grep <- function(source_data, patterns){
    # find multiple patterns in a char vector
    # make a regex to search with
    if (length(patterns) > 1){
        search_pattern <- paste(patterns,
                                collapse = "|") 
    } else {
        search_pattern <- patterns
    }

    matches <- grep(x = source_data,
                    pattern = search_pattern,
                    value = TRUE)
    return(matches)
}

# get only the desired chrom's
new_df <- data.frame()
for(chrmatch in multi_grep(chrom_vec, coverage_df[["chrom"]])){
# print(chrmatch)
new_df <- rbind(new_df, subset(coverage_df, chrom == chrmatch))
}
coverage_df <- droplevels(new_df)
coverage_df[["chrom"]] <- factor(coverage_df[["chrom"]], levels = chrom_vec)

# coverage_df <- coverage_df[!grepl("_", coverage_df[["chrom"]]),]
# coverage_df <- coverage_df[!grepl("chrEBV", coverage_df[["chrom"]]),]
# coverage_df <- coverage_df[!grepl("hs38d1", coverage_df[["chrom"]]),]
# coverage_df <- coverage_df[!grepl("chrM", coverage_df[["chrom"]]),]


# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = "chrom", value.name = "coverage", variable.name = "sample")


# fix chrom order for plot
# coverage_df <- coverage_df[with(coverage_df, order(chrom)), ]

# peel off the coverage stats column and turn into a grouping factor
stat_strings <- strsplit(as.character(coverage_df$coverage), ',')

stats_df <- data.frame(matrix(as.numeric(unlist(stat_strings)), nrow=length(stat_strings), byrow=T))
colnames(stats_df) <- c("average", "std_dev", "count")

coverage_df <- cbind(coverage_df[c("chrom", "sample")], stats_df)

# print(head(coverage_df))
# quit()


# plot by genome
# coverage_df_avg <- subset(coverage_df, statistic == "average")
pdf(file = file.path(outdir, "avg_cov_byGenome-manual.pdf"), height = 8, width = 8)
chrom_plot <- ggplot(coverage_df, aes(x = sample, y = average, fill = factor(chrom)))
chrom_plot <-chrom_plot + geom_bar(stat="identity", position="dodge")
chrom_plot <-chrom_plot + coord_flip()
# chrom_plot <-chrom_plot + scale_x_discrete(limits = rev(levels(coverage_df[["chrom"]])))
chrom_plot <-chrom_plot + labs(title="Average Coverage Per Chromosome\nPer Samples", x="Sample", y = "Average Coverage", fill="Chromosome")
print(chrom_plot)
dev.off()


save.image(file=file.path(outdir, "plot_avg-manual.Rdata"),compress = TRUE)
# load("/home/devsci4/Structural_Variants_CSHL/10x_Read_Simulator/test_output/plot_avg-manual.Rdata")

# plotly
# chrom_plotly <- ggplotly(chrom_plot)
# htmlwidgets::saveWidget(as.widget(chrom_plotly), file.path(outdir, "avg_cov_byGenome-manual.html"))


# coverage_df
# make horizontal stacked grouped barplot
# plot by chrom
# pdf(file = file.path(outdir, "avg_cov_byChrom.pdf"), height = 8, width = 8)
# ggplot(coverage_df, aes(x = chrom, y = coverage, fill = factor(sample))) +
#   geom_bar(stat="identity", position="dodge") + # remove 'position' for stacked plot
#     coord_flip() + 
#     labs(title="Average Coverage Per Chromosome\nPer Samples", x="Chromosome", y = "Average Coverage")
# dev.off()
