#!/usr/bin/env Rscript
# R 3.3

# USAGE: plot_bin_coverages.R data/depth.bin_coverages.txt /path/to/outdir

# DESCRIPTION:
# This script will plot the binned coverage data

library("reshape2")
library("ggplot2")
theme_set(theme_bw())
library("plotly")

# get commands passed to the script
args <- commandArgs(TRUE)
cov_file <- args[1]
outdir <- args[2]

# read in the file
coverage_df <- read.table(cov_file)
save.image(file=file.path(outdir, "plot_bin_start.Rdata"),compress = TRUE)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
colnames(coverage_df)[2] <- "start"
colnames(coverage_df)[3] <- "stop"
colnames(coverage_df)[-(1:3)] <- paste("genome_", seq_along(colnames(coverage_df)[-(1:3)]), sep="")

coverage_df['region'] <- paste(coverage_df[[2]], coverage_df[[3]], sep = '-')
coverage_df['chrom_region'] <- paste(coverage_df[[1]], coverage_df[["region"]], sep = ':')

# get bin size from first row
bin_size <- max(coverage_df[["stop"]] - coverage_df[["start"]])
# drop start stop columns
coverage_df <- coverage_df[-(2:3)]




# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = c("chrom", "region", 'chrom_region'), value.name = "total_coverage", variable.name = "sample")

# fix chrom order for plot
# coverage_df <- coverage_df[with(coverage_df, order(chrom, region, sample)), ]

# peel off the coverage stats column and turn into a grouping factor
stat_strings <- strsplit(as.character(coverage_df$total_coverage), ',')

stats_df <- data.frame(matrix(as.numeric(unlist(stat_strings)), nrow=length(stat_strings), byrow=T))
colnames(stats_df) <- c("average", "std_dev", "count")

# put it back together
coverage_df <- cbind(coverage_df[c("chrom", "region", 'chrom_region', "sample")], stats_df)




pdf(file = file.path(outdir, "cov_byRegion_density.pdf"), height = 8, width = 8, onefile = TRUE)
for(i in seq_along(levels(coverage_df[['chrom']]))){
    ichrom <- levels(coverage_df[['chrom']])[i]
    cov_subset <- subset(coverage_df, chrom == ichrom & count == bin_size)
    # print(head(cov_subset))
    myplot <- ggplot(cov_subset, aes(x = average))
    myplot <- myplot + geom_density(alpha=.5, fill="#FF6666") # geom_histogram(alpha=.5, position="identity") + 
    myplot <- myplot + labs(title=paste0("Chromosome: ", ichrom, "\nAverage Coverage"), x="Average Coverage", y = "Density")
    print(myplot)

}
dev.off()

# ribbon plot
pdf(file = file.path(outdir, "cov_byRegion_line.pdf"), height = 10, width = 10, onefile = TRUE)
for(i in seq_along(levels(coverage_df[['chrom']]))){
    ichrom <- levels(coverage_df[['chrom']])[i]
    cov_subset <- subset(coverage_df, chrom == ichrom & count == bin_size)
    head(cov_subset)
    cov_subset[["region"]] <- as.factor(cov_subset[["region"]])
    cov_subset[["chrom_region"]] <- as.factor(cov_subset[["chrom_region"]])
    # print(head(cov_subset))
    chrom_lines <- ggplot(cov_subset, aes(x=region, y=average, group=sample, colour = sample))
    chrom_lines <- chrom_lines + geom_line() # size = 1
    # chrom_lines <- chrom_lines + geom_ribbon(aes(x=region, ymin=max(0, average - std_dev), ymax=average + std_dev, fill = sample), alpha = 0.05)
    # chrom_lines <- chrom_lines + geom_errorbar(aes(ymin=average - std_dev, ymax=average + std_dev), width=0)
    chrom_lines <- chrom_lines + labs(title=paste0("Chromosome: ", ichrom, "\nAverage Coverage"), x="Region", y = "Average Coverage")
    chrom_lines <- chrom_lines + theme(axis.text.x  = element_text(angle=90, size=6))
    print(chrom_lines)
    
}
dev.off()

# plotly
# chrom_plotly <- ggplotly(chrom_plot)
# htmlwidgets::saveWidget(as.widget(chrom_plotly), file.path(outdir, "avg_cov_byGenome.html"))


# pdf(file = file.path(outdir, "total_cov_byRegion_hist.pdf"), height = 8, width = 8, onefile = TRUE)
# for(i in seq_along(levels(coverage_df[['chrom']]))){
#     ichrom <- levels(coverage_df[['chrom']])[i]
#     cov_subset <- subset(coverage_df, chrom == ichrom)
#     # print(head(cov_subset))
#     myplot <- ggplot(cov_subset, aes(x = average, fill = factor(sample)))
#     myplot <- myplot + geom_histogram(alpha=.5, position="identity")
#     myplot <- myplot + labs(title=paste0("Chromosome: ", ichrom, "\nTotal Coverage"), x="Total Coverage", y = "Count")
#     print(myplot)
# 
# }
# dev.off()

# pdf(file = file.path(outdir, "total_cov_byRegion.pdf"), height = 8, width = 8)
# dev.off()