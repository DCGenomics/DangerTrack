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
args <- commandArgs(TRUE)
cov_file <- args[1]
outdir <- args[2]

# read in the file
coverage_df <- read.table(cov_file)
# save.image(file=file.path(outdir, "plot_avg_start.Rdata"),compress = TRUE)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
colnames(coverage_df)[-1] <- paste("genome_", seq_along(colnames(coverage_df)[-1]), sep="")

# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = "chrom", value.name = "coverage", variable.name = "sample")

# peel off the coverage stats column and turn into a grouping factor
stat_strings <- strsplit(as.character(coverage_df$coverage), ',')

stats_df <- data.frame(matrix(as.numeric(unlist(stat_strings)), nrow=length(stat_strings), byrow=T))
colnames(stats_df) <- c("average", "std_dev", "count")

coverage_df <- cbind(coverage_df[c("chrom", "sample")], stats_df)


# plot by genome
# horizontal barplot per genome
# coverage_df_avg <- subset(coverage_df, statistic == "average")
pdf(file = file.path(outdir, "avg_cov_byGenome.pdf"), height = 8, width = 8)
chrom_plot <- ggplot(coverage_df, aes(x = sample, y = average, fill = factor(chrom)))
chrom_plot <-chrom_plot + geom_bar(stat="identity", position="dodge")
chrom_plot <-chrom_plot + coord_flip()
chrom_plot <-chrom_plot + labs(title="Average Coverage Per Chromosome\nPer Samples", x="Sample", y = "Average Coverage", fill="Chromosome")
print(chrom_plot)
dev.off()


# plotly
chrom_plotly <- ggplotly(chrom_plot)
htmlwidgets::saveWidget(as.widget(chrom_plotly), file.path(outdir, "avg_cov_byGenome.html"))

quit()

# ~~~~~~~~~~~~ # # ~~~~~~~~~~~~ # # ~~~~~~~~~~~~ # 
# ~~~~~~~~~~~~ # # ~~~~~~~~~~~~ # # ~~~~~~~~~~~~ # 
# plotting methods still in development; line + ribbon, p
# ribbon plot
theme_set(theme_bw())
chrom_ribbon <- ggplot(coverage_df, aes(x=chrom, y=average, group=sample, colour = sample))
chrom_ribbon <- chrom_ribbon + geom_line()
chrom_ribbon <- chrom_ribbon + geom_ribbon(aes(x=chrom, ymin=max(0, average - std_dev), ymax=average + std_dev, fill = sample), alpha = 0.05)
print(chrom_ribbon)

# line + errorbars
chrom_lines <- ggplot(coverage_df, aes(x=chrom, y=average, group=sample, colour = sample))
chrom_lines <- chrom_lines + geom_errorbar(aes(ymin=average - std_dev, ymax=average + std_dev), width=.1)
chrom_lines <- chrom_lines + geom_line()
chrom_lines <- chrom_lines + geom_point()
print(chrom_lines)
# save.image(file=file.path(outdir, "plot_avg.Rdata"),compress = TRUE)
# load("plot_avg.Rdata")
print(coverage_sample_df)
quit()



# coverage_df
# make horizontal stacked grouped barplot
# plot by chrom
# pdf(file = file.path(outdir, "avg_cov_byChrom.pdf"), height = 8, width = 8)
# ggplot(coverage_df, aes(x = chrom, y = coverage, fill = factor(sample))) +
#   geom_bar(stat="identity", position="dodge") + # remove 'position' for stacked plot
#     coord_flip() + 
#     labs(title="Average Coverage Per Chromosome\nPer Samples", x="Chromosome", y = "Average Coverage")
# dev.off()
