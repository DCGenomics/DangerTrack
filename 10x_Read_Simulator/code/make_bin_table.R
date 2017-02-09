#!/usr/bin/env Rscript
# R 3.3


library("RColorBrewer")
library("reshape2")
library("ggplot2")
library("gplots")
# library("plotly")


cov_file <- "test_output/depth13Genome.depth_bin-regions_chr21-1Mbp.txt"
outdir <- "test_output"

# read in the file
coverage_df <- read.table(cov_file)
save.image(file=file.path(outdir, "bin_table-manual.Rdata"),compress = TRUE)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
colnames(coverage_df)[2] <- "start"
colnames(coverage_df)[3] <- "stop"

# add sample names as colnames
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
colnames(coverage_df)[-(1:3)]  <- sample_names

# make a df of just avg coverage values
avg_cov_df <- coverage_df[-(1:3)]
avg_cov_df <- apply(X = avg_cov_df, MARGIN = c(1,2), FUN = function(x) as.numeric(gsub(x = as.character(x), pattern = "^(.*),.*,.*$",replacement = "\\1", perl = TRUE)) )

# make chrom region factor levels
coverage_df['region'] <- paste(coverage_df[[2]], coverage_df[[3]], sep = '-')
coverage_df['chrom_region'] <- paste(coverage_df[[1]], coverage_df[["region"]], sep = ':')

# get bin size from first row
bin_size <- max(coverage_df[["stop"]] - coverage_df[["start"]])

# set the avg cov rownames as region labels
row.names(avg_cov_df) <- coverage_df[["region"]]

# write out the avg cov table
write.csv(x = avg_cov_df, file = file.path(outdir, "avg_cov_table_chr21-1Mbp.csv"))

# reorder factor levels
# mydf$task <- factor(mydf$task, levels = c("up", "down", "left", "right", "front", "back"))
coverage_df[["chrom_region"]] <- factor(coverage_df[["chrom_region"]], levels = coverage_df[["chrom_region"]])
coverage_df[["region"]] <- factor(coverage_df[["region"]], levels = coverage_df[["region"]])

# print(levels(coverage_df[["chrom_region"]]))
# print(levels(coverage_df[["region"]]))
# quit()

# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = c("chrom", "region", 'chrom_region', "start", "stop"), value.name = "total_coverage", variable.name = "sample")

# peel off the coverage stats column and turn into a grouping factor
stat_strings <- strsplit(as.character(coverage_df$total_coverage), ',')
stats_df <- data.frame(matrix(as.numeric(unlist(stat_strings)), nrow=length(stat_strings), byrow=T))
colnames(stats_df) <- c("average", "std_dev", "count")

# put it back together
coverage_df <- cbind(coverage_df[c("chrom", "region", 'chrom_region', "sample", "start", "stop")], stats_df)

# make wide table for printing
# dcast(aql, month ~ variable)
# colnames(reshape2::dcast(coverage_df, average*std_dev*count ~sample))
# head(reshape2::dcast(coverage_df, average*std_dev*count ~sample))

# print(coverage_df)

# make line plots to show stats
pdf(file = file.path(outdir, "avg_cov_chr21-1Mbbp.pdf"), height = 8, width = 8, onefile = TRUE)
chrom_lines <- ggplot(coverage_df, aes(x=region, y=average, group=sample, colour = sample))
chrom_lines <- chrom_lines + geom_line() # size = 1
# chrom_lines <- chrom_lines + geom_ribbon(aes(x=region, ymin=max(0, average - std_dev), ymax=average + std_dev, fill = sample), alpha = 0.05)
# chrom_lines <- chrom_lines + geom_errorbar(aes(ymin=average - std_dev, ymax=average + std_dev), width=0)
chrom_lines <- chrom_lines + labs(title=paste0("Chromosome: chr21\nAverage Coverage"), x="Region", y = "Average Coverage")
chrom_lines <- chrom_lines + theme(axis.text.x  = element_text(angle=90, size=6))
print(chrom_lines)


std_lines <- ggplot(coverage_df, aes(x=region, y=std_dev, group=sample, colour = sample))
std_lines <- std_lines + geom_line() # size = 1
std_lines <- std_lines + labs(title=paste0("Chromosome: chr21\nStandard Deviation per Region"), x="Region", y = "Std Dev")
std_lines <- std_lines + theme(axis.text.x  = element_text(angle=90, size=6))
print(std_lines)

count_lines <- ggplot(coverage_df, aes(x=region, y=count, group=sample, colour = sample))
count_lines <- count_lines + geom_line() # size = 1
count_lines <- count_lines + labs(title=paste0("Chromosome: chr21\nCounts per Region"), x="Region", y = "Counts")
count_lines <- count_lines + theme(axis.text.x  = element_text(angle=90, size=6))
print(count_lines)

dev.off()


# make color pallette for heatmap
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
