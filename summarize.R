#!/usr/bin/env Rscript


# R commands to summarize and plot the data


####################


# load relevant libraries
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(GGally, quietly = TRUE)
library(chromPlot, quietly = TRUE)


####################


# find all summary files
summary_files = list.files(path = ".", pattern = "summary.*.5000.txt")
summary_files

# import and merge all summary files (619,150 rows for 5kb bins)
imported_files = lapply(summary_files, read.delim, check.names = FALSE, stringsAsFactors = FALSE)
combined_summary = Reduce(function(x, y) merge(x, y, all = FALSE, by = "#BIN"), imported_files, accumulate = FALSE)
dim(combined_summary)
head(combined_summary)


####################


# plot distributions of values

plot_50mer = ggplot(combined_summary, aes(Mapability50mer)) +
geom_histogram(binwidth = 0.01, fill = "dodgerblue4") +
scale_x_continuous() +
scale_y_sqrt(labels = scales::comma, limits = events_y_limits, breaks = events_y_breaks)

plot_100mer = ggplot(combined_summary, aes(Mapability100mer)) +
geom_histogram(binwidth = 0.01, fill = "dodgerblue4") +
scale_x_continuous() +
scale_y_sqrt(labels = scales::comma, limits = events_y_limits, breaks = events_y_breaks)

events_x_limits = c(-1, 80)
events_y_limits = c(0, 600000)
events_y_breaks = seq(0, 600000, 100000)

plot_giab = ggplot(combined_summary, aes(events_GIAB)) +
geom_histogram(binwidth = 1, fill = "dodgerblue4") +
scale_x_continuous(limits = events_x_limits) +
scale_y_sqrt(labels = scales::comma, limits = events_y_limits, breaks = events_y_breaks)

plot_1kg = ggplot(combined_summary, aes(events_1KG)) +
geom_histogram(binwidth = 1, fill = "dodgerblue4") +
scale_x_continuous(limits = events_x_limits) +
scale_y_sqrt(labels = scales::comma, limits = events_y_limits, breaks = events_y_breaks)

# plot the four distribution plots in one image
plot_grid(plot_50mer, plot_100mer, plot_giab, plot_1kg) +
ggsave("distributions.all.png", width = 12, height = 8, units = "in")


####################


# a matrix of correlation plots with regression lines

# calculate correlation for a random subset of entries to save time (skip bin name column)
correlation_columns = c("events_1KG", "events_GIAB", "Mapability100mer", "Mapability50mer")
correlation_subset = combined_summary[sample(nrow(combined_summary), 10000), ]
correlation_subset = as.matrix(correlation_subset[, correlation_columns])
dim(correlation_subset)
head(correlation_subset)

# regression lines
smoothed_mean = function(data, mapping, ...) {
  p = ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method="loess", fill="red", color="red", ...) +
    geom_smooth(method="gam", fill="blue", color="blue", ...)
  p
}

# plot correlations
png("correlations.png", res = 100, width = 15, height = 10, units = "in")
ggpairs(correlation_subset, lower = list(continuous = smoothed_mean))
dev.off()


####################


# calculate scores

scores = combined_summary
scores[,"chr"] = sub("(.*):.*-.*", "\\1",  scores[,"#BIN"])
scores[,"start_pos"] = as.integer(sub(".*:(.*)-.*", "\\1",  scores[,"#BIN"]))
scores[,"end_pos"] = as.integer(sub(".*:.*-(.*)", "\\1",  scores[,"#BIN"]))

# calculate event percentiles after filtering out bins with 0 events
percentiles_giab = quantile(subset(combined_summary, events_GIAB > 0)[,"events_GIAB"], c(0.95, 0.99, 0.999, 0.9999))
score_cap_giab = percentiles_giab["99%"]
percentiles_1kg = quantile(subset(combined_summary, events_1KG > 0)[,"events_1KG"], c(0.95, 0.99, 0.999, 0.9999))
score_cap_1kg = percentiles_1kg["99%"]

# capped break point events
scores[, "events_GIAB_capped"] = scores[, "events_GIAB"]
scores[scores[, "events_GIAB"] > score_cap_giab,][, "events_GIAB_capped"] = score_cap_giab
scores[, "events_1KG_capped"] = scores[, "events_1KG"]
scores[scores[, "events_1KG"] > score_cap_1kg,][, "events_1KG_capped"] = score_cap_1kg

# break points score
breaks_scores_giab = scores[, "events_GIAB_capped"]
breaks_scores_giab = breaks_scores_giab / max(breaks_scores_giab)
breaks_scores_1kg = scores[, "events_1KG_capped"]
breaks_scores_1kg = breaks_scores_1kg / max(breaks_scores_1kg)
scores[, "breaks_score"] = (breaks_scores_giab + breaks_scores_1kg) / 2

# mapability score
scores[, "map_score"] = scores[, "Mapability50mer"] + scores[, "Mapability100mer"]
scores[, "map_score"] = scores[, "map_score"] / 2
scores[, "map_score"] = 1 - scores[, "map_score"]

# DangerTrack score
scores[, "danger_score"] = scores[, "map_score"] + scores[, "breaks_score"]
scores[, "danger_score"] = scores[, "danger_score"] / 2

# export BED file
scores_bed = scores[, c("chr", "start_pos", "end_pos", "#BIN", "danger_score")]
scores_bed[,"danger_score"] = round(scores_bed[,"danger_score"] * 1000, digits = 0)
write.table(scores_bed, file = "dangertrack.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# export bedGraph file
scores_bedgraph = scores_bed[, c("chr", "start_pos", "end_pos", "danger_score")]
write.table(scores_bedgraph, file = "dangertrack.bedgraph", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


####################


# divide issues into bins

scores[,"GRC_issues_binned"] = "5-95%"
scores[scores[,"GRC_issues"] > .95,][, "GRC_issues_binned"] = "95-100%"
scores[scores[,"GRC_issues"] < .05,][, "GRC_issues_binned"] = "0-5%"

scores[,"ENCODE_DAC_blacklisted_binned"] = "5-95%"
scores[scores[,"ENCODE_DAC_blacklisted"] > .95,][, "ENCODE_DAC_blacklisted_binned"] = "95-100%"
scores[scores[,"ENCODE_DAC_blacklisted"] < .05,][, "ENCODE_DAC_blacklisted_binned"] = "0-5%"

# t-test for DangerTrack score of high and low issue bins
scores_grc_pval = t.test(scores[scores[,"GRC_issues"] < .05,][, "danger_score"],
                         scores[scores[,"GRC_issues"] > .95,][, "danger_score"])
scores_dac_pval = t.test(scores[scores[,"ENCODE_DAC_blacklisted"] < .05,][, "danger_score"],
                         scores[scores[,"ENCODE_DAC_blacklisted"] > .95,][, "danger_score"])

# comparison braket
bracket1 = data.frame(a = c(1, 1:3,3), b = c(1.02, 1.05, 1.05, 1.05, 1.02))

# plot 2/2
ggplot(scores, aes(GRC_issues_binned, danger_score)) +
geom_boxplot(fill = "steelblue3", outlier.color = "gray") +
geom_line(data = bracket1, aes(x = a, y = b)) +
annotate("text", x = 2, y = 1.1, label = "p < 2.2e-16") +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
ggsave("score_GRC_issues.png", width = 8, height = 5, units = "in")

ggplot(scores, aes(ENCODE_DAC_blacklisted_binned, danger_score)) +
geom_boxplot(fill = "steelblue3", outlier.color = "gray") +
geom_line(data = bracket1, aes(x = a, y = b)) +
annotate("text", x = 2, y = 1.1, label = "p < 2.2e-16") +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
ggsave("score_ENCODE_DAC_blacklisted.png", width = 8, height = 5, units = "in")


####################


# chromosome plots

data(hg_gap)

chromplot_score = data.frame(Chrom = scores[, "chr"], Start = scores[, "start_pos"], End = scores[, "end_pos"],
                             score = scores[, "danger_score"])

chromplot_issues = read.table(file = "issues.bed", col.names = c("Chrom", "Start", "End"),
                              check.names = FALSE, stringsAsFactors = FALSE)
chromplot_issues$Group = "Issues"
chromplot_issues$Colors = "dodgerblue4"


chromPlot(chr = c(1:9), figCols = 9, yAxis = FALSE, segLwd = 5, noHist = TRUE, maxSegs = 99000, gaps = hg_gap,
          stat = chromplot_score,
          statCol = "score", colStat="firebrick", statName = "score", statTyp = "l", statSumm = "mean",
          segment = chromplot_issues)


####################



# end
