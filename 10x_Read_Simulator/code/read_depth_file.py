#!/usr/bin/env python
# python 2.7

import sys
import csv
import collections

input_file = sys.argv[1]
# print input_file
# sys.exit()


def get_num_cols(infile):
    # get the number of columns in the file
    with open(infile) as file:
        reader = csv.reader(file, delimiter='\t')
        first_row = next(reader)
        num_cols = len(first_row)
    return num_cols

# number of genomes in depth file 
num_genomes = get_num_cols(input_file) - 2
# print "Number of genomes", num_genomes


def get_avg_coverage(infile):
    # get the average coverage per genome in the file
    # num_genomes = get_num_cols(infile) - 2
    total_coverage = 0 
    line_count = 0
    # calculate the average coverage for a column in the file
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            line_count += 1
            coverage = line[2]
            total_coverage = total_coverage + int(coverage)
    avg_coverage = total_coverage / line_count
    return avg_coverage

avg_coverage = get_avg_coverage(input_file)
# print 'Average Coverage', avg_coverage


def get_binned_stats(infile, avg_coverage):
    # separate the input into binned regions
    bin_size = 100 # size of the bins
    bin_count = 0 # bin iterator
    bin_coverage = 0 # starting coverage vale
    postion_stats = [] # placeholder for bin stats
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            if (bin_count <= bin_size):
                # get values from the file line
                chrom = line[0]
                position = line[1]
                bin_position = int(position) / int(bin_size)
                coverage = line[2]
                bin_coverage = int(coverage) + bin_coverage
                if (bin_count == 0):
                    # start of the position
                    position_values = []
                    position_values.append(chrom)
                    # position_values.append(position)
                    position_values.append(bin_position)
                if (bin_count == bin_size):
                    # end of the position
                    # bin_avg_coverage = bin_coverage / bin_count
                    # position_values.append(position)
                    position_values.append(bin_position)
                    # coverage for the binned region / avg cov whole file
                    position_values.append(bin_coverage / avg_coverage) 
                    postion_stats.append(position_values)
                    # RESET
                    bin_count = 0
                    bin_coverage = 0
                    continue
                bin_count += 1
    return postion_stats

postion_stats = get_binned_stats(infile = input_file, avg_coverage = avg_coverage)

# for line in postion_stats:
#     print('\t'.join(map(str,line)))




