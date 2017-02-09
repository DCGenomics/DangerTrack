#!/usr/bin/env python
# python 2.7

# CHSL Hackathon
# Oct 24th, 2016

# USAGE: bin_coverages.py data/depth_file.txt

# DESCRIPTION:
# This script will parse a depth of coverage file into binned regions
# and calculate the average coverage per region

# INPUT FILE FORMAT:
# chr1  10003   6656    7998    8002    7996    7997    7998    7421    8000    7999    5620    7999    6887    7998
# chr1  10004   7997    7999    8000    7997    7996    7997    7996    8000    8000    7786    8000    7998    7998
# chr1  10005   8000    8000    8005    7998    7999    8002    7999    8002    8001    8000    8001    8001    8000
# chr1  10006   8001    8004    8005    8003    8003    8005    8002    8003    8002    8000    8003    8002    8003
# chr1  10007   8001    8004    8006    8003    8000    8005    8002    8001    7994    7998    8000    7992    8000
# chr1  10008   7996    8000    8002    7994    7990    7999    7992    7985    7991    7995    7997    7991    7989
# chr1  10009   8005    8007    8009    8006    8006    8007    8005    8006    8004    8003    8006    8005    8005

# <chrom>\t<position>\t<coverage1>\t<coverage2>...<coverage_n>

# EXAMPLE OUTPUT:
# <chrom>\t<start>\t<stop>\t<avg-1, std_dev-1, count-1>\t<avg-2, std_dev-2, count-2>...<avg-n, std_dev-n, count-n>

from __future__ import division
import sys
import gzip
import csv
import collections
import average_coverages as av
import math

class OrderedDefaultDict(collections.OrderedDict, collections.defaultdict):
    def __init__(self, default_factory=None, *args, **kwargs):
        #in python3 you can omit the args to super
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory

def make_bin_region(value, bin_size):
    start_bin = int(value) - (int(value) % bin_size)
    end_bin = start_bin + bin_size
    # sanity check - make sure we binned correctly
    test_bins(value, start_bin, end_bin)
    return [start_bin, end_bin]

def test_bins(value, start_bin, end_bin):
    if (int(value) < int(start_bin)) == True:
        print 'ERROR: value {} smaller than Start bin {}'.format(value, start_bin)
        sys.exit()
    if (int(value) > int(end_bin)) == True:
        print 'ERROR: value {} greater than End bin {}'.format(value, end_bin)
        sys.exit()

def std_dev_pop(count_val, sum_val, sumsquares_val):
    # standard deviation for a population
    std_dev = math.sqrt((count_val * sumsquares_val - sum_val * sum_val)/(count_val * (count_val)))
    return std_dev

def std_dev_sample(count_val, sum_val, sumsquares_val):
    # standard deviation for a sample
    std_dev = math.sqrt((count_val * sumsquares_val - sum_val * sum_val)/(count_val * (count_val - 1)))
    return std_dev

def get_stats_list(counts_list, sums_list, sumsquares_list):
    # calculates the average and std dev for a list of values
    values_list = zip(counts_list, sums_list, sumsquares_list)
    stats_list = [ (total/counts, std_dev_pop(counts, total, sumsquare), counts) for counts, total, sumsquare in values_list]
    # stats_list = [ (total/counts, std_dev_sample(counts, total, sumsquare)) for counts, total, sumsquare in values_list]
    return stats_list



def bin_file_regions(input_file, bin_size):
    # ~~ SETUP ~~~ # 
    # track of bin counts while parsing file
    # n, s0
    # chrom_bin_counts = collections.defaultdict(int)
    # hold the coverage values per bin
    # sum(x), s1
    # chrom_total_coverage = OrderedDefaultDict(list)
    # hold the sum of squares of coverage for each genome
    # sum(x*x), s2
    # chrom_SS_coverage = OrderedDefaultDict(list)
    #
    chrom_bin_dict = OrderedDefaultDict(lambda : OrderedDefaultDict(dict))
    # ~~~ READ FILE ~~ # 
    # read tsv.gzip file
    with gzip.open(input_file,'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        for line in reader:
        # read tsv file
        # with open(input_file) as tsvin:
        # tsvin = csv.reader(tsvin, delimiter='\t')
        # for line in tsvin:
            chrom = line.pop(0)
            position = line.pop(0)
            # convert remaining values to int
            line = map(int, line)
            ss_line = [x*x for x in line]
            # ~~~ BIN REGION ~~ # 
            # set position bins
            position_start_bin, position_end_bin = make_bin_region(position, bin_size)
            position_bin = (chrom, position_start_bin, position_end_bin)
            # ~~~ TABULTATE ~~ # 
            # increment counter
            if not chrom_bin_dict[position_bin]["count"]:
                chrom_bin_dict[position_bin]["count"] = 1
            else :
                chrom_bin_dict[position_bin]["count"] += 1 
            # check counter value
            if chrom_bin_dict[position_bin]["count"] <= bin_size:
                # calculate totals
                if not chrom_bin_dict[position_bin]["totals"]:
                    chrom_bin_dict[position_bin]["totals"] = line
                else :
                    chrom_bin_dict[position_bin]["totals"] = [x + y for x, y in zip(chrom_bin_dict[position_bin]["totals"], line)]
                # calculate sums of squares
                if not chrom_bin_dict[position_bin]["sumsquares"]:
                    chrom_bin_dict[position_bin]["sumsquares"] = ss_line
                else : 
                    chrom_bin_dict[position_bin]["sumsquares"] = [x + y for x, y in zip(chrom_bin_dict[position_bin]["sumsquares"], ss_line)]
            else :
                print 'ERROR: bin exceeded!'
                print 'There might be duplicate entries!'
                print chrom, position, position_bin #, chrom_bin_counts[position_bin]
                sys.exit()
        # calculate avg and std dev
        for position_bin in chrom_bin_dict.keys():
            stats_list = get_stats_list(counts_list = [chrom_bin_dict[position_bin]["count"]] * len(chrom_bin_dict[position_bin]["totals"]), 
                sums_list = chrom_bin_dict[position_bin]["totals"], 
                sumsquares_list = chrom_bin_dict[position_bin]["sumsquares"])
            chrom_bin_dict[position_bin]["stats"] = stats_list
        return chrom_bin_dict

input_file = sys.argv[1]
bin_size = 1000000

if __name__ == '__main__':
    chrom_output = bin_file_regions(input_file, bin_size)
    for position_bin in chrom_output.iterkeys():
        # chr start stop avg,sd,count
        print '\t'.join(map(str,position_bin)) + '\t' + '\t'.join(map(str,["%s,%s,%s" % (avg, sd, count) for avg, sd, count in chrom_output[position_bin]["stats"]]))
