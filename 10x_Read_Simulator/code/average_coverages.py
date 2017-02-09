#!/usr/bin/env python
# python 2.7

# CHSL Hackathon
# Oct 24th, 2016

# USAGE: average_coverages.py data/depth_file.txt > output_file

# DESCRIPTION:
# This script will read in a coverage file and output the average coverage & standard deviation
# for chromosome per genome (column) in the file

# INPUT FILE FORMAT:
# chr1	10003	6656	7998	8002	7996	7997	7998	7421	8000	7999	5620	7999	6887	7998
# chr1	10004	7997	7999	8000	7997	7996	7997	7996	8000	8000	7786	8000	7998	7998
# chr1	10005	8000	8000	8005	7998	7999	8002	7999	8002	8001	8000	8001	8001	8000
# chr1	10006	8001	8004	8005	8003	8003	8005	8002	8003	8002	8000	8003	8002	8003
# chr1	10007	8001	8004	8006	8003	8000	8005	8002	8001	7994	7998	8000	7992	8000
# chr1	10008	7996	8000	8002	7994	7990	7999	7992	7985	7991	7995	7997	7991	7989
# chr1	10009	8005	8007	8009	8006	8006	8007	8005	8006	8004	8003	8006	8005	8005

# <chrom>\t<position>\t<coverage1>\t<coverage2>...<coverage_n>

# EXAMPLE OUTPUT: (chr avg,sd,count)
# chr1    358.017,1484.33280699   345.6717,1439.86318396  408.4343,1582.01017231  336.343,1444.68650563   348.8892,1509.07956136  396.3538,1584.93141922  364.8105,1509.45335668  370.6484,1501.2279922   426.4161,1559.15368882  338.9206,1426.47922722  339.8465,1475.07425269  405.7164,1570.28442238  401.4027,1576.27424179
# chr10   225.512,1081.8294437    259.2913,1205.27707878  358.3143,1527.11875462  253.6114,1190.62607967  243.9138,1146.35295096  351.894,1424.49961192   239.5175,1115.93045092  250.8465,1159.68964397  287.7425,1181.79290334  152.9809,763.369967405  201.5209,1006.72654508  210.1142,888.036799664  234.5178,1032.77693501
# chr11   6.2814,3.27899588899    5.6928,2.97123344085    11.0703,4.43519536323   9.3923,4.40302177033    9.5595,4.46636986265    23.2982,8.31388457702   13.8852,4.854608219 14.6912,5.32817441156   25.2844,9.0761069099    14.2905,6.31392981827   8.9247,4.32465373296    27.0815,10.4526866283   21.0386,7.84281263578
# chr12   310.3432,1262.29137081  292.4917,1256.22477182  428.0556,1587.51727194  295.9449,1277.27031542  315.2367,1315.58248509  435.7874,1598.80259532  351.5722,1305.06844625  322.2732,1317.19999763  435.8735,1458.81412143  282.0262,1177.00338602  312.2848,1258.64172769  401.1684,1473.4758554   405.7559,1486.29384064

# <chrom>\t<avg-1,std_dev-1>\t<avg-2,std_dev-2>...<avg-n,std_dev-n>

from __future__ import division
import sys
import csv
import collections
import math

input_file = sys.argv[1]

class OrderedDefaultDict(collections.OrderedDict, collections.defaultdict):
    def __init__(self, default_factory=None, *args, **kwargs):
        #in python3 you can omit the args to super
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory


def get_num_cols(infile):
    # get the number of columns in the file
    with open(infile) as file:
        reader = csv.reader(file, delimiter='\t')
        first_row = next(reader)
        num_cols = len(first_row)
    return num_cols

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

def genome_coverage_stats(infile):
    # get average coverage per chromosome per genome
    # EXAMPLE:
    # {genome1: {chr1: 10000, chr2:5000}, genome2:etc.}
    # ~~~~~ # 
    # ~~ SETUP ~~~ # 
    # dict to hold the calculated values
    genome_bin_dict = OrderedDefaultDict(lambda : OrderedDefaultDict(dict))
    # ~~~ READ FILE ~~ # 
    # calculate the total & SS coverages 
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            chrom = line.pop(0)
            position = line.pop(0)
            # set dict values
            # convert remaining values to int
            line = map(int, line)
            # calc sums of squares
            ss_line = [x * x for x in line]
            # count occurences per chrom
            if not genome_bin_dict[chrom]['counts']:
                genome_bin_dict[chrom]['counts'] = [ item for item in [1] for i in range(len(line))]
            else : 
                genome_bin_dict[chrom]['counts'] = [ item + 1 for item in genome_bin_dict[chrom]['counts']]
            # sum coverages
            if not genome_bin_dict[chrom]['coverages']:
                genome_bin_dict[chrom]['coverages'] = line
            else : 
                genome_bin_dict[chrom]['coverages'] = [ x + y for x, y in zip(genome_bin_dict[chrom]['coverages'], line)]
            # sums of squres
            if not genome_bin_dict[chrom]['sumsquares']:
                genome_bin_dict[chrom]['sumsquares'] = ss_line
            else : 
                genome_bin_dict[chrom]['sumsquares'] = [x + y for x, y in zip(genome_bin_dict[chrom]['sumsquares'], ss_line)]
    # ~~~ STATS ~~ # 
    # calculate the averages & std dev
    for chrom in genome_bin_dict.iterkeys():
        genome_bin_dict[chrom]['stats'] = get_stats_list(counts_list = genome_bin_dict[chrom]['counts'], sums_list = genome_bin_dict[chrom]['coverages'], sumsquares_list = genome_bin_dict[chrom]['sumsquares'])
    return genome_bin_dict

if __name__ == '__main__':
    chrom_output = genome_coverage_stats(input_file)
    for chrom in chrom_output.iterkeys():
        # chr avg,sd,count
        print chrom + '\t' + '\t'.join(map(str,["%s,%s,%s" % (avg, sd, count) for avg, sd, count in chrom_output[chrom]["stats"]]))

