#!/usr/bin/env python
from optparse import OptionParser
import glob, os, sys

################################################################################
# informative_genomes.py
#
# Make a list of genomes that are valid for Glimmer training to place in
# informative_genomes.txt
################################################################################

min_adj = 7

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
phymm_dir = os.path.abspath('%s/../phymm' % scripts_dir)

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    inform_open = open('%s/../data/informative_genomes.txt' % scripts_dir, 'w')
    for gbk_file in glob.glob('%s/.genomeData/*/*.gbk' % phymm_dir):
        genome_pre = gbk_file[:-4]
        informative = True

        if not os.path.isfile('%s.gicm' % genome_pre):
            informative = False

        adjs = 0.0
        for line in open('%s.adj_dist.1.-1.genes.txt' % genome_pre):
            adjs += float(line.split()[1])
        if adjs < min_adj:
            informative = False

        adjs = 0.0
        for line in open('%s.adj_dist.-1.1.genes.txt' % genome_pre):
            adjs += float(line.split()[1])
        if adjs < min_adj:
            informative = False

        if informative:
            (strain, nc_num) = genome_pre.split('/')[-2:]
            print >> inform_open, '%s|%s' % (strain,nc_num)
    inform_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
