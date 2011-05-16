#!/usr/bin/env python
from optparse import OptionParser, SUPPRESS_HELP
from heapq import heappush, heappop
import os, pdb, glob, sys, subprocess, time

################################################################################
# double_icms.py
#
# Make ICMs trained on pairs of genome's genes, using only the most similar
# genomes according to my Phymm distance.
################################################################################

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
bin_dir = os.path.abspath('%s/../bin' % scripts_dir)

genome_dir = os.path.abspath('%s/../phymm/.genomeData' % scripts_dir)

dist_file = os.path.abspath('%s/../data/phymm_dists.dat' % scripts_dir)
inform_file = os.path.abspath('%s/../data/informative_genomes.txt' % scripts_dir)

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-p', dest='proc', type='int', default=2, help='Number of CPUs to utilize [default: %default]')
    parser.add_option('-t','--top', dest='top', type='int', default=20, help='Number of most similar genomes to create double ICMs for')
    parser.add_option('-r','--replace', dest='replace', default=False, action='store_true', help='Replace existing ICMs')

    # help='Run on Condor grid'
    parser.add_option('--condor', dest='condor', default=False, action='store_true', help=SUPPRESS_HELP)

    (options,args) = parser.parse_args()

    cmds = {}

    # get informative genomes
    inf_genomes = set([line.rstrip() for line in open(inform_file)])

    # for each line in the distance matrix, corresponding to a genome sequence
    dist_open = open(dist_file)
    headers = dist_open.readline().split()
    line = dist_open.readline()
    while line:
        a = line.split()
        genome = a[0]
        (strain,nc_num) = genome.split('|')
        #print genome

        # if the genome is informative and has enough genes to have training
        # data from train_features.py
        if genome in inf_genomes:
            # get all informative distances
            dists = []
            for i in range(1,len(a)):
                if headers[i-1] in inf_genomes and headers[i-1] != genome:
                    heappush(dists, (float(a[i]),headers[i-1]))

            j = 0
            while j < options.top:
                # pop minimum distance
                (d,genome2) = heappop(dists)

                # determine lexicographic order
                if genome < genome2:
                    (strain1,nc_num1) = genome.split('|')
                    (strain2,nc_num2) = genome2.split('|')
                else:
                    (strain2,nc_num2) = genome.split('|')
                    (strain1,nc_num1) = genome2.split('|')

                # check for enough training data
                gene_fa1 = '%s/%s/%s.gene.fasta' % (genome_dir,strain1,nc_num1)
                gene_fa2 = '%s/%s/%s.gene.fasta' % (genome_dir,strain2,nc_num2)
                if os.path.isfile(gene_fa1) and os.path.isfile(gene_fa2):

                    # check for or create necessary directories
                    icm_predir = '%s/%s/%s_2' % (genome_dir,strain1,nc_num1)
                    if not os.path.isdir(icm_predir):
                        os.mkdir(icm_predir)

                    icm_dir = '%s/%s' % (icm_predir,strain2)
                    if not os.path.isdir(icm_dir):
                        os.mkdir(icm_dir)

                    # if we're replacing all, or it doesn't exist
                    if options.replace or not os.path.isfile('%s/%s.gene.icm' % (icm_dir,nc_num2)):
                        # make training command
                        tmp_fa = 'tmp.%s_%s.%s_%s.fasta' % (strain1,nc_num1,strain2,nc_num2)                        
                        cat_cmd = 'cat %s/%s/%s.gene.fasta %s/%s/%s.gene.fasta > %s' % (genome_dir,strain1,nc_num1,genome_dir,strain2,nc_num2,tmp_fa)
                        build_cmd = '%s/build-icm -r %s/%s.gene.icm < %s' % (bin_dir,icm_dir,nc_num2,tmp_fa)

                        # hash to avoid duplicates
                        if options.condor:
                            cmds[(strain1,nc_num1,strain2,nc_num2)] = 'runCmd -c "nfs_%s; %s"' % (cat_cmd,build_cmd)
                        else:
                            cmds[(strain1,nc_num1,strain2,nc_num2)] = '%s; %s' % (cat_cmd,build_cmd)

                    j += 1

        line = dist_open.readline()

    # run training commands in parallel
    exec_par(cmds.values(), options.proc)
    #print '\n'.join(cmds.values())

    # clean up temp gene files
    for tmp_fa in glob.glob('tmp.*_*.*_*.fasta'):
        os.remove(tmp_fa)


############################################################
# exec_par
#
# Execute the commands in the list 'cmds' in parallel, but
# only running 'max_proc' at a time.
############################################################
def exec_par(cmds, max_proc, print_cmd=False):
    total = len(cmds)
    finished = 0
    running = 0
    p = []

    while finished + running < total:
        # launch jobs up to max
        while running < max_proc and finished+running < total:
            if print_cmd:
                print cmds[finished+running]
            p.append(subprocess.Popen(cmds[finished+running], shell=True))
            #print 'Running %d' % p[running].pid
            running += 1

        # are any jobs finished
        new_p = []
        for i in range(len(p)):
            if p[i].poll() != None:
                running -= 1
                finished += 1
            else:
                new_p.append(p[i])

        # if none finished, sleep
        if len(new_p) == len(p):
            time.sleep(1)
        p = new_p

    # wait for all to finish
    for i in range(len(p)):
        p[i].wait()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
