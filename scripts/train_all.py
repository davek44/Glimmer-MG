#!/usr/bin/env python
from optparse import OptionParser, SUPPRESS_HELP
import glob, os, time, subprocess, sys

################################################################################
# train_all.py
#
# Run train_gbk in parallel on the gbk files defined by the glob
################################################################################

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
genomes_dir = os.path.abspath('%s/../phymm/.genomeData' % scripts_dir)

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-p', dest='proc', type='int', default=2, help='Number of CPUs to utilize [default: %default]')
    parser.add_option('-l','--min_length', dest='min_length', default=0, help='Minimum length of gene (and ORF) to consider [default: %default]')
    parser.add_option('-o','--max_overlap', dest='max_overlap', default=0, help='Maximum overlap of two genes (or gene and ORF) to consider [default: %default]')
    parser.add_option('-u','--undone', dest='undone', default=False, action='store_true', help='Only train for organisms that are not yet done')

    # run on Condor grid
    parser.add_option('--condor', dest='condor', default=False, action='store_true', help=SUPPRESS_HELP)

    (options,args) = parser.parse_args()

    cmds = []
    for gbk_file in glob.glob('%s/*/*.gbk' % genomes_dir):
        if options.min_length:
            ml = '-l %d ' % options.min_length
        else:
            ml = ''
        if options.max_overlap:
            mo = '-o %d' % options.max_overlap
        else:
            mo = ''

        if not options.undone or not os.path.isfile('%s.lengths.genes.txt' % gbk_file[:-4]):
            cmd = '%s/train_features.py %s%s--gbk %s' % (scripts_dir, ml, mo, gbk_file)

            if options.condor:
                cmds.append('runCmd -c "%s"' % cmd)
            else:
                cmds.append(cmd)

    exec_par(cmds, options.proc, print_cmd=True)


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
