#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, sys

################################################################################
# g3-iterated.py
#
# Run glimmer3 on a single organism.
################################################################################

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
bin_dir = os.path.abspath('%s/../bin' % scripts_dir)
elph_bin = os.path.abspath('%s/../ELPH/sources/elph' % scripts_dir)

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <genome> <tag>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='glimmeropts', default='', help='Additional glimmer3 options')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide <genome> fasta file and <tag> to prefix output files')
    else:
        genome_file = args[0]
        tag = args[1]

    num_steps = 8
    
    # step 1
    # Find long, non-overlapping orfs to use as a training set
    print 'Step 1 of %d: Finding long orfs for training' % num_steps
    cmd = '%s/long-orfs -n -t 1.15 %s %s.longorfs' % (bin_dir,genome_file,tag)
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 2
    # Extract the training sequences from the genome file
    print 'Step 2 of %d: Extracting training sequences' % num_steps
    cmd = '%s/extract -t %s %s.longorfs > %s.train' % (bin_dir,genome_file,tag,tag)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 3
    # Build the icm from the training sequences
    print 'Step 3 of %d: Building ICM' % num_steps
    cmd = '%s/build-icm -r %s.icm < %s.train' % (bin_dir,tag,tag)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 4
    # Run first Glimer
    print 'Step 4 of %d: Running first Glimmer3' % num_steps
    cmd = '%s/glimmer3 %s -u -12 -m %s.icm %s %s.run1' % (bin_dir, options.glimmeropts, tag, genome_file, tag)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 5
    # Retrain models
    print 'Step 5 of %d: Retraining' % num_steps
    cmd = '%s/train_features.py --predict %s.run1.predict --seq %s -f' % (scripts_dir, tag, genome_file)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 6
    # Run second Glimmer
    print 'Step 6 of %d: Running second Glimmer3' % num_steps
    cmd = '%s/glimmer3 %s -F %s.run1.features.txt -b %s.run1.motif -m %s.run1.gicm %s %s.run2' % (bin_dir, options.glimmeropts, tag, tag, tag, genome_file, tag)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 7
    # Retrain models
    print 'Step 7 of %d: Retraining' % num_steps
    cmd = '%s/train_features.py --predict %s.run2.predict --seq %s -f' % (scripts_dir, tag, genome_file)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)

    # step 8
    # Run third Glimmer
    print 'Step 8 of %d: Running third Glimmer3' % num_steps
    cmd = '%s/glimmer3 %s -F %s.run2.features.txt -b %s.run2.motif -m %s.run2.gicm %s %s.run2' % (bin_dir, options.glimmeropts, tag, tag, tag, genome_file, tag)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid,0)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
