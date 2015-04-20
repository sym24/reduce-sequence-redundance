#!/usr/bin/env python

import csv
import argparse
from subprocess import Popen, STDOUT, PIPE
import os
import sys

def run_blastn(query, reference, outdir, percid, outfile):
    """Run blastn under given percent identity """
    print ">>> Run blastn"
    blastnlog = os.path.join(outdir, 'blastn_db_log')
    # make database and run blastn
    ref = Popen(['makeblastdb', '-dbtype', 'nucl', '-in', reference, '-logfile', blastnlog])
    ref.communicate()
    blastn = Popen(['blastn', '-query', query, '-db', reference, '-out', outfile, \
                    '-best_hit_overhang', '0.2', \
                    '-outfmt', '7', \
                    '-perc_identity', percid, '-evalue', '1e-50', '-num_threads', '12'], stderr=PIPE)
    blastn.communicate()

def filter_query(infile, matchfile):
    """Filter out self to self hit and no hit"""
    # use pipe module. Test in seperate script
    print ">>> Filter query self to self hit and no hit"
    print('>>> Run shell cmd "grep -vw ^# *.blastn | awk $1 != $2 > *matchfile*"')
    grep = Popen(['grep', '-vw', '^#', infile], stdout=PIPE)
    print(repr(['grep', '-vw', '^#', infile]))
    awk = Popen(['awk', '$1 != $2'], stdin=grep.stdout, stdout=PIPE)
    grep.stdout.close()
    output = awk.communicate()[0]

    if grep.returncode == 2:
        print grep.returncode
        sys.exit()

    with open(matchfile, 'wb') as ofile:
        print 'Write to file %s' % matchfile
        ofile.write(output)

def filter_unique_hits(matchfile, outdir):
    """Filter blastn results to reduce redundancy"""
    print ">>> Filter blastn results to reduce redundancy"
    print matchfile
    widict = {}
    wodict = {}
    querydict = {}
    reflist = []
    querylist = []
    woreflist = []
    within_specie = os.path.join(outdir,'unique_within_query.tsv')
    without_specie = os.path.join(outdir, 'unique_without_query.tsv')
    # Open matchfile and output files
    with open(matchfile, 'rb') as fi1, \
             open(within_specie, 'wb') as wi_outfile, \
             open(without_specie, 'wb') as wo_outfile:
        for line in csv.reader(fi1, dialect='excel-tab'):
            querylist.append(line[0])
            reflist.append(line[1])
            if line[0] not in querydict:
                querydict[line[0]] = [[line[1], line[7], line[9]]]
            else:
                querydict[line[0]].append([line[1], line[7], line[9]])

        queryset = sorted(set(querylist))
        # open output file writer handler
        wiwriter = csv.writer(wi_outfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        wowriter = csv.writer(wo_outfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        wiwriter.writerow(["query id", "subject id", "q.end", "s.end"])
        wowriter.writerow(["query id", "subject id", "q.end", "s.end"])
        print 'Write to file %s' % within_specie
        print 'Write to file %s' % without_specie
        for key, comp in querydict.iteritems():
            for entry in comp:
                entry.insert(0, key)
                # entry as a list contain [line[0], line[1], line[7], line[9]]
                # corresponds to ["query id", "subject id", "q.end", "s.end"]
                # if ref in query set, then within species sequence overlap
                if entry[1] in queryset:
                    wiwriter.writerow(entry)
                else:
                    wowriter.writerow(entry)
    return within_specie, without_specie

def wifile_rmlist(within_specie, outdir):
    """find within query specie matches"""
    print ">>> Find wihin query specie matches"
    print ">>> If reference file does not contain query species, should output empty files"
    rmlist = []
    rm_file = os.path.join(outdir, 'query_rmlist.txt')
    with open(rm_file, 'wb') as rmfile, \
            open(within_specie, 'rb') as infile:
        # skip the first line as header
        next(infile)
        for line in csv.reader(infile, dialect='excel-tab'):
            # contains duplicate sequnces in output file
            # if query is longer than ref within one specie
            # remove shorter seq, which is ref
            if line[2] >= line[3]:
                rmlist.append(line[1])
            else:
                rmlist.append(line[0])
        sorted(set(rmlist))
        print 'Write to file %s' % rm_file
        for item in sorted(set(rmlist)):
            rmfile.write("%s\n" % item)
        print 'Number of sequence in %s: %d' % (rm_file, len(set(rmlist)))

# runs long time, don't get to finish. Should try a smaller file
def wofile_rmlist(without_specie, outdir, reference):
    """find without species match sequence"""
    print ">>> Find other species in reference file that match to query sequence"
    exter_seq_log = os.path.join(outdir, 'extract_ref_match_seq_log')
    rm_file1 = os.path.join(outdir, 'other_query_rmlist.txt')
    rm_file2 = os.path.join(outdir, 'reference_rmlist.txt')
    querydict = {}
    othermlist = []
    with open(exter_seq_log, 'wb') as logfile, \
            open(rm_file1, 'wb') as rmfile1, \
            open(rm_file2, 'wb') as rmfile2, \
            open(without_specie, 'rb') as infile:
        next(infile)
        for line in csv.reader(infile, dialect='excel-tab'):
            if line[0] not in querydict:          
                querydict[line[0]] = [line[1]] 
            else:
                querydict[line[0]].append(line[1])
        print 'Write to log file %s' % exter_seq_log
        print 'Write to file %s' % rm_file1
        print 'Number of quenence in %s: %d' % (rm_file1, len(querydict.keys()))
        for key, comp in querydict.iteritems():
            logfile.write('%s\n' % key)
            logfile.flush()
            rmfile1.write('%s\n' % key)
            for entry in comp:
                othermlist.append(entry)
                p = Popen(["grep", "-w", str(entry),
                           reference], stdout=logfile, stderr=STDOUT)
                p.communicate()
                logfile.flush()
        print 'Write to file %s' % rm_file2
        # make sure the uniqueness of sequence header in rmfile for number count
        for item in sorted(set(othermlist)):
            rmfile2.write('%s\n' % item)
        print 'Number of sequence in %s: %d' % (rm_file2, len(set(othermlist)))

def get_parser():
    parser = argparse.ArgumentParser(description='Analyze overlap sequence between given '
                                     'query files and reference files')
    # Required arguments
    parser.add_argument('-query', dest='query', metavar='file',
                        type=argparse.FileType('r'), required=True,
                        help='blastn query file')
    parser.add_argument('-ref', dest='ref', metavar='file', 
                        type=argparse.FileType('r'), required=True, 
                        help='blastn reference sequence file')
    parser.add_argument('-percid', dest='percid', metavar='INT', type=int,
                        help='percent identity for blastn alignment. For details check '
                        'blastn man page. Default [%(default)s]', default='100')

    # Optional variables
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='STR', type=str,
                        help='output file path. Default: [%(default)s]',
                        default=os.getcwd())
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    outdir = os.path.abspath(args.outdir)
    base_query = os.path.basename(args.query.name)
    base_ref = os.path.basename(args.ref.name)
    # blastn outfile name
    bloutfile = '%s_%s_perc%d.blastn.outfile' % (os.path.splitext(base_query)[0], \
                                       os.path.splitext(base_ref)[0], \
                                       args.percid)
    path_bloutfile = os.path.join(outdir, bloutfile)

    # filter query seq outfile name
    matchfile = '%s_match_seq.tsv' % (os.path.splitext(base_query)[0])
    path_matchfile = os.path.join(outdir, matchfile)

    # run blastn
    run_blastn(args.query.name, args.ref.name, outdir, str(args.percid), path_bloutfile)    
    # filter blastn output gain only matching information
    filter_query(path_bloutfile, path_matchfile)
    # run filter unique hits and categorize within specie match, without specie match
    wifile, wofile = filter_unique_hits(path_matchfile, outdir)
    wifile_rmlist(wifile, outdir)
    wofile_rmlist(wofile, outdir, args.ref.name)
    
    print 'finish analyzing' + matchfile
    sys.stdout.flush()
if __name__=='__main__':
    main()
