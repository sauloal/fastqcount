#!/usr/bin/env python
import sys
import os
import gzip
import errno
import math
import argparse
from itertools import islice



def getFmt(FMT):
    CMIN    = None
    CMAX    = None
    VMIN    = None
    q2phred = None

    if FMT in fmts:
        CMIN    = fmts[FMT]['CMIN'   ]
        CMAX    = fmts[FMT]['CMAX'   ]
        VMIN    = fmts[FMT]['VMIN'   ]
        q2phred = fmts[FMT]['q2phred']

    else:
        print "unknown format:", FMT
        print "valid formats are: ", ", ".join(sorted(fmts))
        print """
    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
    ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
    ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
    .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
    LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
    |                         |    |        |                              |                     |
    33                        59   64       73                            104                   126
    0........................26...31.......40                                  S
                             -5....0........9.............................40   X
                                   0........9.............................40   I
                                      3.....9.............................40   J
    0.2......................26...31........41                                 L

    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
       with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
       (Note: See discussion above).
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
        """
        sys.exit(1)

    return CMIN, CMAX, VMIN, q2phred

def printPhred():
    print "   " + " ".join(["%-10s"% x for x in sorted(phred2q)])

    for q in xrange(0, 42):
        print "%2d" % q,
        vals = []

        for t in sorted(phred2q):
            try:
                vals.append( "%.6f  " % phred2q[t][q] )
            except:
                vals.append( " "*10 )

        print " ".join( vals )

def openfile(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rb')
    else:
        return open(filename, 'r')

def stats_fastq(fhd):
    n_line = 0
    n_base = 0
    chars  = [0 for x in xrange(256)]

    for quine in islice(fhd, 3, None, 4):
        n_line += 1

        for c in quine[:-1]:
            chars[ ord(c) ] += 1

        if n_line % 1000000 == 0:
            print "{:16,d}".format(n_line)

    n_base = sum(chars)
    return n_line, n_base, chars

def guessFmt(chars):
    Is    = set([ i if v != 0 else None for i,v in enumerate(chars) ]) - set([None])
    min_i = min(Is)
    max_i = max(Is)

    print "Quals:", ", ".join( [str(i) for i in Is] )
    print "min  : %3d" % min_i
    print "max  : %3d" % max_i

    for FMT in fmts:
        cmin = fmts[FMT]['CMIN']
        cmax = fmts[FMT]['CMAX']
        if cmin == min_i and cmax == max_i:
            print "inference 1"
            print "OUR BEST GUESS IS '%s'" % FMT
            return FMT

    fs1 = []
    for FMT in fmts:
        cmin = fmts[FMT]['CMIN']
        cmax = fmts[FMT]['CMAX']
        if cmin == min_i or cmax == max_i:
            fs1.append( FMT )

    if len(fs1) == 1:
        print "inference 2"
        print "OUR BEST GUESS IS '%s'" % fs1[0]
        return fs1[0]

    fs2 = []
    for FMT in fmts:
        cmin = fmts[FMT]['CMIN']
        cmax = fmts[FMT]['CMAX']
        if cmin >= min_i and cmax <= max_i:
            fs2.append( FMT )

    if len(fs2) == 1:
        print "inference 3"
        print "OUR BEST GUESS IS '%s'" % fs2[0]
        return fs2[0]

    return None

def printstats(FMT, chars, n_base):
    if FMT is None:
        FMT = guessFmt(chars)

    CMIN, CMAX, VMIN, q2phred = (None, None, None, None)

    if FMT is not None:
        print "format:", FMT
        CMIN, CMAX, VMIN, q2phred = getFmt(FMT)

    try:
        c_val       = 0
        trustworthy = True

        min_min = min([fmts[x]['CMIN'] for x in fmts] + [fmts[x]['CMAX'] for x in fmts])
        max_max = max([fmts[x]['CMIN'] for x in fmts] + [fmts[x]['CMAX'] for x in fmts])
        for i in xrange(min_min, max_max+1):
            v       = chars[i]
            c_val  += v
            prop    = v     * 1.0 / n_base * 100
            c_prop  = c_val * 1.0 / n_base * 100
            c_propi = 100.0 - c_prop

            qval    = None
            phred   = None

            if FMT is not None and i < CMIN or i > CMAX:
                qvals  = " " * 3
                phreds = " " * 7

                if v != 0:
                    trustworthy = False

            else:
                qval   = i - CMIN + VMIN
                qvals  = "%3s"   % qval
                phreds = "%7.5f" % q2phred[ qval ]

            bar = "["  + "*" * (int(prop+.5)/2) + " "*(50-(int(prop+.5)/2)) + "]"

            sys.stdout.write( " {:<6s} [{:3d}]: {:14,d} ({:7.3f}% | {:7.3f}% | {:7.3f}%) Q {} Phred {} {}\n".format(repr(chr(i)), i, v, prop, c_prop, c_propi, qvals, phreds, bar) )

            sys.stdout.flush()

        if FMT is not None and not trustworthy:
            print "YOUR SCALE '%s' GOES BEYOND THE QUALITY SCALE. PROBABLY THE FORMAT IS WRONG." % FMT
            print "OUR BEST GUESS IS '%s'" % guessFmt(chars)

    except IOError as e:
        if e.errno == errno.EPIPE:
            # Handle error
            pass

def main():
    parser = argparse.ArgumentParser(description='Calculates stats from fastq files with automatic quality format inference')
    parser.add_argument('-infile' , type=str, nargs='*', action='append',                       help='input file. can be used several times')
    parser.add_argument('-infiles', type=str, nargs='?',                                        help='comma separated list of files')
    parser.add_argument('-fmt'    , type=str, nargs='?', default=None,    choices=sorted(fmts), help='format: ' + ', '.join(sorted(fmts)) + '. Default: auto' )
    parser.add_argument('-partial',                      action='store_true',                   help='print partial values' )

    args      = parser.parse_args()

    fhd       = None
    FMT       = args.fmt
    filenames = []

    if args.infile is not None:
        for inf in args.infile:
            for i in inf:
                filenames.append( i )


    if args.infiles is not None:
        filenames.extend( args.infiles )


    if len( filenames ) != 0:
        print "FILENAMES:", ", ".join( filenames )

        for infile in filenames:
            if not os.path.exists( infile ):
                print "infile %s does not exists" % infile
                sys.exit( 1 )

            if os.path.isdir( infile ):
                print "infile %s is a folder" % infile
                sys.exit( 1 )

    else:
        print "no filename given, reading stdin"
        filenames.append( sys.stdin )



    data = []
    for filename in filenames:
        with openfile(filename) as fhd:
            n_line, n_base, chars = stats_fastq(fhd)

            try:
                sys.stdout.flush()

            except IOError as e:
                if e.errno == errno.EPIPE:
                    # Handle error
                    pass


            if args.partial:
                try:
                    print "filename:", filename
                    print "# seqs  : %14d" % n_line
                    print "# bases : %14d" % n_base
                    sys.stdout.flush()

                except IOError as e:
                    if e.errno == errno.EPIPE:
                        # Handle error
                        pass

                printstats(FMT, chars, n_base)

                if len( filenames ) > 1:
                    print "#" * 100

            data.append( (filename, n_line, n_base, chars) )


    if len( filenames ) > 1:
        n_line_g = sum([d[1] for d in data])
        n_base_g = sum([d[2] for d in data])
        chars_g  = data[0][3]

        for d in data[1:]:
            for c,v in enumerate(d[3]):
                chars_g[c] += v

        print "filename: GLOBAL"
        print "# seqs  : %14d" % n_line_g
        print "# bases : %14d" % n_base_g

        printstats(FMT, chars_g, n_base_g)


fmts = {
    'sanger'    : { 'CMIN': 33, 'CMAX':  74, 'VMIN':  0, 'q2phred' :                  [ (10 ** (q / -10.0))                                              for q in xrange( 0,41)] }, # Sanger
    'solexa'    : { 'CMIN': 59, 'CMAX': 104, 'VMIN': -5, 'q2phred' :                  [ (10 ** (( 10.0 * math.log10( 1.0 + 10 ** (Q / 10.0)) ) / -10.0)) for Q in xrange( 0,46)]   }, # Solexa # http://maq.sourceforge.net/qual.shtml
    'illumina13': { 'CMIN': 64, 'CMAX': 104, 'VMIN':  0, 'q2phred' :                  [ (10 ** (q / -10.0))                                              for q in xrange( 0,41)] }, # Illumina 1.3+
    'illumina15': { 'CMIN': 67, 'CMAX': 104, 'VMIN':  3, 'q2phred' : [-1, -1, None] + [ (10 ** (q / -10.0))                                              for q in xrange( 3,41)] }, # Illumina 1.5+
    'illumina18': { 'CMIN': 33, 'CMAX':  73, 'VMIN':  0, 'q2phred' :                  [ (10 ** (q / -10.0))                                              for q in xrange( 0,42)] }  # Illumina 1.8+
}

if __name__ == '__main__':
    main()
