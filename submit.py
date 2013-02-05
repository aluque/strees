""" A python script to submit to the qsub queue. """

import sys
import os
from optparse import OptionParser
from subprocess import call

DEF_QUEUE = 'exe-x86_64'

def encode_envlist(d):
    return ','.join('%s=%s' % (key, item) for key, item in d.iteritems())


def submit(ifile, queue, onlyprint=False):
    mypath = os.path.split(os.path.realpath(sys.argv[0]))[0]
    ipath = os.path.split(os.path.realpath(ifile))[0]
    
    env = encode_envlist(dict(FMM_INPUT_FILE=os.path.join(ipath, ifile),
                              FMM_PATH=mypath))

    runid = os.path.splitext(os.path.basename(ifile))[0]
    script = os.path.join(mypath, 'qrun.sh')
    
    args = {'-N': runid,
            '-j': 'oe',
            '-q': queue,
            '-o': 'localhost:%s' % os.path.join(ipath, runid + '.out'),
            '-v': env}
    cmd = ('qsub %s %s'
           % (' '.join('%s %s' % (key, item)
                       for key, item in args.iteritems()),
              script))
           
    if not onlyprint:
        call(cmd, shell=True)
    else:
        print(cmd)


def completed(ifile):
    ipath = os.path.split(os.path.realpath(ifile))[0]
    ofile = os.path.join(ipath, runid + '.out')

    if not os.path.exists(ofile):
        return False

    itime, otime = [os.stat(f).st_mtime for f in ifile, ofile]
    if itime > otime:
        return False

    return True


def main():
    parser = OptionParser()

    parser.add_option("--queue", "-q", dest="queue", type="str",
                      help="Queue to submit to", default=DEF_QUEUE)

    parser.add_option("--check-completed", "-c", dest="check", 
                      action="store_true",
                      help="Check for a .out file to see if the code has already run", 
                      default=False)

    parser.add_option("--only-print", "-p", dest="onlyprint", 
                      action="store_true",
                      help="Just print commands, do nothing.", 
                      default=False)

    (opts, args) = parser.parse_args()

    for ifile in args:
        if opts.check and not completed(ifile):
            print("Skipping %s due to an existing and newer .out file"
                  % ifile)
            continue

        submit(ifile, queue=opts.queue, onlyprint=opts.onlyprint)



if __name__ == '__main__':
    main()
    
