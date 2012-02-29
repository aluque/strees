""" A python script to submit to the qsub queue. """

import sys
import os.path
from optparse import OptionParser
from subprocess import call

DEF_QUEUE = 'exe-x86_64'

def encode_envlist(d):
    return ','.join('%s=%s' % (key, item) for key, item in d.iteritems())


def submit(ifile, queue):
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
           
    
    call(cmd, shell=True)
    
    
def main():
    parser = OptionParser()

    parser.add_option("--queue", "-q", dest="queue", type="str",
                      help="Queue to submit to", default=DEF_QUEUE)

    (opts, args) = parser.parse_args()

    for ifile in args:
        submit(ifile, queue=opts.queue)



if __name__ == '__main__':
    main()
    
