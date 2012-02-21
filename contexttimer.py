""" Using the with keyword to time code snippets. """

from contextlib import contextmanager
import time
import sys

@contextmanager
def ContextTimer(name="Unknown", outf=sys.stdout):
    t0 = time.time()
    outf.write("[%s ..." % name)
    outf.flush()
    
    yield

    t1 = time.time()

    outf.write(" completed (%f seconds)]\n" % (t1 - t0))


def main():
    with ContextTimer("Long loop"):
        s = 0
        for i in xrange(10000000):
            s += i

    print s


if __name__ == '__main__':
    main()
    
