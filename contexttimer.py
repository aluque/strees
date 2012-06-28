""" Using the with keyword to time code snippets. """

from contextlib import contextmanager
import time
import sys

# @contextmanager
# def ContextTimer(name="Unknown", outf=sys.stdout):
#     t0 = time.time()
#     outf.write("[%s ..." % name)
#     outf.flush()
    
#     yield

#     t1 = time.time()

#     outf.write(" completed (%f seconds)]\n" % (t1 - t0))


class ContextTimer(object):
    def __init__(self, name="Unknown", outf=sys.stdout):
        self.name = name
        self.outf = outf


    def __enter__(self):
        if self.outf is not None:
            self.outf.write("[%s ..." % self.name)
            self.outf.flush()

        self.t0 = time.time()

        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.t1 = time.time()
        self.duration = self.t1 - self.t0
    
        self.outf.write(" completed (%f seconds)]\n" % self.duration)



def main():
    with ContextTimer("Long loop") as ct:
        s = 0
        for i in xrange(10000000):
            s += i

    print s, ct.duration


if __name__ == '__main__':
    main()
    
