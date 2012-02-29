import sys

from readinput import expand_input
import parameters


def main():
    for ifile in sys.argv[1:]:
        expand_input(ifile, parameters)


if __name__ == '__main__':
    main()

