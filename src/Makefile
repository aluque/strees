# Makefile para 'vlasov2'
# Alejandro Luque Estepa, 2004

OBJ = multipol.o misc.o efield.o

#PYTHONDIR = /opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/
#NUMERICDIR = -I/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/numpy/core/include/
NUMERICDIR = -I$(shell python -c 'import numpy; print numpy.get_include()')
PYTHONDIR = $(shell python-config --includes)

#CFLAGS    = $(shell cat icc.opt)
CFLAGS    = -O3 -fomit-frame-pointer -fPIC
#CFLAGS    = -g -fPIC
DFLAGS    = -DCOMPILING_HOSTNAME=\"$(HOSTNAME)\"
INCLUDE_DIRS = -I$(HOME)/include $(PYTHONDIR) $(NUMERICDIR)
CC   = gcc

%.o:	%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDE_DIRS) -o $@ -c $< 

mpolar.so: $(OBJ) mpolarmod.c efield.c setup.py
	./setup.py build -f install -f --install-lib .

clean:  
	rm *.o 2> /dev/null

