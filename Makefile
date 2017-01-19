#
# Compile and run MOSEK examples
#

CC=gcc
IPATHS=-I/home/maxdan94/mosek/8/tools/platform/linux64x86/h
LPATH=-L/home/maxdan94/mosek/8/tools/platform/linux64x86/bin -Wl,-rpath-link,/home/maxdan94/mosek/8/tools/platform/linux64x86/bin '-Wl,-rpath=/home/maxdan94/mosek/8/tools/platform/linux64x86/bin'
LIBS=-lmosek64

all: mclp_cor mclp

%.o: %.c
	$(CC) -c -g $(IPATHS) -o $@ $<

mclp: mclp.o
	$(CC) -g $(LPATH) -o $@ $< $(LIBS)

mclp_cor: mclp_cor.o
	$(CC) -g $(LPATH) -o $@ $< $(LIBS)


clean:
	rm -f mclp mclp_cor
	rm -f mclp.o mclp_cor.o
