#
# Compile and run MOSEK examples
#

CC=gcc
IPATHS=-I/home/maxdan94/mosek/8/tools/platform/linux64x86/h
LPATH=-L/home/maxdan94/mosek/8/tools/platform/linux64x86/bin -Wl,-rpath-link,/home/maxdan94/mosek/8/tools/platform/linux64x86/bin '-Wl,-rpath=/home/maxdan94/mosek/8/tools/platform/linux64x86/bin'
LIBS=-lmosek64


%.o: %.c
	$(CC) -c -g $(IPATHS) -o $@ $<

mincut: lp_mincut.o
	$(CC) -g $(LPATH) -o $@ $< $(LIBS)

.DEFAULT: lp_mincut

clean:
	rm -f lp_mincut
	rm -f lp_mincut.o

