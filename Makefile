CC=gcc
CFLAGS=-Wall -O3 -I/usr/local/include
LDFLAGS=-lm -lgsl -lgslcblas

SUSYQM_hmc.o: SUSYQM_hmc.c
	${CC} ${CFLAGS} -c SUSYQM_hmc.c

SUSYQM_hmc: SUSYQM_hmc.o
	${CC} -o SUSYQM_hmc SUSYQM_hmc.o ${LDFLAGS}

clean:
	rm *.o