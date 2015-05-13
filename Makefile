CC=gcc
want_factor=false
CFLAGS=-Wall -O3 -I/usr/local/include
ifeq ($(strip ${want_factor}), true)
	CFLAGS += -DFACTOR
endif 
LDFLAGS=-lm -lgsl -lgslcblas	

SUSYQM_hmc.o: SUSYQM_hmc.c
	${CC} ${CFLAGS} -c SUSYQM_hmc.c

SUSYQM_hmc: SUSYQM_hmc.o
	${CC} -g -o SUSYQM_hmc SUSYQM_hmc.o ${LDFLAGS}

QM_hmc.o: QM_hmc.c
	${CC} ${CFLAGS} -c QM_hmc.c

QM_hmc: QM_hmc.o
	${CC} -g -o QM_hmc QM_hmc.o ${LDFLAGS}

clean:
	rm *.o