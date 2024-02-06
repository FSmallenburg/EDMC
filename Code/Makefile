# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops 


# define any libraries to link into executable:
LIBS=-lm


md: EDMC_WCA.c  EDMC_WCA.h
	$(CC) $(CFLAGS) EDMC_WCA.c -o md  $(LIBS)
