CC=gcc
CFLAGS=-fPIC -g -pg -O2 -Wall -DDP_NBHD_DIM=6
VANILLA_OBJS=dp_grid.o
LABELLED_OBJS=dp_grid_jose.o
HEADERS=dp_grid.h dp_nbhd.h

all: vanilla

vanilla: $(VANILLA_OBJS)

labelled: $(LABELLED_OBJS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(VANILLA_OBJS) $(LABELLED_OBJS) *.mex
