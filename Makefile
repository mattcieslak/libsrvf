CC=gcc
CFLAGS=-fPIC -g -pg -O2 -Wall -DDP_NBHD_DIM=7
OBJS=dp_grid.o

all: $(OBJS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) dp_mex.mex dp_mex.o
