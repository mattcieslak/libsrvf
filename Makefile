CC=gcc
CFLAGS=-fPIC -g -pg -O0 -Wall -DDP_NBHD_DIM=12
OBJS=dp_grid.o

all: $(OBJS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) dp_mex.mex dp_mex.o
