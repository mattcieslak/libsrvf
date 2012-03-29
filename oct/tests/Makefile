CC=gcc
CFLAGS=-g -pg -O0 -Wall
CPPFLAGS=-I/usr/include/CUnit -I..
LIBS=-lcunit -lm
OBJS=main.o dp_grid_tests.o
SRVFOBJS=../dp_grid.o
EXENAME=test

all: $(OBJS)
	$(CC) $(CFLAGS) -o $(EXENAME) $(OBJS) $(SRVFOBJS) $(LIBS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXENAME)
