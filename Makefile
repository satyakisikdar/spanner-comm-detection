#!/bin/bash

CC=g++
CFLAGS= -ansi -O5 -Wall
CFLAGS2 = -std=c++14
LDFLAGS= -ansi -lm -Wall
EXEC= pre final_cpp


all: $(EXEC)

pre : 
	$(CC) pre_process.cpp -o pre $(CFLAGS) $(CFLAGS2)

final_cpp : 
	$(CC) final_cpp_version.cpp -o final_cpp $(CFLAGS) $(CFLAGS2)

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o *~ $(EXEC)
