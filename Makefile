#!/bin/bash

CC=g++-5
CFLAGS= -ansi -O5 -Wall
CFLAGS2 = -std=c++14
LDFLAGS= -ansi -lm -Wall
EXEC=community convert hierarchy pre final_cpp
OBJ1= graph_binary.o community.o
OBJ2= graph.o

all: $(EXEC)

community : $(OBJ1) main_community.o
	$(CC) -o $@ $^ $(LDFLAGS)

convert : $(OBJ2) main_convert.o
	$(CC) -o $@ $^ $(LDFLAGS)

hierarchy : main_hierarchy.o
	$(CC) -o $@ $^ $(LDFLAGS)

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
