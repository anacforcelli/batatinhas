
CC = gcc
CFLAGS= -Wall 
CLIBS=

all: main

main : qr.h aux.h main.c aux.c qr.c householder.c householder.h
	$(CC)  -o main main.c aux.c qr.c householder.c -lm -g -I -O0
