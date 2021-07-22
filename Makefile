
CC = gcc
CFLAGS= -Wall -g -I
CLIBS=

all: main

main : qr.h aux.h main.c aux.c qr.c householder.c householder.h
	$(CC)  -g -Wall -o main main.c aux.c qr.c householder.c -lm
