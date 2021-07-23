
CC = gcc
CFLAGS= -Wall 
CLIBS=

all: main

main : qr.h aux.h main.c aux.c qr.c householder.c householder.h trelica.h trelica.c
	$(CC)  -o main main.c aux.c qr.c householder.c trelica.c -lm -g -I -O0
