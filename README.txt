Gorman Law
10053193
CPSC 501 Assignment 4 Code Tuning
Dec 8, 2015

Compile with :
	gcc main.c -o main -lm
	gcc fft.c -o fft -lm
	gcc fftTuned.c -o fftTuned -lm

Run with:
	./main GuitarDryTen.wav AutoPark.wav out.wav
	etc.

The compiler optimization used was -O2

After the program is run, it prints out the time it took to run the program in seconds.