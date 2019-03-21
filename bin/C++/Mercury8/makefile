CC = g++
CFLAGS  = -g -Wall -v

ComputeMS1: CMercury8.o ComputeMS1.o FFT.o
	$(CC) $(CFLAGS) -o ComputeMS1 CMercury8.o ComputeMS1.o FFT.o

CMercury8.o:  CMercury8.cpp CMercury8.h FFT.h mercury.h
	$(CC) $(CFLAGS) -c CMercury8.cpp

ComputeMS1.o: ComputeMS1.cpp CMercury8.h
	$(CC) $(CFLAGS) -c ComputeMS1.cpp

FFT.o: FFT.cpp FFT.h
	$(CC) $(CFLAGS) -c FFT.cpp

clean:
	$(RM) ComputeMS1 *.o
