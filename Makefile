CC = gcc
CFLAGS = -Wall -g
TARGET = myprogram

all : $(TARGET)

$(TARGET) : main.o math.o utils.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o math.o utils.o

main.o : main.c math.h utils.h
	$(CC) $(CFLAGS) -c main.c

math.o : math.c math.h
	$(CC) $(CFLAGS) -c math.c

utils.o : utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c

clean :
	rm -f *.o $(TARGET)

