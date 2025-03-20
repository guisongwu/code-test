CC = gcc
CFLAGS = -Wall -g
TARGET = my_program

all: $(TARGET)

$(TARGET): main.o utils.o math.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o utils.o math.o

main.o: main.c math.h utils.h
	$(CC) $(CFLAGS) -c main.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c

math.o: math.c math.h
	$(CC) $(CFLAGS) -c math.c

clean:
	rm -f *.o $(TARGET)
