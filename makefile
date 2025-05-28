CC = gcc

CFLAGS = -lm -Wall

TARGET = malha

all: $(TARGET)

$(TARGET): malha.c
	$(CC) malha.c -o $(TARGET) $(CFLAGS)

clean:
	rm -f $(TARGET)