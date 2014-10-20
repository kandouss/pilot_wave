CC=gcc
CFLAGS=-c -std=c99 -Wall -O3
LDFLAGS=
LDLIBS=-lm
VPATH=src/
SOURCES=main.c pwave_math.c pwave_io.c
OBJECTS=$(SOURCES:.c=.o)
DEPS=$(SOURCES:.c=.h)
EXECUTABLE=pwave

all: $(SOURCES) $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

$(VPATH)/%.o: $(VPATH)/%.c $(VPATH)/$(DEPS)
	$(CC) $(CFLAGS) -MMD -o $@ $<

$(VPATH)/%.c: 
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm *.o

