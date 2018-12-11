#.SUFFIXES:.cc .o
CC=g++
CFLAGS= -Wall
CFLAGS+= -O3
#CFLAGS+= -gstabs+
#CFLAGS+= -pg
CFLAGS+= -lm
CFLAGS+= -lglut
CFLAGS+= -lGLU
CFLAGS+= -lGL

LINKER= $(CC)
LFLAGS= $(CFLAGS)

SRC = main.cxx linkcell.cxx common.cxx spherical.cxx

OBJ = $(SRC:.cpp=.o)

HEADERS= spherical.h Vector3d.h common.h

EXEC = test

%.o:	%.cxx $(HEADERS)
	$(CC) -c $< -o $e $(CFLAGS)

all:	$(HEADERS) $(EXEC)

$(EXEC):	$(OBJ)
		$(LINKER) $(OBJ) -o $@ $(LFLAGS)

clean:
	rm -f *.o; rm -f $(EXEC); rm -f core; rm -f *.out; rm -f *~;

