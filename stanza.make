COMMON_LIBS = -lgeode -lm
ifeq ($(OS2), Darwin)
  LIBS = $(COMMON_LIBS)  -framework GLUT -framework OpenGL
else
  OS := $(strip $(shell uname -o))
  ifeq ($(OS), GNU/Linux)
    LIBS = -L/usr/local/lib $(COMMON_LIBS) -lglut -lGL -lGLU -lm 
  else
    $(error Unknown OS)
  endif
endif


all: star

cad.a: cad.h cad.cpp geom.h geom.cpp
	make -f Makefile

runtime.o: ../stanza-390/runtime/runtime.c
	gcc -O3 -c ../stanza-390/runtime/runtime.c

main.o: ../stanza-390/runtime/main.c
	gcc -O3 -c ../stanza-390/runtime/main.c

star.s: star.stanza jitbot.stanza geom-lslib.lstanza
	../jitpcb/bin/sobc -i star.stanza -s star.s

star: runtime.o main.o star.s cad.a
	g++ -O3  runtime.o main.o star.s $(LIBS) cad.a -o star
