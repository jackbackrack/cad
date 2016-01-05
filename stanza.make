COMMON_LIBS = -lgeode -lm
OS := $(strip $(shell uname))
ifeq ($(OS), Darwin)
  LIBS = $(COMMON_LIBS)  -framework GLUT -framework OpenGL
else
  ifeq ($(OS), Linux)
    LIBS = -L/usr/local/lib $(COMMON_LIBS) -lglut -lGL -lGLU -lm 
  else
    $(error Unknown OS)
  endif
endif


all: star dicer chair four-bar space-frame

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

dicer.s: dicer.stanza jitbot.stanza geom-lslib.lstanza
	../jitpcb/bin/sobc -i dicer.stanza -s dicer.s

dicer: runtime.o main.o dicer.s cad.a
	g++ -O3  runtime.o main.o dicer.s $(LIBS) cad.a -o dicer

chair.s: chair.stanza jitbot.stanza geom-lslib.lstanza
	../jitpcb/bin/sobc -i chair.stanza -s chair.s

chair: runtime.o main.o chair.s cad.a
	g++ -O3  runtime.o main.o chair.s $(LIBS) cad.a -o chair

four-bar.s: four-bar.stanza jitbot.stanza geom-lslib.lstanza
	../jitpcb/bin/sobc -i four-bar.stanza -s four-bar.s

four-bar: runtime.o main.o four-bar.s cad.a
	g++ -O3  runtime.o main.o four-bar.s $(LIBS) cad.a -o four-bar

space-frame.s: space-frame.stanza jitbot.stanza geom-lslib.lstanza
	../jitpcb/bin/sobc -i space-frame.stanza -s space-frame.s

space-frame: runtime.o main.o space-frame.s cad.a
	g++ -O3  runtime.o main.o space-frame.s $(LIBS) cad.a -o space-frame
