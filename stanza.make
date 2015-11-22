all: star

cad.a: cad.h cad.cpp geom.h geom.cpp
	make -f Makefile

runtime.o: /Users/jrb/bar/stanza-390/runtime/runtime.c
	gcc -O3 -c /Users/jrb/bar/stanza-390/runtime/runtime.c

main.o: /Users/jrb/bar/stanza-390/runtime/main.c
	gcc -O3 -c /Users/jrb/bar/stanza-390/runtime/main.c

star.s: star.stanza jitbot.stanza geom-lslib.lstanza
	~/bar/jitpcb/bin/sobc -i star.stanza -s star.s

star: runtime.o main.o star.s cad.a
	g++ -O3  runtime.o main.o  cad.a -lgeode -framework GLUT -framework OpenGL star.s -o star -lm
