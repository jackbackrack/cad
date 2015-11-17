GUI = $(HOME)/gui/src
COMMON_FLAGS = -g -O3 -march=native -mtune=native -funroll-loops -Wall -Winit-self -Woverloaded-virtual -Wsign-compare -fno-strict-aliasing -std=c++11 -Wno-array-bounds -Wno-unknown-pragmas -Wno-deprecated -fPIC -DNDEBUG -DBUILDING_geode -Wno-writable-strings
COMMON_LIBS = $(GUI)/app.a -lportaudio -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_calib3d -lportaudio -lportmidi
COMMON_INCS = -I/usr/local/include/OpenEXR -I. -I/usr/local/include -I$(GUI)
OS2 := $(strip $(shell uname))
ifeq ($(OS2), Darwin)
  FLAGS = $(COMMON_FLAGS) -DMACOSX
  LIBS = $(COMMON_LIBS) -framework GLUT -framework OpenGL -lgeode
  INCS = $(COMMON_INCS) -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include -I. -Ibuild/native/release
  DEFS = -DMACOSX
else
  OS := $(strip $(shell uname -o))
  ifeq ($(OS), GNU/Linux)
    FLAGS = $(COMMON_FLAGS)
    INCS = $(COMMON_INCS) -I/usr/include/python2.7
    LIBS = -L/usr/local/lib $(COMMON_LIBS) -lglut -lGL -lGLU -lm /usr/local/lib/libgeode.so
  else
    $(error Unknown OS)
  endif
endif

all: cad

clean:
	rm -f *.o cad

cad.o: cad.h cad.cpp
	clang++ -c $(FLAGS) $(INCS) cad.cpp

hull.o: cad.h hull.h hull.cpp
	clang++ -c $(FLAGS) $(INCS) hull.cpp

iso-surface.o: cad.h iso-surface.h iso-surface.cpp
	clang++ -c $(FLAGS) $(INCS) iso-surface.cpp

geom.o: cad.h geom.h geom.cpp
	clang++ -c $(FLAGS) $(INCS) geom.cpp

read-eval.o: cad.h read-eval.h read-eval.cpp
	clang++ -c $(FLAGS) $(INCS) read-eval.cpp

app.o: cad.h app.h read-eval.h app.cpp
	clang++ -c $(FLAGS) $(INCS) app.cpp

cad: cad.o hull.o iso-surface.o geom.o read-eval.o app.o 
	clang++ -march=native -g -fPIC -o cad cad.o hull.o iso-surface.o geom.o read-eval.o app.o $(LIBS)
