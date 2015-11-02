GUI = $(HOME)/gui/src
FLAGS = -g -march=native -O3 -DMACOSX -mtune=native -funroll-loops -Wall -Winit-self -Woverloaded-virtual -Wsign-compare -fno-strict-aliasing -std=c++11 -Wno-array-bounds -Wno-unknown-pragmas -Wno-deprecated -fPIC -DNDEBUG -DBUILDING_geode -Wno-writable-strings
COMMON_LIBS = $(GUI)/app.a -lgeode -lportaudio -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_calib3d -lportaudio -lportmidi
INCS = -I/usr/local/include/OpenEXR -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include -I. -Ibuild/native/release -I. -I/usr/local/include -I$(GUI)
OS2 := $(strip $(shell uname))
ifeq ($(OS2), Darwin)
  LIBS = $(COMMON_LIBS) -framework GLUT -framework OpenGL
  DEFS = -DMACOSX
else
  OS := $(strip $(shell uname -o))
  ifeq ($(OS), GNU/Linux)
    LIBS = $(COMMON_LIBS) -L/usr/local/lib -lglut -lGL -lGLU -lm
  else
    $(error Unknown OS)
  endif
endif

all: cad

clean:
	rm *.o cad

cad.o: cad.h cad.cpp
	clang++ -c $(FLAGS) $(INCS) cad.cpp

geom.o: cad.h geom.h geom.cpp
	clang++ -c $(FLAGS) $(INCS) geom.cpp

read-eval.o: cad.h read-eval.h read-eval.cpp
	clang++ -c $(FLAGS) $(INCS) read-eval.cpp

app.o: cad.h app.h read-eval.h app.cpp
	clang++ -c $(FLAGS) $(INCS) app.cpp

cad: cad.o geom.o read-eval.o app.o 
	clang++ -march=native -g -fPIC -o cad cad.o geom.o read-eval.o app.o $(LIBS)
