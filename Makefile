GUI = $(HOME)/gui/src
COMMON_FLAGS = -g -O3 -march=native -mtune=native -funroll-loops -Wall -Winit-self -Woverloaded-virtual -Wsign-compare -fno-strict-aliasing -std=c++11 -Wno-array-bounds -Wno-unknown-pragmas -Wno-deprecated -fPIC -DNDEBUG -DBUILDING_geode -Wno-writable-strings
COMMON_LIBS = $(GUI)/app.a -lportaudio -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_calib3d -lportaudio -lportmidi
COMMON_INCS = -I/usr/local/include/OpenEXR -I. -I/usr/local/include -I$(GUI)
OS2 := $(strip $(shell uname))
ifeq ($(OS2), Darwin)
  C++ := clang++
  FLAGS = $(COMMON_FLAGS) -DMACOSX
  LIBS = $(COMMON_LIBS) -framework GLUT -framework OpenGL -lgeode
  INCS = $(COMMON_INCS) -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include -I. -Ibuild/native/release
  DEFS = -DMACOSX
else
  C++ := g++
  OS := $(strip $(shell uname -o))
  ifeq ($(OS), GNU/Linux)
    FLAGS = $(COMMON_FLAGS)
    INCS = $(COMMON_INCS) -I/usr/include/python2.7
    LIBS = -L/usr/local/lib $(COMMON_LIBS) -lglut -lGL -lGLU -lm /usr/local/lib/libgeode.so
  else
    $(error Unknown OS)
  endif
endif

all: cad cad.a

clean:
	rm -f *.o cad

cad.o: cad.h cad.cpp
	$(C++) -c $(FLAGS) $(INCS) cad.cpp

hull.o: cad.h hull.h hull.cpp
	$(C++) -c $(FLAGS) $(INCS) hull.cpp

iso-surface.o: cad.h iso-surface.h iso-surface.cpp
	$(C++) -c $(FLAGS) $(INCS) iso-surface.cpp

geom.o: cad.h geom.h geom.cpp geom-interface.h
	$(C++) -c $(FLAGS) $(INCS) geom.cpp

expr.o: cad.h expr.h expr.cpp geom-interface.h
	$(C++) -c $(FLAGS) $(INCS) expr.cpp

expr_stub.o: 
	$(C++) -c $(FLAGS) $(INCS) expr_stub.cpp

path.o: cad.h path.h path.cpp
	$(C++) -c $(FLAGS) $(INCS) path.cpp

octree.o: cad.h expr.h octree.h octree.cpp
	$(C++) -c $(FLAGS) $(INCS) octree.cpp

read-eval.o: cad.h geom.h read-eval.h read-eval.cpp expr.h
	$(C++) -c $(FLAGS) $(INCS) read-eval.cpp

app.o: cad.h app.h read-eval.h app.cpp
	$(C++) -c $(FLAGS) $(INCS) app.cpp

cad: cad.o hull.o iso-surface.o geom.o read-eval.o app.o expr.o octree.o path.o
	$(C++) -march=native -g -fPIC -o cad $^ $(LIBS)

cad.a: cad.o hull.o iso-surface.o geom.o octree.o expr.o path.o expr_stub.o
	ar -v -r -u cad.a $^
