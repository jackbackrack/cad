# cad
* simple programmatic 3D solid geometry
* representations for points, polylines, polygons, and meshes
* simple IR for describing operations and is good destination for programming languages

# needs libs
* geode
* gui
* portaudio, portmidi, opencv

# files
* cad.*: cad operations written directly against geode
* geom.*: wrapper around geode reps to be dynamic allow polymorphic construction
* read-eval.*: IR interpreter that read strings into and evals geoms
* app.*: gui front end and viewer

# examples:

* 1D

```cad "line2(vec(-5,5),vec(5,5),vec(5,-5))"```
```cad "line(vec(-5,5,-5),vec(-5,-5,0),vec(5,-5,0),vec(5,5,0))"```
```cad "polyline2(line2(vec(-5,-5),vec(5,-5)),line2(vec(5,5),vec(-5,5)))"```
```cad 'letter("A")'```
```cad 'text("0123456789")'```

* 2D

```cad "square(2)"```
```cad "circle(2)"```
```cad "square(2) \ square(1)"```
```cad "contour(vec(-2,-2),vec(-1,0),vec(-2,2),vec(0,1),vec(2,2),vec(1,0),vec(2,-2),vec(0,-1))"```
```cad "polygon(contour(vec(-2,-2),vec(-1,0),vec(-2,2),vec(0,1),vec(2,2),vec(1,0),vec(2,-2),vec(0,-1)))"```

* 3D

```cad "sphere(2.6)"```
```cad "cube(2)"```
```cad "xrot(15,cube(2))"```
```cad "(xrot(15,cube(2)) | yrot(15,cube(2)) | zrot(15,cube(2)) | cube(2)) \ sphere(2.6)"```
```cad "((sphere(8) | extrude(20,circle(2))) \ sphere(7)) \ extrude(21,circle(1))"```

* 1D -> 2D

```cad 'mag1(4,thicken(0.1,letter("A")))'```
```cad 'mag1(4,offset(0.1,letter("A")))'```

* 1D -> 3D

```cad "thicken(1,line(vec(-5,5,-5),vec(-5,-5,0),vec(5,-5,0),vec(5,5,0)))"```

* 2D -> 3D

```cad "extrude(16,square(2) \ square(1))"```
```cad 'extrude(2,mag1(4,thicken(0.1,letter("A"))))'```
```cad 'extrude(2,mag1(4,offset(0.1,letter("A"))))'```
```cad 'extrude(2,mag1(2,thicken(0.05,text("0123456789"))))'```
```cad 'extrude(2,mag1(2,offset(0.05,text("0123456789"))))'```
```cad "revolve(xmov(6,contour(vec(-2,-2),vec(-1,0),vec(-2,2),vec(0,1),vec(2,2),vec(1,0),vec(2,-2),vec(0,-1))))"```
```cad "revolve(xmov(4,contour(vec(-1,-1),vec(1,-1),vec(0,1))))"```
```cad "revolve(xmov(8,circle(4)))"```
```cad "revolve(xmov(4,square(2) \ square(1)))"```

* 3D -> 2D

```cad 'slice(0,sphere(4) \ cube(2))'```
```cad "slice(0,revolve(xmov(4,square(2) \ square(1))))"```

# TODO

* 3D offsetting
* 3D mesh shadow instead of slice
* 3D->3D revolve/extrude
* SVG reading/writing
* colors
