# cad
# needs libs geode, gui, portaudio, portmidi, opencv
# examples:

```cad "slice(0, revolve(xmov(3,square(2) - square(1))))"```

```cad "revolve(xmov(3,square(2) - square(1)))"```

```cad 'extrude(2, mag1(2,thicken(0.05, text("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))))'```

```cad 'extrude(2, mag1(4, offset(0.1, letter("A"))))'```

```cad 'extrude(2, mag1(4, thicken(0.1, letter("A"))))'```

```cad 'mag1(4, thicken(0.1, letter("A")))'```

```cad 'letter("A")'```

```cad 'slice(0, sphere(4) - cube(2))'```

```cad "revolve(xmov(6, contour(vec(-2, -2), vec(-1, 0), vec(-2, 2), vec(0, 1), vec(2, 2), vec(1, 0), vec(2, -2), vec(0, -1))))"```

```cad "contour(vec(-2, -2), vec(-1, 0), vec(-2, 2), vec(0, 1), vec(2, 2), vec(1, 0), vec(2, -2), vec(0, -1))"```

```cad "revolve(xmov(8, circle(4)))"```

```cad "revolve(xmov(4, contour(vec(-1, -1), vec(1, -1), vec(0, 1))))"```

```cad "((sphere(8) + extrude(20, circle(2))) - sphere(7)) - extrude(21, circle(1))"```

```cad "thicken(1, line(vec(-5, 5, 0), vec(-5, -5, 0), vec(5, -5, 0), vec(5, 5, 0)))"```

```cad "extrude(4, cube(1))"```

```cad "(xrot(15, cube(2)) + yrot(15, cube(2)) + zrot(15, cube(2)) + cube(2)) - sphere(2.6)"```

