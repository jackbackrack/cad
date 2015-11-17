#ifndef __IS_ISO_SURFACE__
#define __IS_ISO_SURFACE__

#include "cad.h"
#include <geode/structure/Hashtable.h>

// BASED ON MATT FISHER'S SIMPLE MARCHING CUBES ALGORITHM
// MODIFIED TO WORK WITH GEODE'S PARTICLE/SIMPLEX TREE CODE AND
//          TO PRODUCE A WATER TIGHT MESH
//          TO ALLOW SPECIFYING A TARGET VALUE

struct GridCell {
  IV3 ijk[8];
  TV3 p[8];	//position of each corner of the grid in world space
  T   val[8];	//value of the function at this grid corner
};

class IsoSurface {
 public:
  void PolygonizeGrids(Array<TV3> points, Array<T> distances, int z);
  RawArray<T> FillGrid(Array<T> distances, int z);
  Mesh IsoApproximate(const TV3 &Start, const TV3 &End, T CellSize, Mesh mesh, T target);
  Mesh IsoApproximate(T BoxSize, T CellSize, Mesh mesh, T target);
  void FreeMemory();
  int polygonise(GridCell &Grid);
  int intersection(GridCell &Grid, int from, int to, int base, int dir);
  int edge_index (int dir, IV3 base);
 private:
  int  _XCount, _YCount, _ZCount;
  TV3  _Start, _End, _Diff;
  T    _CellSize;
  Array<TV3>         _Vertices;
  Array<IV3>         _Triangles;
  Hashtable<IV4,int> _Intersections;
};

extern Mesh offset_mesh(T off, Mesh mesh);

#endif
