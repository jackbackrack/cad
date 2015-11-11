/**
  * Copyright John E. Lloyd, 2004. All rights reserved. Permission to use,
  * copy, modify and redistribute is granted, provided that this copyright
  * notice is retained and the author is given credit whenever appropriate.
  *
  * This  software is distributed "as is", without any warranty, including 
  * any implied warranty of merchantability or fitness for a particular
  * use. The author assumes no responsibility for, and shall not be liable
  * for, any special, indirect, or consequential damages, or any damages
  * whatsoever, arising out of or in connection with the use of this
  * software.
  */

#include "hull.h"

class HalfEdge;
class Vertex;
class FaceList;

static const double DOUBLE_PREC = 2.2204460492503131e-16;

typedef Vector<double,3> V3d;
typedef Vector<int,3> V3i;

class HalfEdge;
class Vertex;
class FaceList;

static const int VISIBLE = 1;
static const int NON_CONVEX = 2;
static const int DELETED = 3;

inline V3d v3d (double x, double y, double z) {
  V3d res(x, y, z);
  return res;
}

inline V3i v3i (int x, int y, int z) {
  V3i res(x, y, z);
  return res;
}

class Face {
 private:
  V3d normal;
  V3d centroid;

  void computeNormalAndCentroid(void);
  void computeNormalAndCentroid(double minArea);
  double areaSquared (HalfEdge* hedge0, HalfEdge* hedge1);

 public:

  HalfEdge* he0;
  double area;
  double planeOffset;
  int mark = VISIBLE;
  int index;
  int numVerts;
  Face* next;
  Vertex* outside;

  V3d computeCentroid (void);
  V3d computeNormal (double minArea);
  V3d computeNormal (void);
  static Face* createTriangle (Vertex* v0, Vertex* v1, Vertex* v2);
  static Face* createTriangle (Vertex* v0, Vertex* v1, Vertex* v2, double minArea);
  Face (void) : he0(NULL), mark(VISIBLE), next(NULL), outside(NULL) { }
  HalfEdge* getEdge(int i);
  HalfEdge* getFirstEdge(void) { return he0; }
  HalfEdge* findEdge (Vertex* vt, Vertex* vh);
  double distanceToPlane (V3d p);
  V3d getNormal (void) { return normal; }
  V3d getCentroid (void) { return centroid; }
  int numVertices(void) { return numVerts; }
  std::string getVertexString (void);
  std::vector< int > getVertexIndices (void);
  Face* connectHalfEdges (HalfEdge* hedgePrev, HalfEdge* hedge);
  void checkConsistency(void);
  std::vector<Face*> mergeAdjacentFace (HalfEdge* hedgeAdj);
  void triangulate (FaceList* newFaces, double minArea);
};
      
class FaceList {
 private:
  Face* head;
  Face* tail;
 public:
  void clear(void) { head = tail = NULL; }
  void add (Face* vtx);
  inline Face* first(void) { return head; }
  inline bool isEmpty () { return head == NULL; }
  FaceList (void) : head(NULL), tail(NULL) { }
};

class HalfEdge {
 public:
  Vertex* vertex;
  Face* face;
  HalfEdge* next;
  HalfEdge* prev;
  // Half-edge associated with the opposite triangle adjacent to this edge.
  HalfEdge* opposite;

  // Constructs a HalfEdge with head vertex v and left-hand triangular face f.
  HalfEdge (Vertex* v, Face* f) : vertex(v), face(f), next(NULL), prev(NULL) { }
  HalfEdge (void) : vertex(NULL), face(NULL), next(NULL), prev(NULL)  { }

  // Sets the value of the next edge adjacent (counter-clockwise) to
  // this one within the triangle.
  void setNext (HalfEdge* edge) { next = edge; }
	
  // Gets the value of the next edge adjacent ccw to this one within the triangle.
  HalfEdge* getNext(void) { return next; }

  // Sets the value of the previous edge adjacent cw to this one within the triangle.
  void setPrev (HalfEdge* edge) { prev = edge; }
  // Gets the value of the previous edge adjacent cw to this one within the triangle.
  HalfEdge* getPrev(void) { return prev; }
  // Returns the triangular face located to the left of this half-edge
  Face* getFace(void) { return face; }
  // Returns the half-edge opposite to this half-edge.
  HalfEdge* getOpposite(void) { return opposite; }
  void setOpposite (HalfEdge* edge) {
    opposite = edge;
    edge->opposite = this;
  }
  // Returns the head vertex associated with this half-edge.
  Vertex* head(void) { return vertex; }
  // Returns the tail vertex associated with this half-edge.
  Vertex* tail(void) { return prev != NULL ? prev->vertex : NULL; }
  // Returns the opposite triangular face associated with this half-edge.
  Face* oppositeFace(void) { return opposite != NULL ? opposite->face : NULL; }

  //  Produces a string identifying this half-edge by the point
  //  index values of its tail and head vertices.
  std::string getVertexString(void);

  // Returns the length of this half-edge.
  double length(void);
  double lengthSquared(void);
};

class Vertex {
 public:
  V3d pnt;
  // back index into an array
  int index;
  Vertex* prev;
  Vertex* next;
  Face* face;
 Vertex(double x, double y, double z, int idx) : pnt(v3d(x, y, z)), index(idx) { }
  Vertex (void) : prev(NULL), next(NULL), face(NULL) { }
};

class VertexList {
 private:
  Vertex* head;
  Vertex* tail;

 public:
  VertexList (void) : head(NULL), tail(NULL) { }
  void clear(void) { head = tail = NULL; }

  void add (Vertex* vtx);
  // Adds a chain of vertices to the end of this list.
  void addAll (Vertex* vtx);
  // Deletes a vertex from this list.
  void del (Vertex* vtx);
  // Deletes a chain of vertices from this list.
  void del (Vertex* vtx1, Vertex* vtx2);
  // Inserts a vertex into this list before another specificed vertex.
  void insertBefore (Vertex* vtx, Vertex* next);
  // Returns the first element in this list.
  Vertex* first(void) { return head; }

  bool isEmpty(void) { return head == NULL; }
};

static const int CLOCKWISE = 0x1;
static const int POINT_RELATIVE = 0x8;
static const double AUTOMATIC_TOLERANCE = -1;
static const int NONCONVEX_WRT_LARGER_FACE = 1;
static const int NONCONVEX = 2;

class QuickHull3D {
 private:
  std::vector<Face*> discardedFaces;
  std::vector<Vertex*> maxVtxs;
  std::vector<Vertex*> minVtxs;
  FaceList* newFaces;
  VertexList* unclaimed;
  VertexList* claimed;

  void addPointToFace (Vertex* vtx, Face* face);
  void removePointFromFace (Vertex* vtx, Face* face);
  Vertex* removeAllPointsFromFace (Face* face);
  HalfEdge* findHalfEdge (Vertex* tail, Vertex* head);
  std::vector< int > getFaceIndices (Face* face, int flags);
  bool doAdjacentMerge (Face* face, int mergeType);
  HalfEdge* addAdjoiningFace (Vertex* eyeVtx, HalfEdge* he);
  void markFaceVertices (Face* face, int mark);

 protected:
  int findIndex = -1;
  bool debug = false;
  double charLength;
  std::vector<Vertex*> pointBuffer;
  std::vector<int> vertexPointIndices;
  std::vector<Face*> faces;
  std::vector<HalfEdge*> horizon;
  int numVertices;
  int numFaces;
  int numPoints;
  double explicitTolerance = AUTOMATIC_TOLERANCE;
  double tolerance;

  void initBuffers (int nump);
  void setPoints (std::vector<double> coords);
  void setPoints (std::vector<V3d> pnts);
  void computeMaxAndMin (void);
  // Creates the initial simplex from which the hull will be built.
  void createInitialSimplex ( void );
  void resolveUnclaimedPoints (FaceList* newFaces);
  void deleteFacePoints (Face* face, Face* absorbingFace);
  double oppFaceDistance (HalfEdge* he);
  void calculateHorizon (V3d eyePnt, HalfEdge* edge0, Face* face, std::vector<HalfEdge*>& horizon);
  void addNewFaces (FaceList* newFaces, Vertex* eyeVtx, std::vector<HalfEdge*>& horizon);
  void addPointToHull(Vertex* eyeVtx);
  void buildHull (void);
  void reindexFacesAndVertices(void);
  bool checkFaceConvexity (Face* face, double tol, std::ostream& ps);
  bool checkFaces(double tol, std::ostream& ps);
  Vertex* nextPointToAdd(void);
  
 public:
  bool getDebug() { return debug; }
  void setDebug (bool enable) { debug = enable; }
  double getDistanceTolerance ( void ) { return tolerance; }
  void setExplicitDistanceTolerance(double tol) { explicitTolerance = tol; }
  double getExplicitDistanceTolerance() { return explicitTolerance; }
  void build (std::vector<V3d> points);
  QuickHull3D (void) { }
  QuickHull3D (std::vector<V3d> points) { build(points); }

  // Triangulates any non-triangular hull faces. In some cases, due to
  // precision issues, the resulting triangles may be very thin or small,
  // and hence appear to be non-convex (this same limitation is present
  // in <a href=http://www.qhull.org>qhull</a>).
  void triangulate (void);
  int getNumVertices() { return numVertices; }

  std::vector<V3d> getVertices(void);
  // Returns an array specifing the index of each hull vertex
  // with respect to the original input points.
  std::vector<int> getVertexPointIndices(void);
  int getNumFaces(void) { return faces.size(); }

  // Returns the faces associated with this hull, where
  // Each face is represented by an integer array which gives the
  // indices of the vertices. These indices are numbered
  // relative to the hull vertices, are zero-based,
  // and are arranged counter-clockwise. More control
  // over the index format can be obtained using
  std::vector< std::vector< int > > getFaces () { return getFaces(0); }

  // Returns the faces associated with this hull.
  //
  // <p>Each face is represented by an integer array which gives the
  // indices of the vertices. By default, these indices are numbered with
  // respect to the hull vertices (as opposed to the input points), are
  // zero-based, and are arranged counter-clockwise. However, this
  // can be changed by setting {@link #POINT_RELATIVE
  // POINT_RELATIVE}, or // {@link #CLOCKWISE CLOCKWISE}
  // in the indexFlags parameter.
  //
  // @param indexFlags specifies index characteristics (0 results
  // in the default)
  // @return array of integer arrays, giving the vertex
  // indices for each face.
  // @see QuickHull3D#getVertices()
  std::vector< std::vector< int > > getFaces (int indexFlags);

  // Checks the correctness of the hull. This is done by making sure that
  // no faces are non-convex and that no points are outside any face.
  // These tests are performed using the distance tolerance <i>tol</i>.
  // Faces are considered non-convex if any edge is non-convex, and an
  // edge is non-convex if the centroid of either adjoining face is more
  // than <i>tol</i> above the plane of the other face. Similarly,
  // a point is considered outside a face if its distance to that face's
  // plane is more than 10 times <i>tol</i>.
  //
  // <p>If the hull has been {@link #triangulate triangulated},
  // then this routine may fail if some of the resulting
  // triangles are very small or thin.
  bool check (std::ostream& ps);
  bool check (std::ostream& ps, double tol);
  Mesh to_mesh (void);
};

extern void error(std::string);

void Face::computeNormalAndCentroid(void) {
  normal = computeNormal ();
  centroid = computeCentroid ();
  planeOffset = dot(normal, centroid);
  int numv = 0;
  HalfEdge* he = he0;
  do {
    numv++;
    he = he->next;
  } while (he != he0);
  if (numv != numVerts) {
    std::stringstream ss;
    ss << "face " + getVertexString() << " numVerts=" << numVerts << " should be " << numv;
    error(ss.str());
  }
}

void Face::computeNormalAndCentroid(double minArea) {
  normal = computeNormal (minArea);
  centroid = computeCentroid ();
  planeOffset = dot(normal, centroid);
}

double Face::areaSquared (HalfEdge* hedge0, HalfEdge* hedge1) {
  // return the squared area of the triangle defined
  // by the half edge hedge0 and the point at the
  // head of hedge1.

  auto p0 = hedge0->tail()->pnt;
  auto p1 = hedge0->head()->pnt;
  auto p2 = hedge1->head()->pnt;

  double dx1 = p1.x - p0.x;
  double dy1 = p1.y - p0.y;
  double dz1 = p1.z - p0.z;

  double dx2 = p2.x - p0.x;
  double dy2 = p2.y - p0.y;
  double dz2 = p2.z - p0.z;

  double x = dy1*dz2 - dz1*dy2;
  double y = dz1*dx2 - dx1*dz2;
  double z = dx1*dy2 - dy1*dx2;

  return x*x + y*y + z*z;	   
}

V3d Face::computeCentroid (void) {
  V3d centroid;
  HalfEdge* he = he0;
  do {
    centroid += he->head()->pnt;
    he = he->next;
  } while (he != he0);
  return 1/(double)numVerts * centroid;
}

V3d Face::computeNormal (double minArea) {
  auto normal = computeNormal();
  if (area < minArea) {
    // make the normal more robust by removing
    // components parallel to the longest edge
    HalfEdge* hedgeMax = NULL;
    double lenSqrMax = 0;
    HalfEdge* hedge = he0;
    do {
      double lenSqr = hedge->lengthSquared();
      if (lenSqr > lenSqrMax) {
        hedgeMax = hedge;
        lenSqrMax = lenSqr;
      }
      hedge = hedge->next;
    } while (hedge != he0);

    V3d p2 = hedgeMax->head()->pnt;
    V3d p1 = hedgeMax->tail()->pnt;
    double lenMax = sqrt(lenSqrMax);
    double ux = (p2.x - p1.x)/lenMax;
    double uy = (p2.y - p1.y)/lenMax;
    double uz = (p2.z - p1.z)/lenMax;	   
    double dot = normal.x*ux + normal.y*uy + normal.z*uz;
    normal.x -= dot*ux;
    normal.y -= dot*uy;
    normal.z -= dot*uz;

    normal.normalize();	      
  }
  return normal;
}

V3d Face::computeNormal (void) {
  HalfEdge* he1 = he0->next;
  HalfEdge* he2 = he1->next;

  V3d p0 = he0->head()->pnt;
  V3d p2 = he1->head()->pnt;

  double d2x = p2.x - p0.x;
  double d2y = p2.y - p0.y;
  double d2z = p2.z - p0.z;

  V3d normal;

  numVerts = 2;

  while (he2 != he0) { 
    double d1x = d2x;
    double d1y = d2y;
    double d1z = d2z;

    p2 = he2->head()->pnt;
    d2x = p2.x - p0.x;
    d2y = p2.y - p0.y;
    d2z = p2.z - p0.z;

    normal.x += d1y*d2z - d1z*d2y;
    normal.y += d1z*d2x - d1x*d2z;
    normal.z += d1x*d2y - d1y*d2x;

    he1 = he2;
    he2 = he2->next;
    numVerts++;
  }
  area = normal.magnitude();
  return (1/area) * normal;
}

Face* Face::createTriangle (Vertex* v0, Vertex* v1, Vertex* v2) {
  return createTriangle (v0, v1, v2, 0);
}

Face* Face::createTriangle (Vertex* v0, Vertex* v1, Vertex* v2, double minArea) {
  Face* face = new Face();
  HalfEdge* he0 = new HalfEdge (v0, face);
  HalfEdge* he1 = new HalfEdge (v1, face);
  HalfEdge* he2 = new HalfEdge (v2, face);

  he0->prev = he2;
  he0->next = he1;
  he1->prev = he0;
  he1->next = he2;
  he2->prev = he1;
  he2->next = he0;

  face->he0 = he0;

  // compute the normal and offset
  face->computeNormalAndCentroid(minArea);
  return face;
}

HalfEdge* Face::getEdge (int i) {
  HalfEdge* he = he0;
  while (i > 0) {
    he = he->next;
    i--;
  }
  while (i < 0) {
    he = he->prev;
    i++;
  }
  return he;
}

// finds the half-edge within this face which has tail vt and head vh
HalfEdge* Face::findEdge (Vertex* vt, Vertex* vh) {
  HalfEdge* he = he0;
  do {
    if (he->head() == vh && he->tail() == vt) {
      return he;
    }
    he = he->next;
  } while (he != he0);
  return NULL;
}

double Face::distanceToPlane (V3d p) {
  return normal.x*p.x + normal.y*p.y + normal.z*p.z - planeOffset;
}

std::string Face::getVertexString (void) {
  std::stringstream ss;
  HalfEdge* he = he0;
  bool is_first = true;
  do {
    if (is_first) {
      ss << he->head()->index;
      is_first = false;
    } else {
      ss << " " << he->head()->index;
    }
    he = he->next;
  } while (he != he0);
  return ss.str();
}

std::vector<int> Face::getVertexIndices (void) {
  HalfEdge* he = he0;
  std::vector<int> indices;
  do {
    indices.push_back(he->head()->index);
    he = he->next;
  } while (he != he0);
  return indices;
}

Face* Face::connectHalfEdges (HalfEdge* hedgePrev, HalfEdge* hedge) {
  Face* discardedFace = NULL;
  if (hedgePrev->oppositeFace() == hedge->oppositeFace()) { 
    // then there is a redundant edge that we can get rid off
    Face* oppFace = hedge->oppositeFace();
    HalfEdge* hedgeOpp;

    if (hedgePrev == he0) {
      he0 = hedge; 
    }
    if (oppFace->numVertices() == 3) {
      // then we can get rid of the opposite face altogether
      hedgeOpp = hedge->getOpposite()->prev->getOpposite();
        
      oppFace->mark = DELETED;
      discardedFace = oppFace;
    } else {
      hedgeOpp = hedge->getOpposite()->next;

      if (oppFace->he0 == hedgeOpp->prev) {
        oppFace->he0 = hedgeOpp; 
      }
      hedgeOpp->prev = hedgeOpp->prev->prev;
      hedgeOpp->prev->next = hedgeOpp;
    }
    hedge->prev = hedgePrev->prev;
    hedge->prev->next = hedge;

    hedge->opposite = hedgeOpp;
    hedgeOpp->opposite = hedge;

    // oppFace was modified, so need to recompute
    oppFace->computeNormalAndCentroid();
  } else {
    hedgePrev->next = hedge;
    hedge->prev = hedgePrev;
  }
  return discardedFace;
}

void Face::checkConsistency(void) {
  // do a sanity check on the face
  HalfEdge* hedge = he0; 
  double maxd = 0;
  int numv = 0;

  if (numVerts < 3) {
    std::stringstream ss;
    ss << "degenerate face: " << getVertexString();
    error(ss.str());
  }
  do {
    HalfEdge* hedgeOpp = hedge->getOpposite();
    if (hedgeOpp == NULL) {
      std::stringstream ss;
      ss << "face " << getVertexString() << ": " << "unreflected half edge " << hedge->getVertexString();
      error(ss.str());
    } else if (hedgeOpp->getOpposite() != hedge) {
      std::stringstream ss;
      ss << "face " << getVertexString() << ": " <<
             "opposite half edge " << hedgeOpp->getVertexString() <<
             " has opposite " << hedgeOpp->getOpposite()->getVertexString();
      error(ss.str());
    }
    if (hedgeOpp->head() != hedge->tail() || hedge->head() != hedgeOpp->tail()) {
      std::stringstream ss;
      ss << "face " << getVertexString() << ": " <<
             "half edge " << hedge->getVertexString() <<
             " reflected by " << hedgeOpp->getVertexString();
      error(ss.str());
    }
    Face* oppFace = hedgeOpp->face;
    if (oppFace == NULL) {
      std::stringstream ss;
      ss << "face " << getVertexString() << ": " << "no face on half edge " << hedgeOpp->getVertexString();
      error (ss.str());
    } else if (oppFace->mark == DELETED) {
      std::stringstream ss;
      ss << "face " << getVertexString() << ": " << "opposite face " << oppFace->getVertexString() << " not on hull";
      error (ss.str());
    }
    double d = fabs(distanceToPlane(hedge->head()->pnt));
    if (d > maxd) {
      maxd = d;
    }
    numv++;
    hedge = hedge->next;
  } while (hedge != he0);

  if (numv != numVerts) {
    std::stringstream ss;
    ss << "face " << getVertexString() << " numVerts=" << numVerts << " should be " << numv;
    error(ss.str());
  }
}

std::vector<Face*> Face::mergeAdjacentFace (HalfEdge* hedgeAdj) {
  Face* oppFace = hedgeAdj->oppositeFace();
  std::vector<Face*> discarded;

  discarded.push_back(oppFace);
  oppFace->mark = DELETED;

  HalfEdge* hedgeOpp = hedgeAdj->getOpposite();

  HalfEdge* hedgeAdjPrev = hedgeAdj->prev;
  HalfEdge* hedgeAdjNext = hedgeAdj->next;
  HalfEdge* hedgeOppPrev = hedgeOpp->prev;
  HalfEdge* hedgeOppNext = hedgeOpp->next;

  while (hedgeAdjPrev->oppositeFace() == oppFace) {
    hedgeAdjPrev = hedgeAdjPrev->prev;
    hedgeOppNext = hedgeOppNext->next;
  }
	   
  while (hedgeAdjNext->oppositeFace() == oppFace) {
    hedgeOppPrev = hedgeOppPrev->prev;
    hedgeAdjNext = hedgeAdjNext->next;
  }

  HalfEdge* hedge;

  for (hedge=hedgeOppNext; hedge!=hedgeOppPrev->next; hedge=hedge->next) {
    hedge->face = this;
  }

  if (hedgeAdj == he0) {
    he0 = hedgeAdjNext; 
  }
	   
  // handle the half edges at the head
  Face* discardedFace;

  discardedFace = connectHalfEdges (hedgeOppPrev, hedgeAdjNext);
  if (discardedFace != NULL) {
    discarded.push_back(discardedFace); 
  }

  // handle the half edges at the tail
  discardedFace = connectHalfEdges (hedgeAdjPrev, hedgeOppNext);
  if (discardedFace != NULL) {
    discarded.push_back(discardedFace); 
  }

  computeNormalAndCentroid ();
  checkConsistency();

  return discarded;
}

void Face::triangulate (FaceList* newFaces, double minArea) {
  HalfEdge* hedge;

  if (numVertices() < 4) {
    return; 
  }

  Vertex* v0 = he0->head();

  hedge = he0->next;
  HalfEdge* oppPrev = hedge->opposite;
  Face* face0 = NULL;

  for (hedge=hedge->next; hedge!=he0->prev; hedge=hedge->next) {
    Face* face = createTriangle (v0, hedge->prev->head(), hedge->head(), minArea);
    face->he0->next->setOpposite (oppPrev);
    face->he0->prev->setOpposite (hedge->opposite);
    oppPrev = face->he0;
    newFaces->add (face);
    if (face0 == NULL) {
      face0 = face; 
    }
  }
  hedge = new HalfEdge (he0->prev->prev->head(), this);
  hedge->setOpposite (oppPrev);

  hedge->prev = he0;
  hedge->prev->next = hedge;

  hedge->next = he0->prev;
  hedge->next->prev = hedge;

  computeNormalAndCentroid (minArea);
  checkConsistency();

  for (Face* face=face0; face!=NULL; face=face->next) {
    face->checkConsistency(); 
  }
}
      
void FaceList::add (Face* vtx) {
  if (head == NULL) {
    head = vtx;
  } else {
    tail->next = vtx; 
  }
  vtx->next = NULL;
  tail = vtx;
}


std::string HalfEdge::getVertexString(void) {
  std::stringstream ss;
  if (tail() != NULL) {
    ss << tail()->index << "-" << head()->index;
  } else {
    ss << "?-" << head()->index;
  }
  return ss.str();
}

// Returns the length of this half-edge.
double HalfEdge::length(void) {
  if (tail() != NULL) {
    return (head()->pnt - tail()->pnt).magnitude();
  } else {
    return -1; 
  }
}

double HalfEdge::lengthSquared(void) {
  if (tail() != NULL) {
    return (head()->pnt - tail()->pnt).sqr_magnitude();
  } else {
    return -1; 
  }
}

void VertexList::add (Vertex* vtx) {  
  if (head == NULL) {
    head = vtx;
  } else {
    tail->next = vtx; 
  }
  vtx->prev = tail;
  vtx->next = NULL;
  tail = vtx;
}

// Adds a chain of vertices to the end of this list.
void VertexList::addAll (Vertex* vtx) { 
  if (head == NULL) {
    head = vtx;
  } else {
    tail->next = vtx; 
  }
  vtx->prev = tail;
  while (vtx->next != NULL) {
    vtx = vtx->next;
  }
  tail = vtx;
}

// Deletes a vertex from this list.
void VertexList::del (Vertex* vtx) {
  if (vtx->prev == NULL) {
    head = vtx->next;
  } else {
    vtx->prev->next = vtx->next; 
  }
  if (vtx->next == NULL) {
    tail = vtx->prev; 
  } else {
    vtx->next->prev = vtx->prev; 
  }
}

// Deletes a chain of vertices from this list.
void VertexList::del (Vertex* vtx1, Vertex* vtx2) {
  if (vtx1->prev == NULL) {
    head = vtx2->next;
  } else {
    vtx1->prev->next = vtx2->next; 
  }
  if (vtx2->next == NULL) {
    tail = vtx1->prev; 
  } else {
    vtx2->next->prev = vtx1->prev; 
  }
}

// Inserts a vertex into this list before another specificed vertex.
void VertexList::insertBefore (Vertex* vtx, Vertex* next) {
  vtx->prev = next->prev;
  if (next->prev == NULL) {
    head = vtx;
  } else {
    next->prev->next = vtx; 
  }
  vtx->next = next;
  next->prev = vtx;
}

void QuickHull3D::addPointToFace (Vertex* vtx, Face* face) {
  vtx->face = face;

  if (face->outside == NULL) {
    claimed->add (vtx);
  } else {
    claimed->insertBefore (vtx, face->outside); 
  }
  face->outside = vtx;
}

void QuickHull3D::removePointFromFace (Vertex* vtx, Face* face) {
  if (vtx == face->outside) {
    if (vtx->next != NULL && vtx->next->face == face) {
      face->outside = vtx->next;
    } else {
      face->outside = NULL; 
    }
  }
  claimed->del (vtx);
}

Vertex* QuickHull3D::removeAllPointsFromFace (Face* face) {
  if (face->outside != NULL) { 
    Vertex* end = face->outside;
    while (end->next != NULL && end->next->face == face) {
      end = end->next;
    }
    claimed->del (face->outside, end);
    end->next = NULL;
    return face->outside;
  } else {
    return NULL; 
  }
}

HalfEdge* QuickHull3D::findHalfEdge (Vertex* tail, Vertex* head) { 
  for (auto face : faces) {
    HalfEdge* he = face->findEdge (tail, head);
    if (he != NULL) {
      return he; 
    }
  }
  return NULL;
}
  
std::vector<int> QuickHull3D::getFaceIndices (Face* face, int flags) { 
  std::vector<int> indices;
  bool ccw = ((flags & CLOCKWISE) == 0);
  bool pointRelative = ((flags & POINT_RELATIVE) != 0);

  HalfEdge* hedge = face->he0;
  do {
    int idx = hedge->head()->index;
    if (pointRelative) {
      idx = vertexPointIndices[idx];
    }
    indices.push_back(idx);
    hedge = (ccw ? hedge->next : hedge->prev);
  } while (hedge != face->he0);	   
  return indices;
}

bool QuickHull3D::doAdjacentMerge (Face* face, int mergeType) {
  HalfEdge* hedge = face->he0;

  bool convex = true;
  do {
    Face* oppFace = hedge->oppositeFace();
    bool merge = false;
    double dist1;

    if (mergeType == NONCONVEX) {
      // then merge faces if they are definitively non-convex
      if (oppFaceDistance (hedge) > -tolerance ||
          oppFaceDistance (hedge->opposite) > -tolerance) {
        merge = true;
      }
    } else { // mergeType == NONCONVEX_WRT_LARGER_FACE 
      // merge faces if they are parallel or non-convex
      // wrt to the larger face; otherwise, just mark
      // the face non-convex for the second pass.
      if (face->area > oppFace->area) {
        if ((dist1 = oppFaceDistance (hedge)) > -tolerance) {
          merge = true;
        } else if (oppFaceDistance (hedge->opposite) > -tolerance) {
          convex = false;
        }
      } else {
        if (oppFaceDistance (hedge->opposite) > -tolerance) {
          merge = true;
        } else if (oppFaceDistance (hedge) > -tolerance) {
          convex = false;
        }
      }
    }

    if (merge) {
      if (debug) {
        std::cout << "  merging " << face->getVertexString() << "  and  " <<
          oppFace->getVertexString() << "\n";
      }

      auto discardedFaces = face->mergeAdjacentFace (hedge);
      for (auto discard : discardedFaces) {
        deleteFacePoints (discard, face);
      }
      if (debug) {
        std::cout << "  result: " << face->getVertexString() << "\n";
      }
      return true;
    }
    hedge = hedge->next;
  }
  while (hedge != face->he0);
  if (!convex) {
    face->mark = NON_CONVEX; 
  }
  return false;
}

HalfEdge* QuickHull3D::addAdjoiningFace (Vertex* eyeVtx, HalfEdge* he) { 
  Face* face = Face::createTriangle (eyeVtx, he->tail(), he->head());
  faces.push_back (face);
  face->getEdge(-1)->setOpposite(he->getOpposite());
  return face->getEdge(0);
}

void QuickHull3D::markFaceVertices (Face* face, int mark) {
  HalfEdge* he0 = face->getFirstEdge();
  HalfEdge* he = he0;
  do {
    he->head()->index = mark;
    he = he->next;
  }
  while (he != he0);
}

void QuickHull3D::initBuffers (int nump) {
  for (int i = 0; i < 4; i++) {
    maxVtxs.push_back(NULL);
    minVtxs.push_back(NULL);
    discardedFaces.push_back(NULL);
  }
  newFaces = new FaceList();
  claimed = new VertexList();
  unclaimed = new VertexList();
  if ((int)pointBuffer.size() < nump) {
    std::vector<Vertex*> newBuffer(nump);
    for (size_t i=0; i<pointBuffer.size(); i++) {
      newBuffer[i] = pointBuffer[i]; 
      vertexPointIndices.push_back(-1);
    }
    for (int i=pointBuffer.size(); i<nump; i++) {
      newBuffer[i] = new Vertex(); 
    }
    pointBuffer = newBuffer;
  }
  faces.clear();
  claimed->clear();
  numFaces = 0;
  numPoints = nump;
}

void QuickHull3D::setPoints (std::vector<double> coords) { 
  int nump = coords.size() / 3;
  for (int i=0; i<nump; i++) { 
    Vertex* vtx = pointBuffer[i];
    vtx->pnt = v3d(coords[i*3+0], coords[i*3+1], coords[i*3+2]);
    vtx->index = i;
  }
}

void QuickHull3D::setPoints (std::vector<V3d> pnts) { 
  int nump = pnts.size();
  for (int i=0; i<nump; i++) { 
    Vertex* vtx = pointBuffer[i];
    vtx->pnt = pnts[i];
    vtx->index = i;
  }
}

void QuickHull3D::computeMaxAndMin (void) {
  V3d max;
  V3d min;

  for (int i=0; i<3; i++) {
    maxVtxs[i] = minVtxs[i] = pointBuffer[0]; 
  }
  max = pointBuffer[0]->pnt;
  min = pointBuffer[0]->pnt;

  for (int i=1; i<numPoints; i++) {
    V3d pnt = pointBuffer[i]->pnt;
    if (pnt.x > max.x) {
      max.x = pnt.x;
      maxVtxs[0] = pointBuffer[i];
    } else if (pnt.x < min.x) {
      min.x = pnt.x;
      minVtxs[0] = pointBuffer[i];
    }
    if (pnt.y > max.y) {
      max.y = pnt.y;
      maxVtxs[1] = pointBuffer[i];
    } else if (pnt.y < min.y) {
      min.y = pnt.y;
      minVtxs[1] = pointBuffer[i];
    }
    if (pnt.z > max.z) {
      max.z = pnt.z;
      maxVtxs[2] = pointBuffer[i];
    } else if (pnt.z < min.z) {
      min.z = pnt.z;
      minVtxs[2] = pointBuffer[i];
    }
  }

  // this epsilon formula comes from QuickHull, and I'm
  // not about to quibble.
  charLength = fmax(max.x-min.x, max.y-min.y);
  charLength = fmax(max.z-min.z, charLength);
  if (explicitTolerance == AUTOMATIC_TOLERANCE) {
    tolerance =
      3*DOUBLE_PREC*(fmax(fabs(max.x),fabs(min.x))+
                     fmax(fabs(max.y),fabs(min.y))+
                     fmax(fabs(max.z),fabs(min.z)));
  } else {
    tolerance = explicitTolerance; 
  }
}

// Creates the initial simplex from which the hull will be built.
void QuickHull3D::createInitialSimplex ( void ) {
  double max = 0;
  int imax = 0;

  for (int i=0; i<3; i++) {
    double diff = maxVtxs[i]->pnt[i]-minVtxs[i]->pnt[i];
    if (diff > max) {
      max = diff;
      imax = i;
    }
  }

  if (max <= tolerance) {
    error("Input points appear to be coincident");
  }
  std::vector<Vertex*> vtx(4);
  // set first two vertices to be those with the greatest
  // one dimensional separation

  vtx[0] = maxVtxs[imax];
  vtx[1] = minVtxs[imax];

  // set third vertex to be the vertex farthest from
  // the line between vtx0 and vtx1
  V3d u01;
  V3d diff02;
  V3d nrml;
  V3d xprod;
  double maxSqr = 0;
  u01 = vtx[1]->pnt - vtx[0]->pnt;
  u01.normalize();
  for (int i=0; i<numPoints; i++) {
    diff02 = pointBuffer[i]->pnt - vtx[0]->pnt;
    xprod = cross (u01, diff02);
    double lenSqr = xprod.sqr_magnitude();
    if (lenSqr > maxSqr &&
        pointBuffer[i] != vtx[0] &&  // paranoid
        pointBuffer[i] != vtx[1]) {
      maxSqr = lenSqr; 
      vtx[2] = pointBuffer[i];
      nrml = xprod;
    }
  }
  if (sqrt(maxSqr) <= 100*tolerance) {
    error("Input points appear to be collinear");
  }
  nrml.normalize();

  double maxDist = 0;
  double d0 = dot(vtx[2]->pnt, nrml);
  for (int i=0; i<numPoints; i++) {
    double dist = fabs (dot(pointBuffer[i]->pnt, nrml) - d0);
    if (dist > maxDist &&
        pointBuffer[i] != vtx[0] &&  // paranoid
        pointBuffer[i] != vtx[1] &&
        pointBuffer[i] != vtx[2]) {
      maxDist = dist;
      vtx[3] = pointBuffer[i];
    }
  }
  if (fabs(maxDist) <= 100*tolerance) {
    error("Input points appear to be coplanar"); 
  }

  if (debug) {
    std::cout << "initial vertices:\n" <<
      vtx[0]->index << ": [" << vtx[0]->pnt.x << "," << vtx[0]->pnt.y << "," << vtx[0]->pnt.z << "]\n" <<
      vtx[1]->index << ": [" << vtx[1]->pnt.x << "," << vtx[1]->pnt.y << "," << vtx[1]->pnt.z << "]\n" <<
      vtx[2]->index << ": [" << vtx[2]->pnt.x << "," << vtx[2]->pnt.y << "," << vtx[2]->pnt.z << "]\n" <<
      vtx[3]->index << ": [" << vtx[3]->pnt.x << "," << vtx[3]->pnt.y << "," << vtx[3]->pnt.z << "]\n";
  }

  std::vector<Face*> tris(4);

  if ((dot(vtx[3]->pnt, nrml) - d0) < 0) {
    tris[0] = Face::createTriangle (vtx[0], vtx[1], vtx[2]);
    tris[1] = Face::createTriangle (vtx[3], vtx[1], vtx[0]);
    tris[2] = Face::createTriangle (vtx[3], vtx[2], vtx[1]);
    tris[3] = Face::createTriangle (vtx[3], vtx[0], vtx[2]);
      
    for (int i=0; i<3; i++) {
      int k = (i+1)%3;
      tris[i+1]->getEdge(1)->setOpposite (tris[k+1]->getEdge(0));
      tris[i+1]->getEdge(2)->setOpposite (tris[0]->getEdge(k));
    }
  } else {
    tris[0] = Face::createTriangle (vtx[0], vtx[2], vtx[1]);
    tris[1] = Face::createTriangle (vtx[3], vtx[0], vtx[1]);
    tris[2] = Face::createTriangle (vtx[3], vtx[1], vtx[2]);
    tris[3] = Face::createTriangle (vtx[3], vtx[2], vtx[0]);

    for (int i=0; i<3; i++) {
      int k = (i+1)%3;
      tris[i+1]->getEdge(0)->setOpposite (tris[k+1]->getEdge(1));
      tris[i+1]->getEdge(2)->setOpposite (tris[0]->getEdge((3-i)%3));
    }
  }

  for (int i=0; i<4; i++) {
    faces.push_back(tris[i]); 
  }

  for (int i=0; i<numPoints; i++) {
    Vertex* v = pointBuffer[i];

    if (v == vtx[0] || v == vtx[1] || v == vtx[2] || v == vtx[3]) {
      continue;
    }

    maxDist = tolerance;
    Face* maxFace = NULL;
    for (int k=0; k<4; k++) {
      double dist = tris[k]->distanceToPlane (v->pnt);
      if (dist > maxDist) {
        maxFace = tris[k];
        maxDist = dist;
      }
    }
    if (maxFace != NULL) {
      addPointToFace (v, maxFace);
    }	      
  }
}
  
void QuickHull3D::resolveUnclaimedPoints (FaceList* newFaces) {
  Vertex* vtxNext = unclaimed->first();
  for (Vertex* vtx=vtxNext; vtx!=NULL; vtx=vtxNext) {
    vtxNext = vtx->next;
	      
    double maxDist = tolerance;
    Face* maxFace = NULL;
    for (Face* newFace=newFaces->first(); newFace != NULL; newFace=newFace->next) { 
      if (newFace->mark == VISIBLE) {
        double dist = newFace->distanceToPlane(vtx->pnt);
        if (dist > maxDist) {
          maxDist = dist;
          maxFace = newFace;
        }
        if (maxDist > 1000*tolerance) {
          break;
        }
      }
    }
    if (maxFace != NULL) { 
      addPointToFace (vtx, maxFace);
      if (debug && vtx->index == findIndex) {
        std::cout << findIndex << " CLAIMED BY " << maxFace->getVertexString() << "\n"; 
      }
    } else {
      if (debug && vtx->index == findIndex) {
        std::cout << findIndex << " DISCARDED\n"; 
      } 
    }
  }
}

void QuickHull3D::deleteFacePoints (Face* face, Face* absorbingFace) {
  Vertex* faceVtxs = removeAllPointsFromFace (face);
  if (faceVtxs != NULL) { 
    if (absorbingFace == NULL) {
      unclaimed->addAll (faceVtxs);
    } else {
      Vertex* vtxNext = faceVtxs;
      for (Vertex* vtx=vtxNext; vtx!=NULL; vtx=vtxNext) {
        vtxNext = vtx->next;
        double dist = absorbingFace->distanceToPlane (vtx->pnt);
        if (dist > tolerance) { 
          addPointToFace (vtx, absorbingFace);
        } else { 
          unclaimed->add (vtx);
        }
      }
    }
  }
}

double QuickHull3D::oppFaceDistance (HalfEdge* he) {
  return he->face->distanceToPlane (he->opposite->face->getCentroid());
}

void QuickHull3D::calculateHorizon (V3d eyePnt, HalfEdge* edge0, Face* face, std::vector<HalfEdge*>& horizon) {
  deleteFacePoints (face, NULL);
  face->mark = DELETED;
  if (debug) {
    std::cout << "  visiting face " << face->getVertexString() << "\n";
  }
  HalfEdge* edge;
  if (edge0 == NULL) {
    edge0 = face->getEdge(0);
    edge = edge0;
  } else {
    edge = edge0->getNext();
  }
  do {
    Face* oppFace = edge->oppositeFace();
    if (oppFace->mark == VISIBLE) {
      if (oppFace->distanceToPlane (eyePnt) > tolerance) {
        calculateHorizon (eyePnt, edge->getOpposite(), oppFace, horizon);
      } else {
        horizon.push_back (edge);
        if (debug) {
          std::cout << "  adding horizon edge " << edge->getVertexString() << "\n";
        }
      }
    }
    edge = edge->getNext();
  } while (edge != edge0);
}

void QuickHull3D::addNewFaces (FaceList* newFaces, Vertex* eyeVtx, std::vector<HalfEdge*>& horizon) { 
  newFaces->clear();

  HalfEdge* hedgeSidePrev = NULL;
  HalfEdge* hedgeSideBegin = NULL;

  for (auto horizonHe : horizon) {
    HalfEdge* hedgeSide = addAdjoiningFace (eyeVtx, horizonHe);
    if (debug) {
      std::cout << "new face: " << hedgeSide->face->getVertexString() << "\n";
    }
    if (hedgeSidePrev != NULL) {
      hedgeSide->next->setOpposite (hedgeSidePrev);		 
    } else {
      hedgeSideBegin = hedgeSide; 
    }
    newFaces->add (hedgeSide->getFace());
    hedgeSidePrev = hedgeSide;
  }
  hedgeSideBegin->next->setOpposite (hedgeSidePrev);
}

Vertex* QuickHull3D::nextPointToAdd() {
  if (!claimed->isEmpty()) {
    Face* eyeFace = claimed->first()->face;
    Vertex* eyeVtx = NULL;
    double maxDist = 0;
    for (Vertex* vtx=eyeFace->outside;
         vtx != NULL && vtx->face==eyeFace;
         vtx = vtx->next) {
      double dist = eyeFace->distanceToPlane(vtx->pnt);
      if (dist > maxDist) {
        maxDist = dist;
        eyeVtx = vtx;
      }
    }
    return eyeVtx;
  } else {
    return NULL;
  }
}
	
void QuickHull3D::addPointToHull(Vertex* eyeVtx) {
  horizon.clear();
  unclaimed->clear();
	      
  if (debug) {
    std::cout << "Adding point: " << eyeVtx->index
              << " which is " << eyeVtx->face->distanceToPlane(eyeVtx->pnt) 
              << " above face " << eyeVtx->face->getVertexString() << "\n";
  }
  removePointFromFace (eyeVtx, eyeVtx->face);
  calculateHorizon (eyeVtx->pnt, NULL, eyeVtx->face, horizon);
  newFaces->clear();
  addNewFaces (newFaces, eyeVtx, horizon);
	     
  // first merge pass ... merge faces which are non-convex
  // as determined by the larger face
	     
  for (Face* face = newFaces->first(); face!=NULL; face=face->next) { 
    if (face->mark == VISIBLE) {
      while (doAdjacentMerge(face, NONCONVEX_WRT_LARGER_FACE))
        ;
    }
  }		 
  // second merge pass ... merge faces which are non-convex
  // wrt either face	     
  for (Face* face = newFaces->first(); face!=NULL; face=face->next) { 
    if (face->mark == NON_CONVEX) {
      face->mark = VISIBLE;
      while (doAdjacentMerge(face, NONCONVEX))
        ;
    }
  }	
  resolveUnclaimedPoints(newFaces);
}

void QuickHull3D::buildHull () {
  int cnt = 0;
  Vertex* eyeVtx;

  computeMaxAndMin ();
  createInitialSimplex ();
  while ((eyeVtx = nextPointToAdd()) != NULL) {
    addPointToHull (eyeVtx);
    cnt++;
    if (debug) {
      std::cout << "iteration " << cnt << " done\n"; 
    }
  }
  reindexFacesAndVertices();
  if (debug) {
    std::cout << "hull done\n";
  }
}

void QuickHull3D::reindexFacesAndVertices() { 
  for (int i=0; i<numPoints; i++) {
    pointBuffer[i]->index = -1; 
  }
  // remove inactive faces and mark active vertices
  std::vector<Face*> newFaces;
  for (auto face : faces) {
    if (face->mark == VISIBLE) {
      markFaceVertices (face, 0);
      newFaces.push_back(face);
    }
  }
  numFaces = newFaces.size();
  faces = newFaces;
  // reindex vertices
  numVertices = 0;
  vertexPointIndices.clear();
  for (int i=0; i<numPoints; i++) {
    Vertex* vtx = pointBuffer[i];
    if (vtx->index == 0) {
      vertexPointIndices.push_back(i);
      vtx->index = numVertices++;
    }
  }
}

bool QuickHull3D::checkFaceConvexity (Face* face, double tol, std::ostream& ps) {
  double dist;
  HalfEdge* he = face->he0;
  do {
    face->checkConsistency();
    // make sure edge is convex
    dist = oppFaceDistance (he);
    if (dist > tol) {
      ps << "Edge " << he->getVertexString() << " non-convex by " << dist;
      return false;
    }
    dist = oppFaceDistance (he->opposite);
    if (dist > tol) {
      ps << "Opposite edge " << he->opposite->getVertexString() << " non-convex by " << dist;
      return false;
    }
    if (he->next->oppositeFace() == he->oppositeFace()) {
      ps << "Redundant vertex " << he->head()->index << " in face " << face->getVertexString();
      return false;
    }
    he = he->next;
  } while (he != face->he0);	   
  return true;
}

bool QuickHull3D::checkFaces(double tol, std::ostream& ps) { 
  // check edge convexity
  bool convex = true;
  for (auto face : faces) {
    if (face->mark == VISIBLE) {
      if (!checkFaceConvexity (face, tol, ps)) {
        convex = false;
      }
    }
  }
  return convex;
}

void QuickHull3D::build (std::vector<V3d> points) {
  int nump = points.size();
  if (nump < 4) {
    error ("Less than four input points specified");
  }
  if ((int)points.size() < nump) {
    error ("Point array too small for specified number of points"); 
  }
  initBuffers (nump);
  setPoints (points);
  buildHull ();
}

void QuickHull3D::triangulate () {
  double minArea = 1000*charLength*DOUBLE_PREC;
  newFaces->clear();
  for (auto face : faces) {
    if (face->mark == VISIBLE) { 
      face->triangulate (newFaces, minArea);
      // splitFace (face);
    }
  }
  for (Face* face=newFaces->first(); face!=NULL; face=face->next) {
    faces.push_back(face);
  }
}

std::vector<V3d> QuickHull3D::getVertices(void) {
  std::vector<V3d> vtxs(numVertices);
  for (int i=0; i<numVertices; i++) {
    vtxs[i] = pointBuffer[vertexPointIndices[i]]->pnt;
  }
  return vtxs;
}

std::vector<int> QuickHull3D::getVertexPointIndices(void) {
  std::vector<int> indices;
  for (int i=0; i<numVertices; i++) {
    indices.push_back(vertexPointIndices[i]);
  }
  return indices;
}

std::vector< std::vector< int > > QuickHull3D::getFaces (int indexFlags) {
  std::vector< std::vector< int > > allFaces;
  for (auto face : faces) {
    allFaces.push_back(getFaceIndices ( face, indexFlags) );
  }
  return allFaces;
}

bool QuickHull3D::check (std::ostream& ps) {
  return check (ps, getDistanceTolerance());
}
  
bool QuickHull3D::check (std::ostream& ps, double tol) {
  // check to make sure all edges are fully connected
  // and that the edges are convex
  double dist;
  double pointTol = 10*tol;

  if (!checkFaces(tolerance, ps)) {
    return false; 
  }

  // check point inclusion
  for (int i=0; i<numPoints; i++) {
    V3d pnt = pointBuffer[i]->pnt;
    for (auto face : faces) {
      if (face->mark == VISIBLE) { 
        dist = face->distanceToPlane (pnt);
        if (dist > pointTol) {
          ps << "Point " << i << " " << dist << " above face " << face->getVertexString();
          return false;
        }
      }
    }
  }
  return true;
}

Mesh QuickHull3D::to_mesh (void) {
  triangulate();
  Array<V3d> new_points;
  auto vertices = getVertices();
  for (auto pt : vertices) 
    new_points.append(pt);
  
  Array<V3i> new_faces;
  auto faces = getFaces();
  for (auto face : faces) {
    ensure(face.size() == 3, "MUST BE TRI FACES AFTER HULLIFYING");
    new_faces.append(v3i(face[0], face[1], face[2]));
  }

  return fab_mesh(new_faces, new_points);
}

Mesh quick_hull (Mesh mesh) {
  std::vector< V3d > points;
  for (auto p : mesh.points) 
    points.push_back(p);
  QuickHull3D hull;
  // hull.setDebug(true);
  hull.build(points);
  return hull.to_mesh();
}
