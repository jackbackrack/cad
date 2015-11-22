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

template <class TV>
class HalfEdge;
template <class TV>
class Vertex;
template <class TV>
class FaceList;

static const double DOUBLE_PREC = 2.2204460492503131e-16;
static const double INF = 1e99;

typedef Vector<double,3> V3d;
typedef Vector<double,2> V2d;
typedef Vector<int,3> V3i;

static const int VISIBLE = 1;
static const int NON_CONVEX = 2;
static const int DELETED = 3;

inline V3d v3d (double x, double y, double z) {
  V3d res(x, y, z);
  return res;
}

inline V2d v2d (double x, double y) {
  V2d res(x, y);
  return res;
}

inline V3i v3i (int x, int y, int z) {
  V3i res(x, y, z);
  return res;
}

template<class TV>
class Face {
 private:
  TV normal;
  TV centroid;

  void computeNormalAndCentroid(void);
  void computeNormalAndCentroid(double minArea);
  double areaSquared (HalfEdge<TV>* hedge0, HalfEdge<TV>* hedge1);

 public:

  HalfEdge<TV>* he0;
  double area;
  double planeOffset;
  int mark = VISIBLE;
  int index;
  int numVerts;
  Face<TV>* next;
  Vertex<TV>* outside;

  TV computeCentroid (void);
  TV computeNormal (double minArea);
  TV computeNormal (void);
  static Face<TV>* createTriangle (Vertex<TV>* v0, Vertex<TV>* v1, Vertex<TV>* v2);
  static Face<TV>* createTriangle (Vertex<TV>* v0, Vertex<TV>* v1, Vertex<TV>* v2, double minArea);
  Face (void) : he0(NULL), mark(VISIBLE), next(NULL), outside(NULL) { }
  HalfEdge<TV>* getEdge(int i);
  HalfEdge<TV>* getFirstEdge(void) { return he0; }
  HalfEdge<TV>* findEdge (Vertex<TV>* vt, Vertex<TV>* vh);
  double distanceToPlane (TV p);
  TV getNormal (void) { return normal; }
  TV getCentroid (void) { return centroid; }
  int numVertices(void) { return numVerts; }
  std::string getVertexString (void);
  std::vector< int > getVertexIndices (void);
  Face<TV>* connectHalfEdges (HalfEdge<TV>* hedgePrev, HalfEdge<TV>* hedge);
  void checkConsistency(void);
  std::vector<Face<TV>*> mergeAdjacentFace (HalfEdge<TV>* hedgeAdj);
  void triangulate (FaceList<TV>* newFaces, double minArea);
};
      
template<class TV>
class FaceList {
 private:
  Face<TV>* head;
  Face<TV>* tail;
 public:
  void clear(void) { head = tail = NULL; }
  void add (Face<TV>* vtx);
  inline Face<TV>* first(void) { return head; }
  inline bool isEmpty () { return head == NULL; }
  FaceList (void) : head(NULL), tail(NULL) { }
};

template<class TV>
class HalfEdge {
 public:
  Vertex<TV>* vertex;
  Face<TV>* face;
  HalfEdge<TV>* next;
  HalfEdge<TV>* prev;
  // Half-edge associated with the opposite triangle adjacent to this edge.
  HalfEdge<TV>* opposite;

  // Constructs a HalfEdge with head vertex v and left-hand triangular face f.
  HalfEdge (Vertex<TV>* v, Face<TV>* f) : vertex(v), face(f), next(NULL), prev(NULL) { }
  HalfEdge (void) : vertex(NULL), face(NULL), next(NULL), prev(NULL)  { }

  // Sets the value of the next edge adjacent (counter-clockwise) to
  // this one within the triangle.
  void setNext (HalfEdge<TV>* edge) { next = edge; }
	
  // Gets the value of the next edge adjacent ccw to this one within the triangle.
  HalfEdge<TV>* getNext(void) { return next; }

  // Sets the value of the previous edge adjacent cw to this one within the triangle.
  void setPrev (HalfEdge<TV>* edge) { prev = edge; }
  // Gets the value of the previous edge adjacent cw to this one within the triangle.
  HalfEdge<TV>* getPrev(void) { return prev; }
  // Returns the triangular face located to the left of this half-edge
  Face<TV>* getFace(void) { return face; }
  // Returns the half-edge opposite to this half-edge.
  HalfEdge<TV>* getOpposite(void) { return opposite; }
  void setOpposite (HalfEdge<TV>* edge) {
    opposite = edge;
    edge->opposite = this;
  }
  // Returns the head vertex associated with this half-edge.
  Vertex<TV>* head(void) { return vertex; }
  // Returns the tail vertex associated with this half-edge.
  Vertex<TV>* tail(void) { return prev != NULL ? prev->vertex : NULL; }
  // Returns the opposite triangular face associated with this half-edge.
  Face<TV>* oppositeFace(void) { return opposite != NULL ? opposite->face : NULL; }

  //  Produces a string identifying this half-edge by the point
  //  index values of its tail and head vertices.
  std::string getVertexString(void);

  // Returns the length of this half-edge.
  double length(void);
  double lengthSquared(void);
};

template<class TV>
class Vertex {
 public:
  TV pnt;
  // back index into an array
  int index;
  Vertex<TV>* prev;
  Vertex<TV>* next;
  Face<TV>* face;
 Vertex(TV pt, int idx) : pnt(pt), index(idx) { }
  Vertex (void) : prev(NULL), next(NULL), face(NULL) { }
};

template<class TV>
class VertexList {
 private:
  Vertex<TV>* head;
  Vertex<TV>* tail;

 public:
  VertexList (void) : head(NULL), tail(NULL) { }
  void clear(void) { head = tail = NULL; }

  void add (Vertex<TV>* vtx);
  // Adds a chain of vertices to the end of this list.
  void addAll (Vertex<TV>* vtx);
  // Deletes a vertex from this list.
  void del (Vertex<TV>* vtx);
  // Deletes a chain of vertices from this list.
  void del (Vertex<TV>* vtx1, Vertex<TV>* vtx2);
  // Inserts a vertex into this list before another specificed vertex.
  void insertBefore (Vertex<TV>* vtx, Vertex<TV>* next);
  // Returns the first element in this list.
  Vertex<TV>* first(void) { return head; }

  bool isEmpty(void) { return head == NULL; }
};

static const int CLOCKWISE = 0x1;
static const int POINT_RELATIVE = 0x8;
static const double AUTOMATIC_TOLERANCE = -1;
static const int NONCONVEX_WRT_LARGER_FACE = 1;
static const int NONCONVEX = 2;

template<class TV>
class QuickHull {
 private:
  std::vector<Face<TV>*> discardedFaces;
  std::vector<Vertex<TV>*> maxVtxs;
  std::vector<Vertex<TV>*> minVtxs;
  FaceList<TV>* newFaces;
  VertexList<TV>* unclaimed;
  VertexList<TV>* claimed;

  void addPointToFace (Vertex<TV>* vtx, Face<TV>* face);
  void removePointFromFace (Vertex<TV>* vtx, Face<TV>* face);
  Vertex<TV>* removeAllPointsFromFace (Face<TV>* face);
  HalfEdge<TV>* findHalfEdge (Vertex<TV>* tail, Vertex<TV>* head);
  std::vector< int > getFaceIndices (Face<TV>* face, int flags);
  bool doAdjacentMerge (Face<TV>* face, int mergeType);
  HalfEdge<TV>* addAdjoiningFace (Vertex<TV>* eyeVtx, HalfEdge<TV>* he);
  void markFaceVertices (Face<TV>* face, int mark);

 protected:
  int findIndex = -1;
  bool debug = false;
  double charLength;
  std::vector<Vertex<TV>*> pointBuffer;
  std::vector<int> vertexPointIndices;
  std::vector<Face<TV>*> faces;
  std::vector<HalfEdge<TV>*> horizon;
  int numVertices;
  int numFaces;
  int numPoints;
  double explicitTolerance = AUTOMATIC_TOLERANCE;
  double tolerance;

  void initBuffers (int nump);
  void setPoints (std::vector<double> coords);
  void setPoints (std::vector<TV> pnts);
  void computeMaxAndMin (void);
  // Creates the initial simplex from which the hull will be built.
  void createInitialSimplex ( void );
  void resolveUnclaimedPoints (FaceList<TV>* newFaces);
  void deleteFacePoints (Face<TV>* face, Face<TV>* absorbingFace);
  double oppFaceDistance (HalfEdge<TV>* he);
  void calculateHorizon (TV eyePnt, HalfEdge<TV>* edge0, Face<TV>* face, std::vector<HalfEdge<TV>*>& horizon);
  void addNewFaces (FaceList<TV>* newFaces, Vertex<TV>* eyeVtx, std::vector<HalfEdge<TV>*>& horizon);
  void addPointToHull(Vertex<TV>* eyeVtx);
  void buildHull (void);
  void reindexFacesAndVertices(void);
  bool checkFaceConvexity (Face<TV>* face, double tol, std::ostream& ps);
  bool checkFaces(double tol, std::ostream& ps);
  Vertex<TV>* nextPointToAdd(void);
  
 public:
  bool getDebug() { return debug; }
  void setDebug (bool enable) { debug = enable; }
  double getDistanceTolerance ( void ) { return tolerance; }
  void setExplicitDistanceTolerance(double tol) { explicitTolerance = tol; }
  double getExplicitDistanceTolerance() { return explicitTolerance; }
  void build (std::vector<TV> points);
  QuickHull (void) { }
  QuickHull (std::vector<TV> points) { build(points); }

  // Triangulates any non-triangular hull faces. In some cases, due to
  // precision issues, the resulting triangles may be very thin or small,
  // and hence appear to be non-convex (this same limitation is present
  // in <a href=http://www.qhull.org>qhull</a>).
  void triangulate (void);
  int getNumVertices() { return numVertices; }

  std::vector<TV> getVertices(void);
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
  // @see QuickHull#getVertices()
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
  Nested<TV2> to_poly (void);
};

extern void error(std::string);

template<class TV>
void Face<TV>::computeNormalAndCentroid(void) {
  normal = computeNormal ();
  centroid = computeCentroid ();
  planeOffset = dot(normal, centroid);
  int numv = 0;
  auto he = he0;
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

template<class TV>
void Face<TV>::computeNormalAndCentroid(double minArea) {
  normal = computeNormal (minArea);
  centroid = computeCentroid ();
  planeOffset = dot(normal, centroid);
}

template<class TV>
double Face<TV>::areaSquared (HalfEdge<TV>* hedge0, HalfEdge<TV>* hedge1) {
  // return the squared area of the triangle defined
  // by the half edge hedge0 and the point at the
  // head of hedge1.

  auto p0 = hedge0->tail()->pnt;
  auto p1 = hedge0->head()->pnt;
  auto p2 = hedge1->head()->pnt;

  auto d10 = p1 - p0;
  auto d20 = p2 - p0;
  auto c   = cross(d10, d20);
  return dot(c, c);
}

template<class TV>
TV Face<TV>::computeCentroid (void) {
  TV centroid;
  auto he = he0;
  do {
    centroid += he->head()->pnt;
    he = he->next;
  } while (he != he0);
  return 1/(double)numVerts * centroid;
}

template<class TV>
TV Face<TV>::computeNormal (double minArea) {
  auto normal = computeNormal();
  if (area < minArea) {
    // make the normal more robust by removing
    // components parallel to the longest edge
    HalfEdge<TV>* hedgeMax = NULL;
    double lenSqrMax = 0;
    HalfEdge<TV>* hedge = he0;
    do {
      double lenSqr = hedge->lengthSquared();
      if (lenSqr > lenSqrMax) {
        hedgeMax = hedge;
        lenSqrMax = lenSqr;
      }
      hedge = hedge->next;
    } while (hedge != he0);

    TV p2 = hedgeMax->head()->pnt;
    TV p1 = hedgeMax->tail()->pnt;
    double lenMax = sqrt(lenSqrMax);
    auto u = (p2 - p1) / lenMax;
    auto dotty = dot(normal, u);
    normal -= dotty * u;
    
    normal.normalize();	      
  }
  return normal;
}

template<class TV>
TV Face<TV>::computeNormal (void) {
  auto he1 = he0->next;
  auto he2 = he1->next;

  auto p0 = he0->head()->pnt;
  auto p2 = he1->head()->pnt;

  auto d2 = p2 - p0;

  TV normal;

  numVerts = 2;

  while (he2 != he0) { 
    auto d1 = d2;
    p2 = he2->head()->pnt;
    d2 = p2 - p0;
    normal += cross(d1, d2);
    he1 = he2;
    he2 = he2->next;
    numVerts++;
  }
  area = normal.magnitude();
  return (1/area) * normal;
}

template<class TV>
Face<TV>* Face<TV>::createTriangle (Vertex<TV>* v0, Vertex<TV>* v1, Vertex<TV>* v2) {
  return createTriangle (v0, v1, v2, 0);
}

template<class TV>
Face<TV>* Face<TV>::createTriangle (Vertex<TV>* v0, Vertex<TV>* v1, Vertex<TV>* v2, double minArea) {
  auto face = new Face<TV>();
  auto he0 = new HalfEdge<TV> (v0, face);
  auto he1 = new HalfEdge<TV> (v1, face);
  auto he2 = new HalfEdge<TV> (v2, face);

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

template<class TV>
HalfEdge<TV>* Face<TV>::getEdge (int i) {
  auto he = he0;
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

template<class TV>
// finds the half-edge within this face which has tail vt and head vh
HalfEdge<TV>* Face<TV>::findEdge (Vertex<TV>* vt, Vertex<TV>* vh) {
  auto he = he0;
  do {
    if (he->head() == vh && he->tail() == vt) {
      return he;
    }
    he = he->next;
  } while (he != he0);
  return NULL;
}

template<class TV>
double Face<TV>::distanceToPlane (TV p) {
  return dot(normal, p) - planeOffset;
}

template<class TV>
std::string Face<TV>::getVertexString (void) {
  std::stringstream ss;
  auto he = he0;
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

template<class TV>
std::vector<int> Face<TV>::getVertexIndices (void) {
  auto he = he0;
  std::vector<int> indices;
  do {
    indices.push_back(he->head()->index);
    he = he->next;
  } while (he != he0);
  return indices;
}

template<class TV>
Face<TV>* Face<TV>::connectHalfEdges (HalfEdge<TV>* hedgePrev, HalfEdge<TV>* hedge) {
  Face<TV>* discardedFace = NULL;
  if (hedgePrev->oppositeFace() == hedge->oppositeFace()) { 
    // then there is a redundant edge that we can get rid off
    auto oppFace = hedge->oppositeFace();
    HalfEdge<TV>* hedgeOpp;

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

template<class TV>
void Face<TV>::checkConsistency(void) {
  // do a sanity check on the face
  auto hedge = he0; 
  double maxd = 0;
  int numv = 0;

  if (numVerts < 3) {
    std::stringstream ss;
    ss << "degenerate face: " << getVertexString();
    error(ss.str());
  }
  do {
    auto hedgeOpp = hedge->getOpposite();
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
    auto oppFace = hedgeOpp->face;
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

template<class TV>
std::vector<Face<TV>*> Face<TV>::mergeAdjacentFace (HalfEdge<TV>* hedgeAdj) {
  auto oppFace = hedgeAdj->oppositeFace();
  std::vector<Face*> discarded;

  discarded.push_back(oppFace);
  oppFace->mark = DELETED;

  auto hedgeOpp = hedgeAdj->getOpposite();

  auto hedgeAdjPrev = hedgeAdj->prev;
  auto hedgeAdjNext = hedgeAdj->next;
  auto hedgeOppPrev = hedgeOpp->prev;
  auto hedgeOppNext = hedgeOpp->next;

  while (hedgeAdjPrev->oppositeFace() == oppFace) {
    hedgeAdjPrev = hedgeAdjPrev->prev;
    hedgeOppNext = hedgeOppNext->next;
  }
	   
  while (hedgeAdjNext->oppositeFace() == oppFace) {
    hedgeOppPrev = hedgeOppPrev->prev;
    hedgeAdjNext = hedgeAdjNext->next;
  }

  HalfEdge<TV>* hedge;

  for (hedge=hedgeOppNext; hedge!=hedgeOppPrev->next; hedge=hedge->next) {
    hedge->face = this;
  }

  if (hedgeAdj == he0) {
    he0 = hedgeAdjNext; 
  }
	   
  // handle the half edges at the head
  Face<TV>* discardedFace;

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

template<class TV>
void Face<TV>::triangulate (FaceList<TV>* newFaces, double minArea) {
  HalfEdge<TV>* hedge;

  if (numVertices() < 4) {
    return; 
  }

  auto v0 = he0->head();

  hedge = he0->next;
  auto oppPrev = hedge->opposite;
  Face* face0 = NULL;

  for (hedge=hedge->next; hedge!=he0->prev; hedge=hedge->next) {
    auto face = createTriangle (v0, hedge->prev->head(), hedge->head(), minArea);
    face->he0->next->setOpposite (oppPrev);
    face->he0->prev->setOpposite (hedge->opposite);
    oppPrev = face->he0;
    newFaces->add (face);
    if (face0 == NULL) {
      face0 = face; 
    }
  }
  hedge = new HalfEdge<TV> (he0->prev->prev->head(), this);
  hedge->setOpposite (oppPrev);

  hedge->prev = he0;
  hedge->prev->next = hedge;

  hedge->next = he0->prev;
  hedge->next->prev = hedge;

  computeNormalAndCentroid (minArea);
  checkConsistency();

  for (auto face=face0; face!=NULL; face=face->next) {
    face->checkConsistency(); 
  }
}
      
template<class TV>
void FaceList<TV>::add (Face<TV>* vtx) {
  if (head == NULL) {
    head = vtx;
  } else {
    tail->next = vtx; 
  }
  vtx->next = NULL;
  tail = vtx;
}


template<class TV>
std::string HalfEdge<TV>::getVertexString(void) {
  std::stringstream ss;
  if (tail() != NULL) {
    ss << tail()->index << "-" << head()->index;
  } else {
    ss << "?-" << head()->index;
  }
  return ss.str();
}

// Returns the length of this half-edge.
template<class TV>
double HalfEdge<TV>::length(void) {
  if (tail() != NULL) {
    return (head()->pnt - tail()->pnt).magnitude();
  } else {
    return -1; 
  }
}

template<class TV>
double HalfEdge<TV>::lengthSquared(void) {
  if (tail() != NULL) {
    return (head()->pnt - tail()->pnt).sqr_magnitude();
  } else {
    return -1; 
  }
}

template<class TV>
void VertexList<TV>::add (Vertex<TV>* vtx) {  
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
template<class TV>
void VertexList<TV>::addAll (Vertex<TV>* vtx) { 
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
template<class TV>
void VertexList<TV>::del (Vertex<TV>* vtx) {
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
template<class TV>
void VertexList<TV>::del (Vertex<TV>* vtx1, Vertex<TV>* vtx2) {
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
template<class TV>
void VertexList<TV>::insertBefore (Vertex<TV>* vtx, Vertex<TV>* next) {
  vtx->prev = next->prev;
  if (next->prev == NULL) {
    head = vtx;
  } else {
    next->prev->next = vtx; 
  }
  vtx->next = next;
  next->prev = vtx;
}

template<class TV>
void QuickHull<TV>::addPointToFace (Vertex<TV>* vtx, Face<TV>* face) {
  vtx->face = face;

  if (face->outside == NULL) {
    claimed->add (vtx);
  } else {
    claimed->insertBefore (vtx, face->outside); 
  }
  face->outside = vtx;
}

template<class TV>
void QuickHull<TV>::removePointFromFace (Vertex<TV>* vtx, Face<TV>* face) {
  if (vtx == face->outside) {
    if (vtx->next != NULL && vtx->next->face == face) {
      face->outside = vtx->next;
    } else {
      face->outside = NULL; 
    }
  }
  claimed->del (vtx);
}

template<class TV>
Vertex<TV>* QuickHull<TV>::removeAllPointsFromFace (Face<TV>* face) {
  if (face->outside != NULL) { 
    auto end = face->outside;
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

template<class TV>
HalfEdge<TV>* QuickHull<TV>::findHalfEdge (Vertex<TV>* tail, Vertex<TV>* head) { 
  for (auto face : faces) {
    auto he = face->findEdge (tail, head);
    if (he != NULL) {
      return he; 
    }
  }
  return NULL;
}
  
template<class TV>
std::vector<int> QuickHull<TV>::getFaceIndices (Face<TV>* face, int flags) { 
  std::vector<int> indices;
  bool ccw = ((flags & CLOCKWISE) == 0);
  bool pointRelative = ((flags & POINT_RELATIVE) != 0);

  auto hedge = face->he0;
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

template<class TV>
bool QuickHull<TV>::doAdjacentMerge (Face<TV>* face, int mergeType) {
  auto hedge = face->he0;

  bool convex = true;
  do {
    auto oppFace = hedge->oppositeFace();
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

template<class TV>
HalfEdge<TV>* QuickHull<TV>::addAdjoiningFace (Vertex<TV>* eyeVtx, HalfEdge<TV>* he) { 
  auto face = Face<TV>::createTriangle (eyeVtx, he->tail(), he->head());
  faces.push_back (face);
  face->getEdge(-1)->setOpposite(he->getOpposite());
  return face->getEdge(0);
}

template<class TV>
void QuickHull<TV>::markFaceVertices (Face<TV>* face, int mark) {
  auto he0 = face->getFirstEdge();
  auto he = he0;
  do {
    he->head()->index = mark;
    he = he->next;
  }
  while (he != he0);
}

template<class TV>
void QuickHull<TV>::initBuffers (int nump) {
  for (int i = 0; i < 4; i++) {
    maxVtxs.push_back(NULL);
    minVtxs.push_back(NULL);
    discardedFaces.push_back(NULL);
  }
  newFaces = new FaceList<TV>();
  claimed = new VertexList<TV>();
  unclaimed = new VertexList<TV>();
  if ((int)pointBuffer.size() < nump) {
    std::vector<Vertex<TV>*> newBuffer(nump);
    for (size_t i=0; i<pointBuffer.size(); i++) {
      newBuffer[i] = pointBuffer[i]; 
      vertexPointIndices.push_back(-1);
    }
    for (int i=pointBuffer.size(); i<nump; i++) {
      newBuffer[i] = new Vertex<TV>(); 
    }
    pointBuffer = newBuffer;
  }
  faces.clear();
  claimed->clear();
  numFaces = 0;
  numPoints = nump;
}

template<class TV>
void QuickHull<TV>::setPoints (std::vector<double> coords) { 
  int nump = coords.size() / 3;
  for (int i=0; i<nump; i++) { 
    auto vtx = pointBuffer[i];
    vtx->pnt = v3d(coords[i*3+0], coords[i*3+1], coords[i*3+2]);
    vtx->index = i;
  }
}

template<class TV>
void QuickHull<TV>::setPoints (std::vector<TV> pnts) { 
  int nump = pnts.size();
  for (int i=0; i<nump; i++) { 
    auto vtx = pointBuffer[i];
    vtx->pnt = pnts[i];
    vtx->index = i;
  }
}

template<class TV>
void QuickHull<TV>::computeMaxAndMin (void) {
  TV max;
  TV min;

  for (int i=0; i<max.dimension; i++) {
    maxVtxs[i] = minVtxs[i] = pointBuffer[0]; 
  }
  max = pointBuffer[0]->pnt;
  min = pointBuffer[0]->pnt;

  for (int i=1; i<numPoints; i++) {
    auto pnt = pointBuffer[i]->pnt;
    for (int d=0; d<max.dimension; d++) {
      if (pnt[d] > max[d]) {
        max[d] = pnt[d];
        maxVtxs[d] = pointBuffer[i];
      } else if (pnt[d] < min[d]) {
        min[d] = pnt[d];
        minVtxs[d] = pointBuffer[i];
      }
    }
  }

  // this epsilon formula comes from QuickHull, and I'm
  // not about to quibble.
  charLength = -INF;
  for (int d=0; d<max.dimension; d++)
    charLength = fmax(charLength, max[d]-min[d]);
  if (explicitTolerance == AUTOMATIC_TOLERANCE) {
    T tmax = 0.0;
    for (int d=0; d<max.dimension; d++)
      tmax += fmax(fabs(max[d]), fabs(min[d]));
    tolerance = max.dimension*DOUBLE_PREC*tmax;
  } else {
    tolerance = explicitTolerance; 
  }
}

// Creates the initial simplex from which the hull will be built.
template<class TV>
void QuickHull<TV>::createInitialSimplex ( void ) {
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
  std::vector<Vertex<TV>*> vtx(4);
  // set first two vertices to be those with the greatest
  // one dimensional separation

  vtx[0] = maxVtxs[imax];
  vtx[1] = minVtxs[imax];

  // set third vertex to be the vertex farthest from
  // the line between vtx0 and vtx1
  TV u01;
  TV diff02;
  TV nrml;
  TV xprod;
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
      vtx[0]->index << ": " << vtx[0]->pnt << "\n" <<
      vtx[1]->index << ": " << vtx[1]->pnt << "\n" <<
      vtx[2]->index << ": " << vtx[2]->pnt << "\n" << 
      vtx[3]->index << ": " << vtx[3]->pnt << "\n";
  }

  std::vector<Face<TV>*> tris(4);

  if ((dot(vtx[3]->pnt, nrml) - d0) < 0) {
    tris[0] = Face<TV>::createTriangle (vtx[0], vtx[1], vtx[2]);
    tris[1] = Face<TV>::createTriangle (vtx[3], vtx[1], vtx[0]);
    tris[2] = Face<TV>::createTriangle (vtx[3], vtx[2], vtx[1]);
    tris[3] = Face<TV>::createTriangle (vtx[3], vtx[0], vtx[2]);
      
    for (int i=0; i<3; i++) {
      int k = (i+1)%3;
      tris[i+1]->getEdge(1)->setOpposite (tris[k+1]->getEdge(0));
      tris[i+1]->getEdge(2)->setOpposite (tris[0]->getEdge(k));
    }
  } else {
    tris[0] = Face<TV>::createTriangle (vtx[0], vtx[2], vtx[1]);
    tris[1] = Face<TV>::createTriangle (vtx[3], vtx[0], vtx[1]);
    tris[2] = Face<TV>::createTriangle (vtx[3], vtx[1], vtx[2]);
    tris[3] = Face<TV>::createTriangle (vtx[3], vtx[2], vtx[0]);

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
    auto v = pointBuffer[i];

    if (v == vtx[0] || v == vtx[1] || v == vtx[2] || v == vtx[3]) {
      continue;
    }

    maxDist = tolerance;
    Face<TV>* maxFace = NULL;
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
  
template<class TV>
void QuickHull<TV>::resolveUnclaimedPoints (FaceList<TV>* newFaces) {
  auto vtxNext = unclaimed->first();
  for (auto vtx=vtxNext; vtx!=NULL; vtx=vtxNext) {
    vtxNext = vtx->next;
	      
    double maxDist = tolerance;
    Face<TV>* maxFace = NULL;
    for (auto newFace=newFaces->first(); newFace != NULL; newFace=newFace->next) { 
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

template<class TV>
void QuickHull<TV>::deleteFacePoints (Face<TV>* face, Face<TV>* absorbingFace) {
  auto faceVtxs = removeAllPointsFromFace (face);
  if (faceVtxs != NULL) { 
    if (absorbingFace == NULL) {
      unclaimed->addAll (faceVtxs);
    } else {
      auto vtxNext = faceVtxs;
      for (auto vtx=vtxNext; vtx!=NULL; vtx=vtxNext) {
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

template<class TV>
double QuickHull<TV>::oppFaceDistance (HalfEdge<TV>* he) {
  return he->face->distanceToPlane (he->opposite->face->getCentroid());
}

template<class TV>
void QuickHull<TV>::calculateHorizon (TV eyePnt, HalfEdge<TV>* edge0, Face<TV>* face, std::vector<HalfEdge<TV>*>& horizon) {
  deleteFacePoints (face, NULL);
  face->mark = DELETED;
  if (debug) {
    std::cout << "  visiting face " << face->getVertexString() << "\n";
  }
  HalfEdge<TV>* edge;
  if (edge0 == NULL) {
    edge0 = face->getEdge(0);
    edge = edge0;
  } else {
    edge = edge0->getNext();
  }
  do {
    auto oppFace = edge->oppositeFace();
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

template<class TV>
void QuickHull<TV>::addNewFaces (FaceList<TV>* newFaces, Vertex<TV>* eyeVtx, std::vector<HalfEdge<TV>*>& horizon) { 
  newFaces->clear();

  HalfEdge<TV>* hedgeSidePrev = NULL;
  HalfEdge<TV>* hedgeSideBegin = NULL;

  for (auto horizonHe : horizon) {
    auto hedgeSide = addAdjoiningFace (eyeVtx, horizonHe);
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

template<class TV>
Vertex<TV>* QuickHull<TV>::nextPointToAdd() {
  if (!claimed->isEmpty()) {
    auto eyeFace = claimed->first()->face;
    Vertex<TV>* eyeVtx = NULL;
    double maxDist = 0;
    for (auto vtx=eyeFace->outside;
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
	
template<class TV>
void QuickHull<TV>::addPointToHull(Vertex<TV>* eyeVtx) {
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
	     
  for (auto face = newFaces->first(); face!=NULL; face=face->next) { 
    if (face->mark == VISIBLE) {
      while (doAdjacentMerge(face, NONCONVEX_WRT_LARGER_FACE))
        ;
    }
  }		 
  // second merge pass ... merge faces which are non-convex
  // wrt either face	     
  for (auto face = newFaces->first(); face!=NULL; face=face->next) { 
    if (face->mark == NON_CONVEX) {
      face->mark = VISIBLE;
      while (doAdjacentMerge(face, NONCONVEX))
        ;
    }
  }	
  resolveUnclaimedPoints(newFaces);
}

template<class TV>
void QuickHull<TV>::buildHull () {
  int cnt = 0;
  Vertex<TV>* eyeVtx;

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

template<class TV>
void QuickHull<TV>::reindexFacesAndVertices() { 
  for (int i=0; i<numPoints; i++) {
    pointBuffer[i]->index = -1; 
  }
  // remove inactive faces and mark active vertices
  std::vector<Face<TV>*> newFaces;
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
    auto vtx = pointBuffer[i];
    if (vtx->index == 0) {
      vertexPointIndices.push_back(i);
      vtx->index = numVertices++;
    }
  }
}

template<class TV>
bool QuickHull<TV>::checkFaceConvexity (Face<TV>* face, double tol, std::ostream& ps) {
  double dist;
  auto he = face->he0;
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

template<class TV>
bool QuickHull<TV>::checkFaces(double tol, std::ostream& ps) { 
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

template<class TV>
void QuickHull<TV>::build (std::vector<TV> points) {
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

template<class TV>
void QuickHull<TV>::triangulate () {
  double minArea = 1000*charLength*DOUBLE_PREC;
  newFaces->clear();
  for (auto face : faces) {
    if (face->mark == VISIBLE) { 
      face->triangulate (newFaces, minArea);
      // splitFace (face);
    }
  }
  for (auto face=newFaces->first(); face!=NULL; face=face->next) {
    faces.push_back(face);
  }
}

template<class TV>
std::vector<TV> QuickHull<TV>::getVertices(void) {
  std::vector<TV> vtxs(numVertices);
  for (int i=0; i<numVertices; i++) {
    vtxs[i] = pointBuffer[vertexPointIndices[i]]->pnt;
  }
  return vtxs;
}

template<class TV>
std::vector<int> QuickHull<TV>::getVertexPointIndices(void) {
  std::vector<int> indices;
  for (int i=0; i<numVertices; i++) {
    indices.push_back(vertexPointIndices[i]);
  }
  return indices;
}

template<class TV>
std::vector< std::vector< int > > QuickHull<TV>::getFaces (int indexFlags) {
  std::vector< std::vector< int > > allFaces;
  for (auto face : faces) {
    allFaces.push_back(getFaceIndices ( face, indexFlags) );
  }
  return allFaces;
}

template<class TV>
bool QuickHull<TV>::check (std::ostream& ps) {
  return check (ps, getDistanceTolerance());
}
  
template<class TV>
bool QuickHull<TV>::check (std::ostream& ps, double tol) {
  // check to make sure all edges are fully connected
  // and that the edges are convex
  double dist;
  double pointTol = 10*tol;

  if (!checkFaces(tolerance, ps)) {
    return false; 
  }

  // check point inclusion
  for (int i=0; i<numPoints; i++) {
    TV pnt = pointBuffer[i]->pnt;
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

template<class TV>
Mesh QuickHull<TV>::to_mesh (void) {
  triangulate();
  Array<TV> new_points;
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

Mesh quick_hull_mesh (Mesh mesh) {
  std::vector< V3d > points;
  for (auto p : mesh.points) 
    points.push_back(p);
  QuickHull<V3d> hull;
  // hull.setDebug(true);
  hull.build(points);
  return hull.to_mesh();
}

template<class TV>
Nested<TV2> QuickHull<TV>::to_poly (void) {
  auto mesh = to_mesh();
  auto normals = mesh.soup->element_normals(RawArray<TV3>(mesh.points));
  Array<IV3> bot_triangles;

  // find bottom facing triangles
  for(int i=0; i < normals.size(); i++) {
    if (normals[i].z < 0) 
      bot_triangles.append(mesh.soup->elements[i]);
  }

  // create mesh of only bottom facing triangles
  auto bot_mesh = fab_mesh(bot_triangles, mesh.points);
  auto bot_topo = new_<TriangleTopology>(bot_mesh.soup->elements);
  // find boundary contour
  auto loops = bot_topo->boundary_loops();

  Nested<TV2,false> poly;
  Array<TV2> contour;
  ensure(loops.size() == 1, "ONLY ONE BOUNDARY LOOP IN BOTTOM MESH OF CONVEX HULL");
  // project into plane
  for (auto e : loops[0]) {
    auto ends = bot_topo->vertices(e);
    auto pt = mesh.points[(int)ends[0]];
    contour.append(vec(pt.x, pt.y));
  }
  poly.append(contour);
  poly.freeze();
  
  return poly;
}

Nested<TV2> quick_hull_poly (Nested<V2d> poly) {
  std::vector< V3d > points;
  double t = 0.0;
  // Assign unique z coordinates
  for (auto c : poly) {
    for (auto p : c) {
      // auto v = v3d(p.x, p.y, dot(p, p));
      auto v = v3d(p.x, p.y, t);
      // printf("V %f,%f,%f\n", v.x, v.y, v.z);
      points.push_back(v);
      t += 1.0;
    }
  }
  QuickHull<V3d> hull;
  // hull.setDebug(true);
  hull.build(points);
  return hull.to_poly();
}
