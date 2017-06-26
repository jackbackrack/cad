#ifndef __CAD__
#define __CAD__

#include <geode/mesh/io.h>
#include <geode/vector/Matrix3x3.h>
#include <geode/vector/Matrix4x4.h>
#include <geode/vector/Matrix4x4.h>
#include <geode/vector/Rotation.h>
#include <geode/geometry/platonic.h>
#include <geode/geometry/Triangle3d.h>
#include <geode/geometry/polygon.h>
#include <geode/geometry/offset_mesh.h>
#include <geode/exact/mesh_csg.h>
// #include <geode/exact/delaunay.h>
#include <geode/exact/polygon_csg.h>
#include <geode/mesh/SegmentSoup.h>
#include <geode/mesh/TriangleTopology.h>

#include <cstdio>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

#ifdef MACOSX
#include <GLUT/glut.h> 
#else
#include <GL/glut.h> // Linux, Windows
#endif

using namespace geode;

extern void error (std::string msg);
extern void error (std::string msg, std::string arg);
extern void ensure (bool val, std::string msg);

typedef real T;
typedef Vector<T,2> TV2;
typedef Vector<T,3> TV3;
typedef Vector<int,2> IV2;
typedef Vector<int,3> IV3;
typedef Vector<int,4> IV4;
// typedef Tuple<Ref<const TriangleSoup>, Array<TV3>> Mesh;

// cons compatibility
typedef Triangle<TV3> Tri;
typedef T flo_t;
typedef TV3 V3d;
typedef TV3 Vec;

inline Tri tri (TV3 p0, TV3 p1, TV3 p2) { return Triangle<TV3>(p0, p1, p2); }

const T INFTY = (T)INFINITY;

extern Vec rnd_vec_of (const Vec &lo, const Vec &hi);
extern Vec rnd_vec_of (const flo_t lo, const flo_t hi);

class Interval {
 public:
  flo_t lo;
  flo_t hi;
  inline bool operator <  (Interval o) { return hi < o.lo; }
  inline bool operator <= (Interval o) { return hi <= o.lo; }
  inline bool operator >  (Interval o) { return lo > o.hi; }
  inline bool operator >= (Interval o) { return lo >= o.hi; }
 Interval(void) : lo(0.0), hi(0.0) /*, Geom(geom_interval_kind) */ { }
 Interval(flo_t lo, flo_t hi) : lo(lo), hi(hi) /*, Geom(geom_interval_kind) */ { }
};

inline Interval interval(flo_t lo, flo_t hi) {
  const Interval res(lo, hi); return res;
}

class Boxy {
 public:
  Vec lo;
  Vec hi;
  inline Vec center( void ) const { return (lo + hi)*0.5; }
  inline Vec rnd_vec( void ) const { return rnd_vec_of(lo, hi); }
  inline Vec dims (void) const { return vec(hi.x - lo.x, hi.y - lo.y, hi.z - lo.z); }
  inline bool is_inside (const Vec &p) const {
    return lo.x < p.x && lo.y < p.y && lo.z < p.z &&
           p.x < hi.x && p.y < hi.y && p.z < hi.z;
  }
  inline bool is_overlap (const Boxy &b) const { 
    for (int i = 0; i < 3; i++)
      if (b.hi[i] < lo[i] || b.lo[i] > hi[i])
        return false;
    return true;
  }
 Boxy(Vec lo, Vec hi) : lo(lo), hi(hi) /*, Geom(geom_box_kind) */ { }
};

inline Boxy box(Vec lo, Vec hi)  {
  const Boxy res(lo, hi); return res;
}
inline Boxy add(Boxy a, Boxy b)  {
  return box(vec(min(a.lo.x, b.lo.x), min(a.lo.y, b.lo.y), min(a.lo.z, b.lo.z)),
             vec(max(a.hi.x, b.hi.x), max(a.hi.y, b.hi.y), max(a.hi.z, b.hi.z)));
}
inline Boxy add(Boxy b, Vec p)  {
  return box(vec(min(b.lo.x, p.x), min(b.lo.y, p.y), min(b.lo.z, p.z)),
             vec(max(b.hi.x, p.x), max(b.hi.y, p.y), max(b.hi.z, p.z)));
}

inline Segment<TV3> segment(TV3 from, TV3 to) {
  Segment<TV3> res(from, to); return res;
}
inline Segment<TV2> segment(TV2 from, TV2 to) {
  Segment<TV2> res(from, to); return res;
}

inline TV3 get_color(TV3 n) {
  return vec((n.x > 0.0 ? n.x : 0.0) + (n.y < 0.0 ? -0.5*n.y : 0.0) + (n.z < 0.0 ? -0.5*n.z : 0.0),
             (n.y > 0.0 ? n.y : 0.0) + (n.z < 0.0 ? -0.5*n.z : 0.0) + (n.x < 0.0 ? -0.5*n.x : 0.0),
             (n.z > 0.0 ? n.z : 0.0) + (n.x < 0.0 ? -0.5*n.x : 0.0) + (n.y < 0.0 ? -0.5*n.y : 0.0));
}

// std::string to_str(TV3 v) { return "vec(" + std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z) + ")"; }
// std::string to_str(Segment<TV3> l) { return "LineSegment(" + to_str(l.x0) + "," + to_str(l.x1) + ")"; }
// std::string to_str(Line<TV3> l) { return "Line(" + to_str(l.pos) + "," + to_str(l.dir) + ")"; }
// std::string to_str(Triangle<TV3> t) { return "tri(" + to_str(t.p0) + "," + to_str(t.p1) + "," + to_str(t.p2) ")"; }

struct Mesh {
public:
  Ref<const TriangleSoup> soup;
  Array<TV3> points;
  Mesh(Ref<const TriangleSoup> soup, Array<TV3> points) : soup(soup), points(points) { }
};

extern Mesh fab_mesh (Array<IV3> faces, Array<TV3> points);

extern Matrix<T,4> to_matrix44(Matrix<T,3> M);
extern Matrix<T,4> rotation_matrix(TV3 angles);
extern Matrix<T,4> rotation_matrix(TV3 from, TV3 to);
extern Matrix<T,4> rotation_matrix(TV2 from, TV2 to);
extern Matrix<T,4> scale_matrix(TV3 v);
extern Matrix<T,4> translation_matrix(TV3 v);
extern Matrix<T,4> translation_matrix(TV2 v);

template<class ET>
Nested<ET> line_to_polyline(Array<ET> line) {
  Nested<ET,false> polyline;
  polyline.append(line);
  polyline.freeze();
  return polyline;
}

extern double is_clockwise (Array<TV2> contour);

extern Ref<const TriangleSoup> const_soup(Ref<TriangleSoup> val);
extern Mesh const_mesh(Tuple<Ref<TriangleSoup>, Array<TV3>> val);

extern Mesh simplify_mesh(Mesh mesh);
extern Mesh cleanup_mesh(Mesh mesh);

extern Mesh invert_mesh(Mesh mesh);

extern Array<TV2> invert_contour(Array<TV2> contour);
extern Nested<TV2> invert_poly(Nested<TV2> poly);

extern Array<TV2> cleanup_contour(RawArray<TV2> contour);
extern Nested<TV2> cleanup_poly(Nested<TV2> poly);

extern Nested<TV2> union_add(Nested<TV2> c0, Nested<TV2> c1);
extern Nested<TV2> union_all(Nested<TV2> cs);
extern Nested<TV2> intersection(Nested<TV2> c0, Nested<TV2> c1);
extern Nested<TV2> difference(Nested<TV2> c0, Nested<TV2> c1);

extern Nested<TV2> offset(T a, Nested<TV2> c);

extern Mesh union_add(Mesh mesh0, Mesh mesh1, bool is_cleanup = true);

extern void pretty_print_num(T pt);
extern void pretty_print_v2d(TV2 pt);
extern void pretty_print_v3d(TV3 pt);
extern void pretty_print_v3i(IV3 pt);
extern void pretty_print_array_v2d(Array<TV2> line);
extern void pretty_print_array_v3d(Array<TV3> line);
extern void pretty_print_array_v3i(Array<IV3> line);
extern void pretty_print_matrix(Matrix<T,4> M);
extern void pretty_print_nested_v3d(Nested<TV3> polyline);
extern void pretty_print_nested_v2d(Nested<TV2> polyline);
extern void pretty_print_poly(Nested<TV2> poly);
extern void pretty_print_mesh(Mesh soup);
extern void obj_mesh(Mesh soup);

extern std::string num_to_str (T num);
extern std::string v2d_to_str (TV2 pt);
extern std::string v3d_to_str (TV3 pt);
extern std::string v3i_to_str (IV3 pt);
extern std::string array_v2d_to_str(Array<TV2> line);
extern std::string array_v3d_to_str(Array<TV3> line);
extern std::string array_v3i_to_str(Array<const IV3> line);
extern std::string mesh_to_str(Mesh mesh);
extern std::string array_v2d_to_str(Array<TV2> contour);
extern std::string poly_to_str(Nested<TV2> poly);
extern std::string nested_v3d_to_str(Nested<TV3> polyline);
extern std::string nested_v2d_to_str(Nested<TV2> polyline);
extern std::string matrix_to_str(Matrix<T,4> M);

extern void save_polygon(std::string filename, Nested<TV2> polygon);
extern void save_polyline(std::string filename, Nested<TV2> polyline);

extern Mesh intersection(Mesh mesh0, Mesh mesh1, bool is_cleanup = true);

extern Mesh mesh_from(int start, Mesh soup);

extern Nested<TV2> slice(T z, Mesh soup);

extern Mesh difference(Mesh mesh0, Mesh mesh1, bool is_cleanup = true);

template<class ET> Nested<ET> array_to_nested(Array<ET> contour);

template<class ET> Array<ET> nested_elt(Nested<ET> poly, int idx);

extern Mesh all_mesh(void);

extern Mesh none_mesh(void);

extern Nested<TV2> square_poly(TV2 min, TV2 max);

extern Nested<TV2> square_poly(T rad);

extern Nested<TV2> all_poly(void);

extern Nested<TV2> none_poly(void);

extern Nested<TV2> circle_poly(T rad, int n);

extern Nested<TV2> star_poly(T rad_min, T rad_max, int n);


extern Mesh triangulate (Nested<TV3> poly);
extern Mesh triangulate (Nested<TV2> poly);

extern TV3 mul(Matrix<T,4> m, TV3 pt);

extern TV2 mul(Matrix<T,4> m, TV2 pt);

template<class E> Array<E> mul(Matrix<T,4> m, Array<E> pts) {
  Array<E> res;
  for (auto p : pts)
    res.append(mul(m, p));
  return res;
}
template<class E> Nested<E> mul(Matrix<T,4> m, Nested<E> poly) {
  Nested<E,false> pres;
  for (auto contour : poly) {
    Array<E> cres;
    for (auto p : contour)
      cres.append(mul(m, p));
    pres.append(cres);
  }
  pres.freeze();
  return pres;
}

extern Nested<TV2> mul_poly(Matrix<T,4> m, Nested<TV2> poly, bool is_invert = false);

extern Array<TV2> mul_contour(Matrix<T,4> m, Array<TV2> contour, bool is_invert = false);

extern Mesh mul(Matrix<T,4> m, Mesh soup, bool is_invert = false);

extern Mesh cone_mesh(T len, Array<TV2> poly);

extern Mesh cone_mesh(T len, Nested<TV2> contours);

extern Mesh taper_mesh(T len, T r0, T r1, Array<TV2> poly);

extern Mesh taper_mesh(T len, T r0, T r1, Nested<TV2> contours);

extern Mesh check_mesh(Mesh mesh);

extern Mesh extrude(T len, Nested<TV2> poly);

extern Mesh revolve(int n, Array<TV2> poly);

extern Mesh revolve(int n, Nested<TV2> contours);

extern Nested<TV2> fat_square_edge(int n, T rad, TV2 from, TV2 to);

extern Mesh fat_edge(int n, T rad, TV3 from, TV3 to);

extern Nested<TV2> fat_edge(int n, T rad, TV2 from, TV2 to);

extern Nested<TV2> fat_dot(int n, T rad, TV2 pt);

extern Nested<TV2> offset_poly(int n, T rad, Nested<TV2> poly);

extern Nested<TV2> offset_polyline(int n, T rad, Nested<TV2> polyline);

extern Mesh thicken(int n, T rad, Mesh mesh);

extern Mesh thicken(int n, T rad, Nested<TV3> polyline);

extern Nested<TV2> thicken(int n, T rad, Nested<TV2> line);

// extern Mesh offset_mesh(int n, T rad, Mesh mesh);

extern Nested<TV2> stroke_char (char letter);

extern Nested<TV2> stroke_text (std::string txt);
      
extern Mesh dither_mesh(Mesh mesh, double delta);

extern Mesh fab_mesh (Array<IV3> faces, Array<TV3> points);

extern void init_cad ( void );

inline T degrees_to_radians (T d) {
  return d * M_PI / 180.0;
}

#endif
