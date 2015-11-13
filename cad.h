#ifndef __CAD__
#define __CAD__

#include <geode/mesh/io.h>
#include <geode/vector/Matrix3x3.h>
#include <geode/vector/Matrix4x4.h>
#include <geode/vector/Matrix4x4.h>
#include <geode/vector/Rotation.h>
#include <geode/geometry/platonic.h>
#include <geode/geometry/BoxVector.h>

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
// typedef Tuple<Ref<const TriangleSoup>, Array<TV3>> Mesh;
struct Mesh {
public:
  Ref<const TriangleSoup> soup;
  Array<TV3> points;
  Mesh(Ref<const TriangleSoup> soup, Array<TV3> points) : soup(soup), points(points) { }
};

extern Mesh fab_mesh (Array<IV3> faces, Array<TV3> points);

extern Matrix<T,4> to_matrix44(Matrix<T,3> M);
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
extern Mesh real_simplify_mesh(Mesh mesh);

extern Mesh invert_mesh(Mesh mesh);

extern Array<TV2> invert_contour(Array<TV2> contour);
extern Nested<TV2> invert_poly(Nested<TV2> poly);

extern Array<TV2> simplify_contour(RawArray<TV2> contour);
extern Nested<TV2> simplify_poly(Nested<TV2> poly);

extern Nested<TV2> union_add(Nested<TV2> c0, Nested<TV2> c1);
extern Nested<TV2> intersection(Nested<TV2> c0, Nested<TV2> c1);
extern Nested<TV2> difference(Nested<TV2> c0, Nested<TV2> c1);

extern Nested<TV2> offset(T a, Nested<TV2> c);

extern Mesh union_add(Mesh mesh0, Mesh mesh1, bool is_simplify = true);

extern void pretty_print_mesh(Mesh soup);
extern void pretty_print_line2(Array<TV2> line);
extern void pretty_print_line3(Array<TV3> line);
extern void pretty_print_contour(Array<TV2> contour);
extern void pretty_print_poly(Nested<TV2> poly);
extern void pretty_print_matrix(Matrix<T,4> M);
extern void pretty_print_polyline3(Nested<TV3> polyline);
extern void pretty_print_polyline2(Nested<TV2> polyline);

extern std::string mesh_to_str(Mesh soup);
extern std::string contour_to_str(Array<TV2> contour);
extern std::string line3_to_str(Array<TV3> contour);
extern std::string line2_to_str(Array<TV2> contour);
extern std::string poly_to_str(Nested<TV2> poly);
extern std::string polyline3_to_str(Nested<TV3> polyline);
extern std::string polyline2_to_str(Nested<TV2> polyline);
extern std::string matrix_to_str(Matrix<T,4> M);

extern Mesh intersection(Mesh mesh0, Mesh mesh1, bool is_simplify = true);

extern Mesh mesh_from(int start, Mesh soup);

extern Nested<TV2> slice(T z, Mesh soup);

extern Mesh difference(Mesh mesh0, Mesh mesh1, bool is_simplify = true);

template<class ET> Nested<ET> contour_to_poly(Array<ET> contour);

template<class ET> Array<ET> poly_to_contour(Nested<ET> poly, int i = 0);

extern Mesh all_mesh(void);

extern Mesh none_mesh(void);

extern Nested<TV2> square_poly(TV2 min, TV2 max);

extern Nested<TV2> square_poly(T rad);

extern Nested<TV2> all_poly(void);

extern Nested<TV2> none_poly(void);

extern Nested<TV2> circle_poly(T rad, int n);

extern Nested<TV2> star_poly(T rad_min, T rad_max, int n);


extern Mesh triangulate (Nested<TV3> poly);

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

extern Mesh offset_mesh(int n, T rad, Mesh mesh);

extern Nested<TV2> stroke_char (char letter);

extern Nested<TV2> stroke_text (std::string txt);
      
extern void init_cad ( void );

#endif
