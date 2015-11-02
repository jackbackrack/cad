#include "cad.h"
#include "geom.h"

#include <cstdio>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

double rndd () {
  return (double)((double)rand() / (double)RAND_MAX);
}

double rndd (double mn, double mx) {
  return rndd() * (mx-mn) + mn;
}

double g_val(Geom* g) {
  ensure(g->k == float_kind, "NOT FLOAT");
  return ((FloatGeom*)g)->val;
}

std::string g_string(Geom* g) {
  ensure(g->k == string_kind, "NOT STRING");
  return ((StringGeom*)g)->val;
}

TV2 g_vec2(Geom* g) {
  ensure(g->k == vec2_kind, "NOT VEC2");
  return ((Vec2Geom*)g)->val;
}

TV g_vec3(Geom* g) {
  ensure(g->k == vec3_kind, "NOT VEC3");
  return ((Vec3Geom*)g)->val;
}

TV g_vec(Geom* g) { return g_vec3(g); }

Matrix<T,4> g_mat(Geom* g) {
  ensure(g->k == mat_kind, "NOT MAT");
  return ((MatGeom*)g)->val;
}

Array<TV2> g_line2(Geom* g) {
  ensure(g->k == line2_kind, "NOT LINE2");
  return ((Line2Geom*)g)->val;
}

Array<TV> g_line3(Geom* g) {
  ensure(g->k == line3_kind, "NOT LINE3");
  return ((Line3Geom*)g)->val;
}

Array<IV> g_faces(Geom* g) {
  ensure(g->k == line3_kind, "NOT FACES");
  auto points =  ((Line3Geom*)g)->val;
  Array<IV> faces;
  for (auto p : points)
    faces.append(vec((int)p.x, (int)p.y, (int)p.z));
  return faces;
}

bool is_polyline2(Geom* g) {
  return g->k == polyline2_kind || g->k == line2_kind;
}

Nested<TV2> g_polyline2(Geom* g) {
  ensure(is_polyline2(g), "NOT POLYLINE3");
  if (g->k == line2_kind)
    return line_to_polyline(g_line2(g));
  else
    return ((PolyLine2Geom*)g)->val;
}

bool is_polyline3(Geom* g) {
  return g->k == polyline3_kind || g->k == line3_kind;
}

Nested<TV> g_polyline3(Geom* g) {
  ensure(is_polyline3(g), "NOT POLYLINE3");
  if (g->k == line3_kind)
    return line_to_polyline(g_line3(g));
  else
    return ((PolyLine3Geom*)g)->val;
}

Array<TV2> g_contour(Geom* g) {
  ensure(g->k == contour_kind, "NOT CONTOUR");
  return ((ContourGeom*)g)->val;
}

bool is_poly(Geom* g) {
  return g->k == contour_kind || g->k == poly_kind;
}

Nested<TV2> g_poly(Geom* g) {
  ensure(is_poly(g), "NOT POLY"); 
  if (g->k == contour_kind)
    return contour_to_poly(g_contour(g));
  else // if (g->k == poly_kind)
    return ((PolyGeom*)g)->val;
}

Mesh g_mesh(Geom* g) {
  ensure(g->k == mesh_kind, "NOT MESH");
  return ((MeshGeom*)g)->val;
}

Geom* g_load(Geom* s) { 
  // TODO: LOAD SVG
  return new MeshGeom(read_soup(g_string(s)));
}

Geom* g_save(Geom* s, Geom* g) { 
  // TODO: LOAD SVG
  auto mesh = g_mesh(g);
  write_mesh(g_string(s), mesh.soup, mesh.points);
  return g;
}

Geom* g_print(Geom* g) { 
  if (g->k == mesh_kind)
    print_soup(g_mesh(g));
  else if (g->k == poly_kind)
    print_poly(g_poly(g));
  else if (g->k == contour_kind)
    print_contour(g_contour(g));
  else if (g->k == mat_kind)
    print_matrix(g_mat(g));
  else if (g->k == line2_kind)
    print_line2(g_line2(g));
  else if (g->k == line3_kind)
    print_line3(g_line3(g));
  else if (g->k == polyline2_kind)
    print_polyline2(g_polyline2(g));
  else if (g->k == polyline3_kind)
    print_polyline3(g_polyline3(g));
  return g;
}

Geom* g_pretty_print(Geom* g) { 
  if (g->k == mesh_kind)
    pretty_print_soup(g_mesh(g));
  else if (g->k == poly_kind)
    pretty_print_poly(g_poly(g));
  else if (g->k == contour_kind)
    pretty_print_contour(g_contour(g));
  else if (g->k == mat_kind)
    pretty_print_matrix(g_mat(g));
  else if (g->k == line2_kind)
    pretty_print_line2(g_line2(g));
  else if (g->k == line3_kind)
    pretty_print_line3(g_line3(g));
  else if (g->k == polyline2_kind)
    pretty_print_polyline2(g_polyline2(g));
  else if (g->k == polyline3_kind)
    pretty_print_polyline3(g_polyline3(g));
  return g;
}

Geom* g_float(double a) { return new FloatGeom(a); }
Geom* g_pi(void) { return new FloatGeom(M_PI); }
Geom* g_none2(void) { return new PolyGeom(none_poly()); }
Geom* g_all2(void) { return new PolyGeom(all_poly()); }
Geom* g_none(void) { return new MeshGeom(none_mesh()); }
Geom* g_all(void) { return new MeshGeom(all_mesh()); }
Geom* g_circle(Geom* a) { return new PolyGeom(circle_poly(g_val(a), 16)); }
Geom* g_square(Geom* a) { return new PolyGeom(square_poly(g_val(a))); }
Geom* g_square(Geom* lo, Geom* hi) { return new PolyGeom(square_poly(g_vec2(lo), g_vec2(hi))); }
Geom* g_letter(Geom* a) {
  char c = g_string(a)[0];
  auto ol = stroke_char(c);
  return new PolyLine2Geom(ol);
}
Geom* g_text(Geom* a) {
  auto txt = g_string(a);
  auto res = stroke_text(txt);
  return new PolyLine2Geom(res);
}
Geom* g_elt(Geom* g, Geom* i) {
  if (g->k == poly_kind)
    return new ContourGeom(poly_to_contour(g_poly(g), (int)g_val(i)));
  else {
    error("Bad arg for elt"); return NULL;
  }
}
Geom* g_xxx(Matrix<T,4> m, Geom* g) { 
  if (g->k == mesh_kind)
    return new MeshGeom(xxx(m, g_mesh(g)));
  else if (g->k == poly_kind)
    return new PolyGeom(xxx(m, g_poly(g)));
  else if (g->k == contour_kind)
    return new ContourGeom(xxx(m, g_contour(g)));
  else if (g->k == line2_kind)
    return new Line2Geom(xxx(m, g_line2(g)));
  else if (g->k == line3_kind)
    return new Line3Geom(xxx(m, g_line3(g)));
  else if (g->k == polyline2_kind)
    return new PolyLine2Geom(xxx(m, g_polyline2(g)));
  else if (g->k == polyline3_kind)
    return new PolyLine3Geom(xxx(m, g_polyline3(g)));
  else {
    error("Bad xxx kind"); return NULL;
  }
}
Geom* g_mag(Geom* v, Geom* g) { 
  return g_xxx(scale_matrix(g_vec3(v)), g);
}
Geom* g_mag1(Geom* a, Geom* g) { 
  return g_xxx(scale_matrix(vec(g_val(a), g_val(a), g_val(a))), g);
}
Geom* g_xmag(Geom* a, Geom* g) { 
  return g_xxx(scale_matrix(vec(g_val(a), 1.0, 1.0)), g);
}
Geom* g_ymag(Geom* a, Geom* g) { 
  return g_xxx(scale_matrix(vec(1.0, g_val(a), 1.0)), g);
}
Geom* g_zmag(Geom* a, Geom* g) { 
  return g_xxx(scale_matrix(vec(1.0, 1.0, g_val(a))), g);
}
Geom* g_mov(Geom* v, Geom* g) {
  return g_xxx(translation_matrix(g_vec3(v)), g);
}
Geom* g_xmov(Geom* a, Geom* g) {
  return g_xxx(translation_matrix(vec(g_val(a), 0.0, 0.0)), g);
}
Geom* g_ymov(Geom* a, Geom* g) {
  return g_xxx(translation_matrix(vec(0.0, g_val(a), 0.0)), g);
}
Geom* g_zmov(Geom* a, Geom* g) {
  return g_xxx(translation_matrix(vec(0.0, 0.0, g_val(a))), g);
}
Geom* g_rot(Geom* from, Geom* to, Geom* g) {
  return g_xxx(rotation_matrix(g_vec3(from), g_vec3(to)), g);
}
Geom* g_xrot(Geom* a, Geom* g) {
  T c=cos(g_val(a)),s=sin(g_val(a));
  return g_xxx(Matrix<T,4>(1,0,0,0,0,c,s,0,0,-s,c,0,0,0,0,1), g);
}
Geom* g_yrot(Geom* a, Geom* g) {
  T c=cos(g_val(a)),s=sin(g_val(a));
  return g_xxx(Matrix<T,4>(c,0,-s,0,0,1,0,0,s,0,c,0,0,0,0,1), g);
}
Geom* g_zrot(Geom* a, Geom* g) {
  T c=cos(g_val(a)),s=sin(g_val(a));
  return g_xxx(Matrix<T,4>(c,s,0,0,-s,c,0,0,0,0,1,0,0,0,0,1), g);
}
Geom* g_reflect_x(Geom* g) {
  return g_xxx(Matrix<T,4>(-1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1), g);
}
Geom* g_reflect_y(Geom* g) { 
  return g_xxx(Matrix<T,4>(1,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,1), g);
}
Geom* g_reflect_xy(Geom* g) {
  return g_xxx(Matrix<T,4>(-1,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,1), g);
}
Geom* g_reflect_z(Geom* g) {
  return g_xxx(Matrix<T,4>(1,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,1), g);
}
Geom* g_reflect_yz(Geom* g) {
  return g_xxx(Matrix<T,4>(1,0,0,0,-1,0,0,0,-1,0,0,0,0,0,0,1), g);
}
Geom* g_reflect_xz(Geom* g) {
  return g_xxx(Matrix<T,4>(-1,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,1), g);
}
Geom* g_add(Geom* a, Geom* b) {
  if (a->k == mesh_kind && b->k == mesh_kind)
    return new MeshGeom(add(g_mesh(a), g_mesh(b)));
  else if (is_poly(a) && is_poly(b))
    return new PolyGeom(add(g_poly(a), g_poly(b)));
  else {
    error("Bad args for add"); return NULL;
  }
}
Geom* g_mul(Geom* a, Geom* b) { 
  if (a->k == mat_kind)
    return g_xxx(g_mat(a), b);
  else if (a->k == mesh_kind && b->k == mesh_kind)
    return new MeshGeom(mul(g_mesh(a), g_mesh(b)));
  else if (is_poly(a) && is_poly(b))
    return new PolyGeom(mul(g_poly(a), g_poly(b)));
  else {
    error("Bad args for mul"); return NULL;
  }
}
Geom* g_sub(Geom* a, Geom* b) {
  if (a->k == mesh_kind && b->k == mesh_kind)
    return new MeshGeom(sub(g_mesh(a), g_mesh(b)));
  else if (is_poly(a) && is_poly(b))
    return new PolyGeom(sub(g_poly(a), g_poly(b)));
  else {
    error("Bad args for sub"); return NULL;
  }
}
Geom* g_not(Geom* a) {
  if (a->k == mesh_kind)
    return new MeshGeom(invert_mesh(g_mesh(a)));
  else if (is_poly(a))
    return new PolyGeom(invert_poly(g_poly(a)));
  else {
    error("Bad args for not"); return NULL;
  }
}
Geom* g_offset(Geom* a, Geom* g) {
  if (g->k == mesh_kind)
    return new MeshGeom(offset_mesh(1, g_val(a), g_mesh(g)));
  else if (is_poly(g))
    return new PolyGeom(offset_poly(16, g_val(a), g_poly(g)));
  else if (is_polyline2(g))
    return new PolyGeom(offset_polyline(16, g_val(a), g_polyline2(g)));
  else {
    error("Bad args for offset"); return NULL;
  }
}
Geom* g_simplify(Geom* g) { return new MeshGeom(simplify_mesh(g_mesh(g))); }
Geom* g_slice(Geom* a, Geom* g) { return new PolyGeom(slice(g_val(a), g_mesh(g))); }
Geom* g_extrude(Geom* a, Geom* p) { return new MeshGeom(extrude(g_val(a), g_poly(p))); }
Geom* g_thicken(Geom* a, Geom* l) {
  if (is_polyline2(l))
    return new PolyGeom(thicken(1, g_val(a), g_polyline2(l)));
  else //  (is_polyline3(l))
    return new MeshGeom(thicken(1, g_val(a), g_polyline3(l)));
}
Geom* g_sphere(Geom* a) { return new MeshGeom(sphere_mesh(2, vec(0.0, 0.0, 0.0), g_val(a))); }
Geom* g_cube(Geom* a) { auto r = g_val(a); return new MeshGeom(cube_mesh(vec(-r, -r, -r), vec(r, r, r))); }
Geom* g_cube(Geom* lo, Geom* hi) { return new MeshGeom(cube_mesh(g_vec3(lo), g_vec3(hi))); }
Geom* g_cone(Geom* a, Geom* p) { return new MeshGeom(cone_mesh(g_val(a), g_poly(p))); }
Geom* g_revolve(Geom* p) { return new MeshGeom(revolve(16, g_poly(p))); }
// Geom* g_hollow(Geom* a, Geom* m) { return new MeshGeom(hollow(g_val(a), g_mesh(m))); }
// Geom* g_shear_x_z(Geom* z0, Geom* z1, Geom* dx0, Geom* dx1, Geom* m) {
//   return new MeshGeom(shear_x_z(g_val(z0), g_val(z1), g_val(dx0), g_val(dx1), g_mesh(m))); }
Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p) {
  return new MeshGeom(taper_mesh(g_val(l), g_val(r0), g_val(r1), g_poly(p))); }
