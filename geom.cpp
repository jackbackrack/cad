#include "cad.h"
#include "hull.h"
#include "geom.h"

#include <cstdio>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

extern Mesh quick_hull(Mesh mesh);

double g_num_val(Geom* g) {
  ensure(g->k == float_kind, "NOT FLOAT");
  return ((FloatGeom*)g)->val;
}

std::string g_string_val(Geom* g) {
  ensure(g->k == string_kind, "NOT STRING");
  return ((StringGeom*)g)->val;
}

TV2 g_vec2_val(Geom* g) {
  ensure(g->k == vec2_kind, "NOT VEC2");
  return ((Vec2Geom*)g)->val;
}

TV g_vec3_val(Geom* g) {
  ensure(g->k == vec3_kind, "NOT VEC3");
  return ((Vec3Geom*)g)->val;
}

Matrix<T,4> g_mat_val(Geom* g) {
  ensure(g->k == mat_kind, "NOT MAT");
  return ((MatGeom*)g)->val;
}

Array<TV2> g_line2_val(Geom* g) {
  ensure(g->k == line2_kind, "NOT LINE2");
  return ((Line2Geom*)g)->val;
}

Array<TV> g_line3_val(Geom* g) {
  ensure(g->k == line3_kind, "NOT LINE3");
  return ((Line3Geom*)g)->val;
}

Array<IV> g_faces_val(Geom* g) {
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

Nested<TV2> g_polyline2_val(Geom* g) {
  ensure(is_polyline2(g), "NOT POLYLINE3");
  if (g->k == line2_kind)
    return line_to_polyline(g_line2_val(g));
  else
    return ((PolyLine2Geom*)g)->val;
}

bool is_polyline3(Geom* g) {
  return g->k == polyline3_kind || g->k == line3_kind;
}

Nested<TV> g_polyline3_val(Geom* g) {
  ensure(is_polyline3(g), "NOT POLYLINE3");
  if (g->k == line3_kind)
    return line_to_polyline(g_line3_val(g));
  else
    return ((PolyLine3Geom*)g)->val;
}

Array<TV2> g_contour_val(Geom* g) {
  ensure(g->k == contour_kind, "NOT CONTOUR");
  return ((ContourGeom*)g)->val;
}

bool is_poly(Geom* g) {
  return g->k == contour_kind || g->k == poly_kind;
}

Nested<TV2> g_poly_val(Geom* g) {
  ensure(is_poly(g), "NOT POLY"); 
  if (g->k == contour_kind)
    return contour_to_poly(g_contour_val(g));
  else // if (g->k == poly_kind)
    return ((PolyGeom*)g)->val;
}

Mesh g_mesh_val(Geom* g) {
  ensure(g->k == mesh_kind, "NOT MESH");
  return ((MeshGeom*)g)->val;
}

Geom* g_load(Geom* s) { 
  // TODO: LOAD SVG
  return new MeshGeom(read_soup(g_string_val(s)));
}

Geom* g_save(Geom* s, Geom* g) { 
  // TODO: LOAD SVG
  auto mesh = g_mesh_val(g);
  write_mesh(g_string_val(s), mesh.soup, mesh.points);
  return g;
}

Geom* g_print(Geom* g) { 
  if (g->k == mesh_kind)
    print_mesh(g_mesh_val(g));
  else if (g->k == poly_kind)
    print_poly(g_poly_val(g));
  else if (g->k == contour_kind)
    print_contour(g_contour_val(g));
  else if (g->k == mat_kind)
    print_matrix(g_mat_val(g));
  else if (g->k == line2_kind)
    print_line2(g_line2_val(g));
  else if (g->k == line3_kind)
    print_line3(g_line3_val(g));
  else if (g->k == polyline2_kind)
    print_polyline2(g_polyline2_val(g));
  else if (g->k == polyline3_kind)
    print_polyline3(g_polyline3_val(g));
  return g;
}

Geom* g_pretty_print(Geom* g) { 
  if (g->k == mesh_kind)
    pretty_print_mesh(g_mesh_val(g));
  else if (g->k == poly_kind)
    pretty_print_poly(g_poly_val(g));
  else if (g->k == contour_kind)
    pretty_print_contour(g_contour_val(g));
  else if (g->k == mat_kind)
    pretty_print_matrix(g_mat_val(g));
  else if (g->k == line2_kind)
    pretty_print_line2(g_line2_val(g));
  else if (g->k == line3_kind)
    pretty_print_line3(g_line3_val(g));
  else if (g->k == polyline2_kind)
    pretty_print_polyline2(g_polyline2_val(g));
  else if (g->k == polyline3_kind)
    pretty_print_polyline3(g_polyline3_val(g));
  return g;
}

Geom* g_num(double a) { return new FloatGeom(a); }
Geom* g_string(std::string s) { return new StringGeom(s); }
Geom* g_vec2(TV2 v) { return new Vec2Geom(v); }
Geom* g_vec3(TV v) { return new Vec3Geom(v); }
Geom* g_mat(Matrix<T,4> mat)  { return new MatGeom(mat); }
Geom* g_line2(Array<TV2> line) { return new Line2Geom(line); }
Geom* g_line3(Array<TV> line) { return new Line3Geom(line); }
Geom* g_faces(Array<TV> faces) { return new Line3Geom(faces); }
Geom* g_polyline2(Nested<TV2> polyline) { return new PolyLine2Geom(polyline); }
Geom* g_polyline3(Nested<TV> polyline) { return new PolyLine3Geom(polyline); }
Geom* g_contour(Array<TV2> contour) { return new ContourGeom(contour); }
Geom* g_poly(Nested<TV2> poly) { return new PolyGeom(poly); }
Geom* g_mesh(Mesh mesh) { return new MeshGeom(mesh); }
Geom* g_pi(void) { return new FloatGeom(M_PI); }
Geom* g_none2(void) { return new PolyGeom(none_poly()); }
Geom* g_all2(void) { return new PolyGeom(all_poly()); }
Geom* g_none(void) { return new MeshGeom(none_mesh()); }
Geom* g_all(void) { return new MeshGeom(all_mesh()); }
Geom* g_circle(Geom* a) { return new PolyGeom(circle_poly(g_num_val(a), 16)); }
Geom* g_square(Geom* a) { return new PolyGeom(square_poly(g_num_val(a))); }
Geom* g_square(Geom* lo, Geom* hi) { return new PolyGeom(square_poly(g_vec2_val(lo), g_vec2_val(hi))); }
Geom* g_letter(Geom* a) {
  char c = g_string_val(a)[0];
  auto ol = stroke_char(c);
  return new PolyLine2Geom(ol);
}
Geom* g_text(Geom* a) {
  auto txt = g_string_val(a);
  auto res = stroke_text(txt);
  return new PolyLine2Geom(res);
}
Geom* g_elt(Geom* g, Geom* i) {
  if (g->k == poly_kind)
    return new ContourGeom(poly_to_contour(g_poly_val(g), (int)g_num_val(i)));
  else {
    error("Bad arg for elt"); return NULL;
  }
}
Geom* g_add(Geom* a, Geom* b) {
  // TODO: FILL IN FOR NUMS AND VECS
  error("Bad args for add"); return NULL;
}
Geom* g_sub(Geom* a, Geom* b) {
  // TODO: FILL IN FOR NUMS AND VECS
  error("Bad args for sub"); return NULL;
}
Geom* do_g_mul(Matrix<T,4> m, Geom* g, bool is_invert = false) { 
  if (g->k == mesh_kind)
    return new MeshGeom(mul(m, g_mesh_val(g), is_invert));
  else if (g->k == poly_kind)
    return new PolyGeom(mul_poly(m, g_poly_val(g), is_invert));
  else if (g->k == contour_kind)
    return new ContourGeom(mul_contour(m, g_contour_val(g), is_invert));
  else if (g->k == line2_kind)
    return new Line2Geom(mul(m, g_line2_val(g)));
  else if (g->k == line3_kind)
    return new Line3Geom(mul(m, g_line3_val(g)));
  else if (g->k == polyline2_kind)
    return new PolyLine2Geom(mul(m, g_polyline2_val(g)));
  else if (g->k == polyline3_kind)
    return new PolyLine3Geom(mul(m, g_polyline3_val(g)));
  else {
    error("Bad mul kind"); return NULL;
  }
}
Geom* g_mul(Geom* a, Geom* b) { 
  // TODO: FILL IN FOR NUMS AND VECS
  if (a->k == mat_kind)
    return do_g_mul(g_mat_val(a), b);
  else {
    error("Bad args for mul"); return NULL;
  }
}
Geom* g_mag(Geom* v, Geom* g) { 
  return do_g_mul(scale_matrix(g_vec3_val(v)), g);
}
Geom* g_mag1(Geom* a, Geom* g) { 
  return do_g_mul(scale_matrix(vec(g_num_val(a), g_num_val(a), g_num_val(a))), g);
}
Geom* g_xmag(Geom* a, Geom* g) { 
  return do_g_mul(scale_matrix(vec(g_num_val(a), 1.0, 1.0)), g);
}
Geom* g_ymag(Geom* a, Geom* g) { 
  return do_g_mul(scale_matrix(vec(1.0, g_num_val(a), 1.0)), g);
}
Geom* g_zmag(Geom* a, Geom* g) { 
  return do_g_mul(scale_matrix(vec(1.0, 1.0, g_num_val(a))), g);
}
Geom* g_mov(Geom* v, Geom* g) {
  return do_g_mul(translation_matrix(g_vec3_val(v)), g);
}
Geom* g_xmov(Geom* a, Geom* g) {
  return do_g_mul(translation_matrix(vec(g_num_val(a), 0.0, 0.0)), g);
}
Geom* g_ymov(Geom* a, Geom* g) {
  return do_g_mul(translation_matrix(vec(0.0, g_num_val(a), 0.0)), g);
}
Geom* g_zmov(Geom* a, Geom* g) {
  return do_g_mul(translation_matrix(vec(0.0, 0.0, g_num_val(a))), g);
}
Geom* g_rot(Geom* from, Geom* to, Geom* g) {
  return do_g_mul(rotation_matrix(g_vec3_val(from), g_vec3_val(to)), g);
}
Geom* g_xrot(Geom* a, Geom* g) {
  T c=cos(g_num_val(a)),s=sin(g_num_val(a));
  return do_g_mul(Matrix<T,4>(1,0,0,0,0,c,s,0,0,-s,c,0,0,0,0,1), g);
}
Geom* g_yrot(Geom* a, Geom* g) {
  T c=cos(g_num_val(a)),s=sin(g_num_val(a));
  return do_g_mul(Matrix<T,4>(c,0,-s,0,0,1,0,0,s,0,c,0,0,0,0,1), g);
}
Geom* g_zrot(Geom* a, Geom* g) {
  T c=cos(g_num_val(a)),s=sin(g_num_val(a));
  return do_g_mul(Matrix<T,4>(c,s,0,0,-s,c,0,0,0,0,1,0,0,0,0,1), g);
}
Geom* g_reflect_x(Geom* g) {
  return do_g_mul(Matrix<T,4>(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1), g, true);
}
Geom* g_reflect_y(Geom* g) { 
  return do_g_mul(Matrix<T,4>(1,0,0,0, 0,-1,0,0, 0,0,1,0, 0,0,0,1), g, true);
}
Geom* g_reflect_xy(Geom* g) {
  return do_g_mul(Matrix<T,4>(-1,0,0,0, 0,-1,0,0, 0,0,1,0, 0,0,0,1), g, false); // TODO: NORMALS BROKEN FOR POLYS
}
Geom* g_reflect_z(Geom* g) {
  return do_g_mul(Matrix<T,4>(1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,1), g, true);
}
Geom* g_reflect_yz(Geom* g) {
  return do_g_mul(Matrix<T,4>(1,0,0,0, 0,-1,0,0, 0,0,-1,0,0,0,0,1), g, false);
}
Geom* g_reflect_xz(Geom* g) {
  return do_g_mul(Matrix<T,4>(-1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,1), g, false);
}
Geom* g_union(Geom* a, Geom* b) {
  if (a->k == mesh_kind && b->k == mesh_kind)
    return new MeshGeom(union_add(g_mesh_val(a), g_mesh_val(b)));
  else if (is_poly(a) && is_poly(b))
    return new PolyGeom(union_add(g_poly_val(a), g_poly_val(b)));
  else {
    error("Bad args for union"); return NULL;
  }
}
Geom* g_intersection(Geom* a, Geom* b) { 
  if (a->k == mesh_kind && b->k == mesh_kind)
    return new MeshGeom(intersection(g_mesh_val(a), g_mesh_val(b)));
  else if (is_poly(a) && is_poly(b))
    return new PolyGeom(intersection(g_poly_val(a), g_poly_val(b)));
  else {
    error("Bad args for intersection"); return NULL;
  }
}
Geom* g_difference(Geom* a, Geom* b) {
  if (a->k == mesh_kind && b->k == mesh_kind)
    return new MeshGeom(difference(g_mesh_val(a), g_mesh_val(b)));
  else if (is_poly(a) && is_poly(b))
    return new PolyGeom(difference(g_poly_val(a), g_poly_val(b)));
  else {
    error("Bad args for difference"); return NULL;
  }
}
Geom* g_not(Geom* a) {
  if (a->k == mesh_kind)
    return new MeshGeom(invert_mesh(g_mesh_val(a)));
  else if (is_poly(a))
    return new PolyGeom(invert_poly(g_poly_val(a)));
  else {
    error("Bad args for not"); return NULL;
  }
}
Geom* g_offset(Geom* a, Geom* g) {
  if (g->k == mesh_kind)
    return new MeshGeom(offset_mesh(1, g_num_val(a), g_mesh_val(g)));
  else if (is_poly(g))
    return new PolyGeom(offset_poly(16, g_num_val(a), g_poly_val(g)));
  else if (is_polyline2(g))
    return new PolyGeom(offset_polyline(16, g_num_val(a), g_polyline2_val(g)));
  else {
    error("Bad args for offset"); return NULL;
  }
}
Geom* g_simplify(Geom* g) { return new MeshGeom(real_simplify_mesh(g_mesh_val(g))); }
Geom* g_slice(Geom* a, Geom* g) { return new PolyGeom(slice(g_num_val(a), g_mesh_val(g))); }
Geom* g_extrude(Geom* a, Geom* p) { return new MeshGeom(extrude(g_num_val(a), g_poly_val(p))); }
Geom* g_thicken(Geom* a, Geom* l) {
  if (is_polyline2(l))
    return new PolyGeom(thicken(1, g_num_val(a), g_polyline2_val(l)));
  else //  (is_polyline3(l))
    return new MeshGeom(thicken(1, g_num_val(a), g_polyline3_val(l)));
}
Geom* g_sphere(Geom* a) { return new MeshGeom(sphere_mesh(1, vec(0.0, 0.0, 0.0), g_num_val(a))); }
Geom* g_cube(Geom* a) { auto r = g_num_val(a); return new MeshGeom(cube_mesh(vec(-r, -r, -r), vec(r, r, r))); }
Geom* g_cube(Geom* lo, Geom* hi) { return new MeshGeom(cube_mesh(g_vec3_val(lo), g_vec3_val(hi))); }
Geom* g_cone(Geom* a, Geom* p) { return new MeshGeom(cone_mesh(g_num_val(a), g_poly_val(p))); }
Geom* g_revolve(Geom* p) { return new MeshGeom(revolve(16, g_poly_val(p))); }
Geom* g_hull(Geom* m) { return new MeshGeom(quick_hull(g_mesh_val(m))); }
// Geom* g_hollow(Geom* a, Geom* m) { return new MeshGeom(hollow(g_num_val(a), g_mesh(m))); }
// Geom* g_shear_x_z(Geom* z0, Geom* z1, Geom* dx0, Geom* dx1, Geom* m) {
//   return new MeshGeom(shear_x_z(g_num_val(z0), g_num_val(z1), g_num_val(dx0), g_num_val(dx1), g_mesh(m))); }
Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p) {
  return new MeshGeom(taper_mesh(g_num_val(l), g_num_val(r0), g_num_val(r1), g_poly_val(p))); }
