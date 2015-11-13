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

Geom* g_args(void) { return new ArgsGeom(); }
Geom* g_args_add(Geom* g) { g_args_val(g).push_back(g); return g; }
int g_args_len(Geom* g) { return g_args_val(g).size(); }
std::vector<Geom*> g_args_val(Geom* g) {
  ensure(g->k == args_kind, "NOT ARGS");
  return ((ArgsGeom*)g)->val;
}

double g_num_val(Geom* g) {
  ensure(g->k == float_kind, "NOT FLOAT");
  return ((FloatGeom*)g)->val;
}

std::string g_string_val(Geom* g) {
  ensure(g->k == string_kind, "NOT STRING");
  return ((StringGeom*)g)->val;
}

Geom* g_string_fab(char* str) {
  std::string s(str);
  return g_string(s);
}

int g_string_len(Geom* g) {
  return g_string_val(g).size();
}

char* g_string_c_str(Geom* g) {
  return (char*)g_string_val(g).c_str();
}

TV2 g_vec2_val(Geom* g) {
  ensure(g->k == vec2_kind, "NOT VEC2");
  return ((Vec2Geom*)g)->val;
}

Geom* g_vec2(T x, T y) { return g_vec2(vec(x, y)); }
T g_vec2_elt(Geom* g, int idx)  { return g_vec2_val(g)[idx]; }
T g_vec2_x(Geom* g) { return g_vec2_val(g).x; }
T g_vec2_y(Geom* g) { return g_vec2_val(g).y; }

TV3 g_vec3_val(Geom* g) {
  ensure(g->k == vec3_kind, "NOT VEC3");
  return ((Vec3Geom*)g)->val;
}

Geom* g_vec3(T x, T y, T z) { return g_vec3(vec(x, y, z)); }
T g_vec3_elt(Geom* g, int idx)  { return g_vec3_val(g)[idx]; }
T g_vec3_x(Geom* g) { return g_vec3_val(g).x; }
T g_vec3_y(Geom* g) { return g_vec3_val(g).y; }
T g_vec3_z(Geom* g) { return g_vec3_val(g).z; }

Matrix<T,4> g_mat_val(Geom* g) {
  ensure(g->k == mat_kind, "NOT MAT");
  return ((MatGeom*)g)->val;
}

Geom* g_mat(T i00, T i01, T i02, T i03, T i10, T i11, T i12, T i13,
            T i20, T i21, T i22, T i23, T i30, T i31, T i32, T i33) {
  Matrix<T,4> m(i00, i01, i02, i03, i10, i11, i12, i13, i20, i21, i22, i23, i30, i31, i32, i33);
  return g_mat(m);
}

T g_mat_elt(Geom* g, int i, int j) {
  return g_mat_val(g).x[i][j];
}

Array<TV2> g_line2_val(Geom* g) {
  ensure(g->k == line2_kind, "NOT LINE2");
  return ((Line2Geom*)g)->val;
}

Geom* g_line2_fab(Geom* args) {
  Array<TV2> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_vec2_val(arg));
  return g_line2(v);
}
Geom* g_line2_elt(Geom* g, int idx) { return g_vec2(g_line2_val(g)[idx]); }
int g_line2_len(Geom* g) { return g_line2_val(g).size(); }

Array<TV3> g_line3_val(Geom* g) {
  ensure(g->k == line3_kind, "NOT LINE3");
  return ((Line3Geom*)g)->val;
}

Geom* g_line3_fab(Geom* args) {
  Array<TV3> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_vec3_val(arg));
  return g_line3(v);
}
Geom* g_line3_elt(Geom* g, int idx) { return g_vec3(g_line3_val(g)[idx]); }
int g_line3_len(Geom* g) { return g_line3_val(g).size(); }

Array<IV3> g_faces_val(Geom* g) {
  ensure(g->k == faces_kind, "NOT FACES");
  return ((FacesGeom*)g)->val;
}

Geom* g_faces_fab(Geom* args) {
  Array<IV3> faces;
  for (auto arg : g_args_val(args)) {
    auto v = g_vec3_val(arg);
    faces.append(vec((int)(v.x), (int)(v.y), (int)(v.z)));
  }
  return g_faces(faces);
}
Geom* g_faces_elt(Geom* g, int idx) {
  auto v = g_faces_val(g)[idx];
  return g_vec3(vec((T)(v.x), (T)(v.y), (T)(v.z)));
}
int g_faces_len(Geom* g) { return g_faces_val(g).size(); }

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

Geom* g_polyline2_fab(Geom* args) {
  Nested<TV2, false> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_line2_val(arg));
  v.freeze();
  return g_polyline2(v);
}
Geom* g_polyline2_elt(Geom* g, int idx) { return g_line2(g_polyline2_val(g)[idx]); }
int g_polyline2_len(Geom* g) { return g_polyline2_val(g).size(); }

bool is_polyline3(Geom* g) {
  return g->k == polyline3_kind || g->k == line3_kind;
}

Nested<TV3> g_polyline3_val(Geom* g) {
  ensure(is_polyline3(g), "NOT POLYLINE3");
  if (g->k == line3_kind)
    return line_to_polyline(g_line3_val(g));
  else
    return ((PolyLine3Geom*)g)->val;
}

Geom* g_polyline3_fab(Geom* args) {
  Nested<TV3, false> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_line3_val(arg));
  v.freeze();
  return g_polyline3(v);
}
Geom* g_polyline3_elt(Geom* g, int idx) { return g_line3(g_polyline3_val(g)[idx]); }
int g_polyline3_len(Geom* g) { return g_polyline3_val(g).size(); }

Array<TV2> g_contour_val(Geom* g) {
  ensure(g->k == contour_kind, "NOT CONTOUR");
  return ((ContourGeom*)g)->val;
}

Geom* g_contour_fab(Geom* args) {
  Array<TV2> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_vec2_val(arg));
  return g_contour(v);
}
Geom* g_contour_elt(Geom* g, int idx) { return g_vec2(g_contour_val(g)[idx]); }
int g_contour_len(Geom* g) { return g_contour_val(g).size(); }

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

Geom* g_poly_fab(Geom* args) {
  Nested<TV2, false> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_contour_val(arg));
  v.freeze();
  return g_poly(v);
}
Geom* g_poly_elt(Geom* g, int idx) { return g_contour(g_poly_val(g)[idx]); }
int g_poly_len(Geom* g) { return g_poly_val(g).size(); }

Mesh g_mesh_val(Geom* g) {
  ensure(g->k == mesh_kind, "NOT MESH");
  return ((MeshGeom*)g)->val;
}

Geom* g_mesh_fab(Geom* faces, Geom* vertices) {
  return g_mesh(fab_mesh(g_faces_val(faces), g_line3_val(vertices)));
}
Geom* g_mesh_faces(Geom* g) {
  auto faces = g_mesh_val(g).soup->elements;
  Array<IV3> a;
  for (auto e : faces)
    a.append(vec(e.x, e.y, e.z));
  return g_faces(a); 
}
Geom* g_mesh_points(Geom* g) {
  return g_line3(g_mesh_val(g).points);
}

Box<TV3> bbox(Array<TV3> points) {
  return bounding_box(points);
  
}
Box<TV2> bbox(Array<TV2> points) {
  return bounding_box(points);
}
Box<TV3> bbox(Nested<TV3> points) {
  Box<TV3> box = Box<TV3>::empty_box();
  for (auto line : points)
    box.enlarge(bounding_box(line));
  return box;
}
Box<TV2> bbox(Nested<TV2> points) {
  Box<TV2> box = Box<TV2>::empty_box();
  for (auto line : points)
    box.enlarge(bounding_box(line));
  return box;
}
  
Geom* g_bbox2(Box<TV2> box) {
  Array<TV2> bb; bb.append(box.min); bb.append(box.max);
  return g_line2(bb);
}
Geom* g_bbox3(Box<TV3> box) {
  Array<TV3> bb; bb.append(box.min); bb.append(box.max);
  return g_line3(bb);
}

Geom* g_bbox2_min(Geom* g) { return g_line2_elt(g, 0); }
Geom* g_bbox2_max(Geom* g) { return g_line2_elt(g, 1); }
Geom* g_bbox3_min(Geom* g) { return g_line3_elt(g, 0); }
Geom* g_bbox3_max(Geom* g) { return g_line3_elt(g, 1); }

Geom* g_bbox(Geom* g) {
  if (g->k == mesh_kind)
    return g_bbox3(bbox(g_mesh_val(g).points));
  else if (g->k == poly_kind)
    return g_bbox2(bbox(g_poly_val(g)));
  else if (g->k == contour_kind)
    return g_bbox2(bbox(g_contour_val(g)));
  else if (g->k == line2_kind)
    return g_bbox2(bbox(g_line2_val(g)));
  else if (g->k == line3_kind)
    return g_bbox3(bbox(g_line3_val(g)));
  else if (g->k == polyline2_kind)
    return g_bbox2(bbox(g_polyline2_val(g)));
  else if (g->k == polyline3_kind)
    return g_bbox3(bbox(g_polyline3_val(g)));
  else {
    error("BBOX UNDEFINED"); return NULL;
  }
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

std::string g_to_str_val(Geom* g) { 
  if (g->k == mesh_kind)
    return mesh_to_str(g_mesh_val(g));
  else if (g->k == poly_kind)
    return poly_to_str(g_poly_val(g));
  else if (g->k == contour_kind)
    return contour_to_str(g_contour_val(g));
  else if (g->k == mat_kind)
    return matrix_to_str(g_mat_val(g));
  else if (g->k == line2_kind)
    return line2_to_str(g_line2_val(g));
  else if (g->k == line3_kind)
    return line3_to_str(g_line3_val(g));
  else if (g->k == polyline2_kind)
    return polyline2_to_str(g_polyline2_val(g));
  else if (g->k == polyline3_kind)
    return polyline3_to_str(g_polyline3_val(g));
  return "-1";
}

Geom* g_to_str(Geom* g) {
  return g_string(g_to_str_val(g));
}

Geom* g_print(Geom* g) { 
  printf("%s\n", g_to_str_val(g).c_str());
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
Geom* g_vec3(TV3 v) { return new Vec3Geom(v); }
Geom* g_mat(Matrix<T,4> mat)  { return new MatGeom(mat); }
Geom* g_line2(Array<TV2> line) { return new Line2Geom(line); }
Geom* g_line2(RawArray<TV2> line) { return g_line2(line.copy()); }
Geom* g_line3(Array<TV3> line) { return new Line3Geom(line); }
Geom* g_line3(RawArray<TV3> line) { return g_line3(line.copy()); }
Geom* g_faces(Array<IV3> faces) { return new FacesGeom(faces); }
Geom* g_polyline2(Nested<TV2> polyline) { return new PolyLine2Geom(polyline); }
Geom* g_polyline3(Nested<TV3> polyline) { return new PolyLine3Geom(polyline); }
Geom* g_contour(Array<TV2> contour) { return new ContourGeom(contour); }
Geom* g_contour(RawArray<TV2> contour) { return g_contour(contour.copy()); }
Geom* g_poly(Nested<TV2> poly) { return new PolyGeom(poly); }
Geom* g_mesh(Mesh mesh) { return new MeshGeom(mesh); }
Geom* g_pi(void) { return new FloatGeom(M_PI); }
Geom* g_none2(void) { return new PolyGeom(none_poly()); }
Geom* g_all2(void) { return new PolyGeom(all_poly()); }
Geom* g_none(void) { return new MeshGeom(none_mesh()); }
Geom* g_all(void) { return new MeshGeom(all_mesh()); }
Geom* g_circle(Geom* a) { return new PolyGeom(circle_poly(g_num_val(a), 16)); }
Geom* g_square(Geom* a) { return new PolyGeom(square_poly(g_num_val(a))); }
Geom* g_square_lo_hi(Geom* lo, Geom* hi) { return new PolyGeom(square_poly(g_vec2_val(lo), g_vec2_val(hi))); }
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
Geom* g_cube_lo_hi(Geom* lo, Geom* hi) { return new MeshGeom(cube_mesh(g_vec3_val(lo), g_vec3_val(hi))); }
Geom* g_cone(Geom* a, Geom* p) { return new MeshGeom(cone_mesh(g_num_val(a), g_poly_val(p))); }
Geom* g_revolve(Geom* p) { return new MeshGeom(revolve(16, g_poly_val(p))); }
Geom* g_hull(Geom* m) { return new MeshGeom(quick_hull(g_mesh_val(m))); }
// Geom* g_hollow(Geom* a, Geom* m) { return new MeshGeom(hollow(g_num_val(a), g_mesh(m))); }
// Geom* g_shear_x_z(Geom* z0, Geom* z1, Geom* dx0, Geom* dx1, Geom* m) {
//   return new MeshGeom(shear_x_z(g_num_val(z0), g_num_val(z1), g_num_val(dx0), g_num_val(dx1), g_mesh(m))); }
Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p) {
  return new MeshGeom(taper_mesh(g_num_val(l), g_num_val(r0), g_num_val(r1), g_poly_val(p))); }
