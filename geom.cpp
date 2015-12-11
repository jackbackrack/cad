#include "cad.h"
#include "hull.h"
#include "iso-surface.h"
#include "geom.h"

#include <cstdio>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

int g_kind(Geom* g) { return g->k; }

extern Mesh quick_hull(Mesh mesh);

Geom* g_args_fab(void) { return new ArgsGeom(); }
Geom* g_args_add(Geom* g, Geom* e) {
  ((ArgsGeom*)g)->val.push_back(e);
  return g;
}
int g_args_len(Geom* g) { return g_args_val(g).size(); }
std::vector<Geom*> g_args_val(Geom* g) {
  ensure(g->k == args_kind, "NOT ARGS");
  return ((ArgsGeom*)g)->val;
}

inline float int_as_float(int xi) {
  float x =  *(float*)&xi;
  return x;
}

inline int float_as_int(float xf) {
  int x =  *(int*)&xf;
  return x;
}

Geom* g_num_fab(int a) {
  return g_num((T)int_as_float(a));
}
double g_num_val(Geom* g) {
  ensure(g->k == num_kind, "NOT NUM");
  return ((NumGeom*)g)->val;
}

int g_num_value(Geom* g) {
  return float_as_int(g_num_val(g));
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

Geom* g_v2d_fab(int x, int y) {
  auto res = g_v2d(vec((T)int_as_float(x), (T)int_as_float(y)));
  return res;
}
TV2 g_v2d_val(Geom* g) {
  ensure(g->k == v2d_kind, "NOT V2D");
  return ((V2dGeom*)g)->val;
}
int g_v2d_elt(Geom* g, int idx)  { return float_as_int(g_v2d_val(g)[idx]); }
int g_v2d_x(Geom* g) { return float_as_int(g_v2d_val(g).x); }
int g_v2d_y(Geom* g) { return float_as_int(g_v2d_val(g).y); }

Geom* g_v3d_fab(int x, int y, int z) {
  return g_v3d(vec((T)int_as_float(x), (T)int_as_float(y), (T)int_as_float(z)));
}
TV3 g_v3d_val(Geom* g) {
  ensure(g->k == v3d_kind, "NOT V3D");
  return ((V3dGeom*)g)->val;
}
int g_v3d_elt(Geom* g, int idx)  { return float_as_int(g_v3d_val(g)[idx]); }
int g_v3d_x(Geom* g) { return float_as_int(g_v3d_val(g).x); }
int g_v3d_y(Geom* g) { return float_as_int(g_v3d_val(g).y); }
int g_v3d_z(Geom* g) { return float_as_int(g_v3d_val(g).z); }

Geom* g_v3i_fab(int x, int y, int z) { return g_v3i(vec(x, y, z)); }
IV3 g_v3i_val(Geom* g) {
  ensure(g->k == v3i_kind || g->k == v3d_kind, "NOT V3I");
  if (g->k == v3d_kind) {
    auto v = g_v3d_val(g);
    return vec((int)(v.x), (int)(v.y), (int)(v.z));
  } else 
    return ((V3iGeom*)g)->val;
}
int g_v3i_elt(Geom* g, int idx)  { return g_v3i_val(g)[idx]; }
int g_v3i_x(Geom* g) { return g_v3i_val(g).x; }
int g_v3i_y(Geom* g) { return g_v3i_val(g).y; }
int g_v3i_z(Geom* g) { return g_v3i_val(g).z; }

Geom* g_mat_fab(int i00, int i01, int i02, int i03, int i10, int i11, int i12, int i13,
                int i20, int i21, int i22, int i23, int i30, int i31, int i32, int i33) {
  Matrix<T,4> m(int_as_float(i00), int_as_float(i01), int_as_float(i02), int_as_float(i03),
                int_as_float(i10), int_as_float(i11), int_as_float(i12), int_as_float(i13),
                int_as_float(i20), int_as_float(i21), int_as_float(i22), int_as_float(i23),
                int_as_float(i30), int_as_float(i31), int_as_float(i32), int_as_float(i33));
  return g_mat(m);
}
Matrix<T,4> g_mat_val(Geom* g) {
  ensure(g->k == mat_kind, "NOT MAT");
  return ((MatGeom*)g)->val;
}

int g_mat_elt(Geom* g, int i, int j) {
  return float_as_int(g_mat_val(g).x[i][j]);
}

Array<TV2> g_array_v2d_val(Geom* g) {
  ensure(g->k == array_v2d_kind, "NOT ARRAY_V2D");
  return ((ArrayV2dGeom*)g)->val;
}

Geom* g_array_v2d_fab(Geom* args) {
  Array<TV2> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_v2d_val(arg));
  return g_array_v2d(v);
}
Geom* g_array_v2d_elt(Geom* g, int idx) { return g_v2d(g_array_v2d_val(g)[idx]); }
int g_array_v2d_len(Geom* g) { return g_array_v2d_val(g).size(); }

Array<TV3> g_array_v3d_val(Geom* g) {
  ensure(g->k == array_v3d_kind, "NOT ARRAY_V3D");
  return ((ArrayV3dGeom*)g)->val;
}

Geom* g_array_v3d_fab(Geom* args) {
  Array<TV3> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_v3d_val(arg));
  return g_array_v3d(v);
}
Geom* g_array_v3d_elt(Geom* g, int idx) { return g_v3d(g_array_v3d_val(g)[idx]); }
int g_array_v3d_len(Geom* g) { return g_array_v3d_val(g).size(); }

Array<IV3> g_array_v3i_val(Geom* g) {
  ensure(g->k == array_v3i_kind || g->k == array_v3d_kind, "NOT ARRAY_V3I");
  if (g->k == array_v3d_kind) {
    Array<IV3> array_v3i;
    for (auto v : g_array_v3d_val(g)) 
      array_v3i.append(vec((int)(v.x), (int)(v.y), (int)(v.z)));
    return array_v3i;
  }
  return ((ArrayV3iGeom*)g)->val;
}

Geom* g_array_v3i_fab(Geom* args) {
  Array<IV3> array_v3i;
  for (auto arg : g_args_val(args)) 
    array_v3i.append(g_v3i_val(arg));
  return g_array_v3i(array_v3i);
}
Geom* g_array_v3i_elt(Geom* g, int idx) {
  return g_v3i(g_array_v3i_val(g)[idx]);
}
int g_array_v3i_len(Geom* g) { return g_array_v3i_val(g).size(); }


bool is_nested_v2d(Geom* g) {
  return g->k == nested_v2d_kind || g->k == array_v2d_kind;
}
Nested<TV2> g_nested_v2d_val(Geom* g) {
  ensure(is_nested_v2d(g), "NOT NESTED_V3D");
  if (g->k == array_v2d_kind)
    return array_to_nested(g_array_v2d_val(g));
  else
    return ((NestedV2dGeom*)g)->val;
}

Geom* g_nested_v2d_fab(Geom* args) {
  Nested<TV2, false> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_array_v2d_val(arg));
  v.freeze();
  return g_nested_v2d(v);
}
Geom* g_nested_v2d_elt(Geom* g, int idx) { return g_array_v2d(g_nested_v2d_val(g)[idx]); }
int g_nested_v2d_len(Geom* g) { return g_nested_v2d_val(g).size(); }

bool is_poly(Geom* g) {
  return g->k == poly_kind || g->k == array_v2d_kind;
}

Nested<TV2> g_poly_val(Geom* g) {
  ensure(is_poly(g), "NOT POLY");
  if (g->k == array_v2d_kind)
    return array_to_nested(g_array_v2d_val(g));
  else
    return ((PolyGeom*)g)->val;
}

Geom* g_poly_fab(Geom* args) {
  Nested<TV2, false> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_array_v2d_val(arg));
  v.freeze();
  return g_poly(v);
}
Geom* g_poly_elt(Geom* g, int idx) { return g_array_v2d(g_poly_val(g)[idx]); }
int g_poly_len(Geom* g) { return g_poly_val(g).size(); }

bool is_nested_v3d(Geom* g) {
  return g->k == nested_v3d_kind || g->k == array_v3d_kind;
}

Nested<TV3> g_nested_v3d_val(Geom* g) {
  ensure(is_nested_v3d(g), "NOT NESTED_V3D");
  if (g->k == array_v3d_kind)
    return line_to_polyline(g_array_v3d_val(g));
  else
    return ((NestedV3dGeom*)g)->val;
}

Geom* g_nested_v3d_fab(Geom* args) {
  Nested<TV3, false> v;
  for (auto arg : g_args_val(args)) 
    v.append(g_array_v3d_val(arg));
  v.freeze();
  return g_nested_v3d(v);
}
Geom* g_nested_v3d_elt(Geom* g, int idx) { return g_array_v3d(g_nested_v3d_val(g)[idx]); }
int g_nested_v3d_len(Geom* g) { return g_nested_v3d_val(g).size(); }


Mesh g_mesh_val(Geom* g) {
  ensure(g->k == mesh_kind, "NOT MESH");
  return ((MeshGeom*)g)->val;
}

Geom* g_mesh_fab(Geom* vertices, Geom* faces) {
  return g_mesh(fab_mesh(g_array_v3i_val(faces), g_array_v3d_val(vertices)));
}
Geom* g_mesh_faces(Geom* g) {
  auto faces = g_mesh_val(g).soup->elements;
  Array<IV3> a;
  for (auto e : faces)
    a.append(vec(e.x, e.y, e.z));
  return g_array_v3i(a); 
}
Geom* g_mesh_points(Geom* g) {
  return g_array_v3d(g_mesh_val(g).points);
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
  return g_array_v2d(bb);
}
Geom* g_bbox3(Box<TV3> box) {
  Array<TV3> bb; bb.append(box.min); bb.append(box.max);
  return g_array_v3d(bb);
}

Geom* g_bbox2_min(Geom* g) { return g_array_v2d_elt(g, 0); }
Geom* g_bbox2_max(Geom* g) { return g_array_v2d_elt(g, 1); }
Geom* g_bbox3_min(Geom* g) { return g_array_v3d_elt(g, 0); }
Geom* g_bbox3_max(Geom* g) { return g_array_v3d_elt(g, 1); }

Geom* g_bbox(Geom* g) {
  if (g->k == mesh_kind)
    return g_bbox3(bbox(g_mesh_val(g).points));
  else if (g->k == poly_kind)
    return g_bbox2(bbox(g_poly_val(g)));
  else if (g->k == array_v2d_kind)
    return g_bbox2(bbox(g_array_v2d_val(g)));
  else if (g->k == array_v3d_kind)
    return g_bbox3(bbox(g_array_v3d_val(g)));
  else if (g->k == nested_v2d_kind)
    return g_bbox2(bbox(g_nested_v2d_val(g)));
  else if (g->k == nested_v3d_kind)
    return g_bbox3(bbox(g_nested_v3d_val(g)));
  else {
    error("BBOX UNDEFINED"); return NULL;
  }
}
Geom* g_dims(Geom* g) {
  auto bb = g_bbox(g);
  return g_sub(g_elt(bb, g_num(1.0)), g_elt(bb, g_num(0.0)));
}
Geom* g_center(Geom* g) {
  auto bb = g_bbox(g);
  return g_mul(g_num(0.5), g_add(g_elt(bb, g_num(0.0)), g_elt(bb, g_num(1.0))));
}
Geom* g_centering(Geom* g) {
  auto c = g_center(g);
  return g_mov(g_mul(g_num(-1), c), g);
}
Geom* g_load(Geom* s) { 
  // TODO: LOAD SVG
  return new MeshGeom(read_soup(g_string_val(s)));
}

Geom* g_save(Geom* s, Geom* g) { 
  // TODO: LOAD SVG
  if (g->k == mesh_kind) {
    auto mesh = g_mesh_val(g);
    write_mesh(g_string_val(s), mesh.soup, mesh.points);
  } else if (g->k == poly_kind) {
    save_polygon(g_string_val(s), g_poly_val(g));
  } else if (g->k == nested_v2d_kind) {
    save_polyline(g_string_val(s), g_nested_v2d_val(g));
  }
  return g;
}

std::string g_to_str_val(Geom* g) { 
  if (g->k == num_kind)
    return num_to_str(g_num_val(g));
  else if (g->k == v2d_kind)
    return v2d_to_str(g_v2d_val(g));
  else if (g->k == v3d_kind)
    return v3d_to_str(g_v3d_val(g));
  else if (g->k == v3i_kind)
    return v3i_to_str(g_v3i_val(g));
  else if (g->k == mat_kind)
    return matrix_to_str(g_mat_val(g));
  else if (g->k == array_v2d_kind)
    return array_v2d_to_str(g_array_v2d_val(g));
  else if (g->k == array_v3d_kind)
    return array_v3d_to_str(g_array_v3d_val(g));
  else if (g->k == nested_v2d_kind)
    return nested_v2d_to_str(g_nested_v2d_val(g));
  else if (g->k == nested_v3d_kind)
    return nested_v3d_to_str(g_nested_v3d_val(g));
  else if (g->k == poly_kind)
    return poly_to_str(g_poly_val(g));
  else if (g->k == mesh_kind)
    return mesh_to_str(g_mesh_val(g));
  else
    return "-1";
}

Geom* g_to_str(Geom* g) {
  return g_string(g_to_str_val(g));
}

Geom* g_print(Geom* g) { 
  printf("%s\n", g_to_str_val(g).c_str());
  return g;
}
Geom* g_check(Geom* g) { 
  if (g->k == mesh_kind) 
    return g_mesh(check_mesh(g_mesh_val(g)));
  else {
    error("BAD CHECK"); return NULL;
  }
}

Geom* g_pretty_print(Geom* g) { 
  if (g->k == num_kind)
    pretty_print_num(g_num_val(g));
  else if (g->k == v2d_kind)
    pretty_print_v2d(g_v2d_val(g));
  else if (g->k == v3d_kind)
    pretty_print_v3d(g_v3d_val(g));
  else if (g->k == v3i_kind)
    pretty_print_v3i(g_v3i_val(g));
  else if (g->k == poly_kind)
    pretty_print_poly(g_poly_val(g));
  else if (g->k == mat_kind)
    pretty_print_matrix(g_mat_val(g));
  else if (g->k == array_v2d_kind)
    pretty_print_array_v2d(g_array_v2d_val(g));
  else if (g->k == array_v3d_kind)
    pretty_print_array_v3d(g_array_v3d_val(g));
  else if (g->k == nested_v2d_kind)
    pretty_print_nested_v2d(g_nested_v2d_val(g));
  else if (g->k == nested_v3d_kind)
    pretty_print_nested_v3d(g_nested_v3d_val(g));
  else if (g->k == mesh_kind)
    pretty_print_mesh(g_mesh_val(g));
  return g;
}

Geom* g_num(T a) { return new NumGeom(a); }
Geom* g_string(std::string s) { return new StringGeom(s); }
Geom* g_v2d(TV2 v) { return new V2dGeom(v); }
Geom* g_v3d(TV3 v) { return new V3dGeom(v); }
Geom* g_v3i(IV3 v) { return new V3iGeom(v); }
Geom* g_mat(Matrix<T,4> mat)  { return new MatGeom(mat); }
Geom* g_array_v2d(Array<TV2> line) { return new ArrayV2dGeom(line); }
Geom* g_array_v2d(RawArray<TV2> line) { return g_array_v2d(line.copy()); }
Geom* g_array_v3d(Array<TV3> line) { return new ArrayV3dGeom(line); }
Geom* g_array_v3d(RawArray<TV3> line) { return g_array_v3d(line.copy()); }
Geom* g_array_v3i(Array<IV3> line) { return new ArrayV3iGeom(line); }
Geom* g_array_v3i(RawArray<IV3> line) { return g_array_v3i(line.copy()); }
Geom* g_nested_v2d(Nested<TV2> polyline) { return new NestedV2dGeom(polyline); }
Geom* g_nested_v3d(Nested<TV3> polyline) { return new NestedV3dGeom(polyline); }
Geom* g_poly(Nested<TV2> poly) { return new PolyGeom(poly); }
Geom* g_mesh(Mesh mesh) { return new MeshGeom(mesh); }

bool all_args_kind (std::vector<Geom*> args, int kind) {
  bool res = true;
  for (auto arg : args) 
    res = res && arg->k == kind;
  return res;
}

Geom* g_array_v3d (std::vector<Geom*> args) {
  Array< TV3 > points;
  for (auto arg : args)
    points.append(g_v3d_val(arg));
  return new ArrayV3dGeom(points);
}

Geom* g_array_v3i (std::vector<Geom*> args) {
  Array< IV3 > points;
  for (auto arg : args)
    points.append(g_v3i_val(arg));
  return new ArrayV3iGeom(points);
}

Geom* g_array_v2d (std::vector<Geom*> args) {
  Array< TV2 > points;
  for (auto arg : args) {
    points.append(g_v2d_val(arg));
  }
  return new ArrayV2dGeom(points);
}

Geom* g_nested_v3d (std::vector<Geom*> args) {
  Nested< TV3,false > lines;
  for (auto arg : args)
    lines.append(g_array_v3d_val(arg));
  lines.freeze();
  return new NestedV3dGeom(lines);
}

Geom* g_nested_v2d (std::vector<Geom*> args) {
  Nested< TV2,false > lines;
  for (auto arg : args)
    lines.append(g_array_v2d_val(arg));
  lines.freeze();
  return new NestedV2dGeom(lines);
}

Geom* g_poly (std::vector<Geom*> args) {
  Nested< TV2,false > lines;
  if (all_args_kind(args, array_v2d_kind)) {
    for (auto arg : args)
      lines.append(g_array_v2d_val(arg));
  } else if (all_args_kind(args, v2d_kind)) {
    Array< TV2 > points;
    for (auto arg : args)
      points.append(g_v2d_val(arg));
    lines.append(points);
  } else {
    error("Bad Poly Args"); return NULL;
  }
  lines.freeze();
  return new PolyGeom(lines);
}

Geom* g_pi(void) { return new NumGeom(M_PI); }
Geom* g_none2(void) { return new PolyGeom(none_poly()); }
Geom* g_all2(void) { return new PolyGeom(all_poly()); }
Geom* g_none(void) { return new MeshGeom(none_mesh()); }
Geom* g_all(void) { return new MeshGeom(all_mesh()); }
Geom* g_circle(Geom* a) { return new PolyGeom(circle_poly(g_num_val(a), 16)); }
Geom* g_square(Geom* a) { return new PolyGeom(square_poly(g_num_val(a))); }
Geom* g_square_lo_hi(Geom* lo, Geom* hi) { return new PolyGeom(square_poly(g_v2d_val(lo), g_v2d_val(hi))); }
Geom* g_letter(Geom* a) {
  char c = g_string_val(a)[0];
  auto ol = stroke_char(c);
  return new NestedV2dGeom(ol);
}
Geom* g_text(Geom* a) {
  auto txt = g_string_val(a);
  auto res = stroke_text(txt);
  return new NestedV2dGeom(res);
}
Geom* g_elt(Geom* g, Geom* i) {
  if (g->k == poly_kind)
    return new ArrayV2dGeom(nested_elt(g_poly_val(g), (int)g_num_val(i)));
  else if (g->k == array_v2d_kind)
    return g_v2d(g_array_v2d_val(g)[(int)g_num_val(i)]);
  else if (g->k == array_v3d_kind)
    return g_v3d(g_array_v3d_val(g)[(int)g_num_val(i)]);
  else if (g->k == v3d_kind)
    return g_num(g_v3d_val(g)[(int)g_num_val(i)]);
  else if (g->k == v2d_kind)
    return g_num(g_v2d_val(g)[(int)g_num_val(i)]);
  else {
    error("Bad arg for elt"); return NULL;
  }
}
Geom* g_add(Geom* a, Geom* b) {
  if (a->k == num_kind && b->k == num_kind)
    return g_num(g_num_val(a) + g_num_val(b));
  else if (a->k == v3d_kind && b->k == v3d_kind)
    return g_v3d(g_v3d_val(a) + g_v3d_val(b));
  else if (a->k == v2d_kind && b->k == v2d_kind)
    return g_v2d(g_v2d_val(a) + g_v2d_val(b));
  else {
    error("Bad args for add"); return NULL;
  }
}
Geom* g_sub(Geom* a, Geom* b) {
  if (a->k == num_kind && b->k == num_kind)
    return g_num(g_num_val(a) - g_num_val(b));
  else if (a->k == v3d_kind && b->k == v3d_kind)
    return g_v3d(g_v3d_val(a) - g_v3d_val(b));
  else if (a->k == v2d_kind && b->k == v2d_kind)
    return g_v2d(g_v2d_val(a) - g_v2d_val(b));
  else {
    error("Bad args for sub"); return NULL;
  }
}
Geom* g_dither(Geom* g) { 
  if (g->k == mesh_kind)
    return new MeshGeom(dither_mesh(g_mesh_val(g), 1e-16));
  else {
    error("Bad args for dither"); return NULL;
  }
}
Geom* do_g_mul(Matrix<T,4> m, Geom* g, bool is_invert = false) { 
  if (g->k == mesh_kind)
    return new MeshGeom(mul(m, g_mesh_val(g), is_invert));
  else if (g->k == poly_kind)
    return g_poly(mul_poly(m, g_poly_val(g), is_invert));
  else if (g->k == array_v2d_kind)
    return g_array_v2d(mul(m, g_array_v2d_val(g)));
  else if (g->k == array_v3d_kind)
    return g_array_v3d(mul(m, g_array_v3d_val(g)));
  else if (g->k == nested_v2d_kind)
    return g_nested_v2d(mul(m, g_nested_v2d_val(g)));
  else if (g->k == nested_v3d_kind)
    return g_nested_v3d(mul(m, g_nested_v3d_val(g)));
  else {
    error("Bad mul kind"); return NULL;
  }
}
Geom* g_mul(Geom* a, Geom* b) { 
  if (a->k == mat_kind)
    return do_g_mul(g_mat_val(a), b);
  else if (a->k == num_kind && b->k == num_kind)
    return g_num(g_num_val(a) * g_num_val(b));
  else if (a->k == num_kind && b->k == v3d_kind)
    return g_v3d(g_num_val(a) * g_v3d_val(b));
  else if (a->k == num_kind && b->k == v2d_kind)
    return g_v2d(g_num_val(a) * g_v2d_val(b));
  else {
    error("Bad args for mul"); return NULL;
  }
}

Geom* g_normalize(Geom* a) {
  if (a->k == v3d_kind)
    return g_v3d(g_v3d_val(a).normalized());
  else if (a->k == v2d_kind)
    return g_v2d(g_v2d_val(a).normalized());
  else {
    error("Bad args for normalize"); return NULL;
  }
  
}

Geom* g_magnitude(Geom* a) {
  if (a->k == v3d_kind)
    return g_num(g_v3d_val(a).magnitude());
  else if (a->k == v2d_kind)
    return g_num(g_v2d_val(a).magnitude());
  else {
    error("Bad args for magnitude"); return NULL;
  }
  
}
Geom* g_cross(Geom* a, Geom* b) {
  if (a->k == v3d_kind && b->k == v3d_kind) {
    return g_v3d(cross(g_v3d_val(a), g_v3d_val(b)));
  } else {
    error("Bad args for cross"); return NULL;
  }
}
Geom* g_dot(Geom* a, Geom* b) {
  if (a->k == v3d_kind && b->k == v3d_kind)
    return g_num(dot(g_v3d_val(a), g_v3d_val(b)));
  else if (a->k == num_kind && b->k == v2d_kind)
    return g_num(dot(g_v2d_val(a), g_v2d_val(b)));
  else {
    error("Bad args for dot"); return NULL;
  }
}

Geom* g_div(Geom* a, Geom* b) { 
  if (a->k == num_kind && b->k == num_kind)
    return g_num(g_num_val(a) * g_num_val(b));
  else if (a->k == v3d_kind && b->k == num_kind)
    return g_v3d(g_v3d_val(a) / g_num_val(b));
  else if (a->k == v2d_kind && b->k == num_kind)
    return g_v2d(g_v2d_val(a) * g_num_val(b));
  else {
    error("Bad args for div"); return NULL;
  }
}
Geom* g_mag(Geom* v, Geom* g) { 
  return do_g_mul(scale_matrix(g_v3d_val(v)), g);
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
  return do_g_mul(translation_matrix(g_v3d_val(v)), g);
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
inline T degrees_to_radians (T d) {
  return d * M_PI / 180.0;
}
Geom* g_rot(Geom* v, Geom* g) {
  printf("ROT: ");
  g_pretty_print(v);
  TV3 rv = g_v3d_val(v);
  TV3 av = vec(degrees_to_radians(rv.x), degrees_to_radians(rv.y), degrees_to_radians(rv.z));
  return do_g_mul(rotation_matrix(av), g);
}
Geom* g_rot_from_to(Geom* from, Geom* to, Geom* g) {
  return do_g_mul(rotation_matrix(g_v3d_val(from), g_v3d_val(to)), g);
}
Geom* g_xrot(Geom* d, Geom* g) {
  T a = degrees_to_radians(g_num_val(d));
  T c=cos(a),s=sin(a);
  return do_g_mul(Matrix<T,4>(1,0,0,0,0,c,s,0,0,-s,c,0,0,0,0,1), g);
}
Geom* g_yrot(Geom* d, Geom* g) {
  T a = degrees_to_radians(g_num_val(d));
  T c=cos(a),s=sin(a);
  return do_g_mul(Matrix<T,4>(c,0,-s,0,0,1,0,0,s,0,c,0,0,0,0,1), g);
}
Geom* g_zrot(Geom* d, Geom* g) {
  T a = degrees_to_radians(g_num_val(d));
  T c=cos(a),s=sin(a);
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
    return g_poly(difference(g_poly_val(a), g_poly_val(b)));
  else {
    error("Bad args for difference"); return NULL;
  }
}
Geom* g_not(Geom* a) {
  if (a->k == mesh_kind)
    return new MeshGeom(invert_mesh(g_mesh_val(a)));
  else if (is_poly(a))
    return g_poly(invert_poly(g_poly_val(a)));
  else {
    error("Bad args for not"); return NULL;
  }
}
Geom* g_offset(Geom* a, Geom* g) {
  if (g->k == mesh_kind) {
    // return new MeshGeom(offset_mesh(1, g_num_val(a), g_mesh_val(g)));
    return new MeshGeom(offset_mesh(g_num_val(a), g_mesh_val(g)));
  } else if (g->k == nested_v2d_kind) {
    return g_poly(offset_polyline(16, g_num_val(a), g_nested_v2d_val(g)));
  } else if (is_poly(g)) {
    return g_poly(offset_poly(16, g_num_val(a), g_poly_val(g)));
  } else {
    error("Bad args for offset"); return NULL;
  }
}
Geom* g_hollow(Geom* a, Geom* m) { return g_difference(m, g_offset(g_num(-g_num_val(a)), m)); }
Geom* g_simplify(Geom* g) { return new MeshGeom(simplify_mesh(g_mesh_val(g))); }
Geom* g_cleanup(Geom* g) { return new MeshGeom(cleanup_mesh(g_mesh_val(g))); }
Geom* g_slice(Geom* a, Geom* g) { return new PolyGeom(slice(g_num_val(a), g_mesh_val(g))); }
Geom* g_extrude(Geom* a, Geom* p) { return new MeshGeom(extrude(g_num_val(a), g_poly_val(p))); }
Geom* g_thicken(Geom* a, Geom* l) {
  if (is_nested_v2d(l))
    return new PolyGeom(thicken(1, g_num_val(a), g_nested_v2d_val(l)));
  else //  (is_nested_v3d(l))
    return new MeshGeom(thicken(1, g_num_val(a), g_nested_v3d_val(l)));
}

Geom* g_sphere(Geom* a) { return new MeshGeom(sphere_mesh(1, vec(0.0, 0.0, 0.0), g_num_val(a))); }
Geom* g_cube(Geom* a) { auto r = g_num_val(a); return new MeshGeom(cube_mesh(vec(-r, -r, -r), vec(r, r, r))); }
Geom* g_cube_lo_hi(Geom* lo, Geom* hi) { return new MeshGeom(cube_mesh(g_v3d_val(lo), g_v3d_val(hi))); }
Geom* g_cone(Geom* a, Geom* p) { return new MeshGeom(cone_mesh(g_num_val(a), g_poly_val(p))); }
Geom* g_revolve(Geom* p) { return new MeshGeom(revolve(16, g_poly_val(p))); }
Geom* g_hull(Geom* g) {
  if (g->k == mesh_kind)
    return new MeshGeom(quick_hull_mesh(g_mesh_val(g)));
  else if (g->k == poly_kind)
    return new PolyGeom(quick_hull_poly(g_poly_val(g)));
  else {
    error("Bad HULL type"); return NULL;
  }
}
// Geom* g_shear_x_z(Geom* z0, Geom* z1, Geom* dx0, Geom* dx1, Geom* m) {
//   return new MeshGeom(shear_x_z(g_num_val(z0), g_num_val(z1), g_num_val(dx0), g_num_val(dx1), g_mesh(m))); }
Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p) {
  return new MeshGeom(taper_mesh(g_num_val(l), g_num_val(r0), g_num_val(r1), g_poly_val(p))); }
