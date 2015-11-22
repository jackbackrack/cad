#ifndef __IS_GEOM__
#define __IS_GEOM__

enum GeomKind {
  args_kind,
  num_kind,
  string_kind,
  v2d_kind,
  v3d_kind,
  v3i_kind,
  mat_kind,
  array_v2d_kind,
  array_v3d_kind,
  array_v3i_kind,
  nested_v2d_kind,
  nested_v3d_kind,
  poly_kind,
  mesh_kind };

class Geom;

typedef Geom* (*zero_geom_op_t)(void);
typedef Geom* (*unary_geom_op_t)(Geom*);
typedef Geom* (*binary_geom_op_t)(Geom*, Geom*);
typedef Geom* (*triple_geom_op_t)(Geom*, Geom*, Geom*);
typedef Geom* (*quad_geom_op_t)(Geom*, Geom*, Geom*, Geom*);
typedef Geom* (*quint_geom_op_t)(Geom*, Geom*, Geom*, Geom*, Geom*);
typedef Geom* (*hex_geom_op_t)(Geom*, Geom*, Geom*, Geom*, Geom*, Geom*);
typedef Geom* (*seven_geom_op_t)(Geom*, Geom*, Geom*, Geom*, Geom*, Geom*, Geom*);

class Geom {
 public:
  GeomKind k;
 Geom(GeomKind k) : k(k) { }
};

class ArgsGeom : public Geom {
 public:
  std::vector<Geom*> val;
 ArgsGeom(void) : Geom(args_kind) { }
};

class NumGeom : public Geom {
 public:
  double val;
 NumGeom(double val) : Geom(num_kind), val(val) { }
};

class StringGeom : public Geom {
 public:
  std::string val;
 StringGeom(std::string val) : Geom(string_kind), val(val) { }
};

class V2dGeom : public Geom {
 public:
  TV2 val;
 V2dGeom(TV2 val) : Geom(v2d_kind), val(val) { }
};

class V3dGeom : public Geom {
 public:
  TV3 val;
 V3dGeom(TV3 val) : Geom(v3d_kind), val(val) { }
};

class V3iGeom : public Geom {
 public:
  IV3 val;
 V3iGeom(IV3 val) : Geom(v3i_kind), val(val) { }
};

class MatGeom : public Geom {
 public:
  Matrix<T,4> val;
 MatGeom(Matrix<T,4> val) : Geom(mat_kind), val(val) { }
};

class ArrayV2dGeom : public Geom {
 public:
  Array<TV2> val;
 ArrayV2dGeom(Array<TV2> val) : Geom(array_v2d_kind), val(val) { }
};

class ArrayV3dGeom : public Geom {
 public:
  Array<TV3> val;
 ArrayV3dGeom(Array<TV3> val) : Geom(array_v3d_kind), val(val) { }
};

class ArrayV3iGeom : public Geom {
 public:
  Array<IV3> val;
 ArrayV3iGeom(Array<IV3> val) : Geom(array_v3i_kind), val(val) { }
};

class NestedV2dGeom : public Geom {
 public:
  Nested<TV2> val;
 NestedV2dGeom(Nested<TV2> val) : Geom(nested_v2d_kind), val(val) { }
};

class NestedV3dGeom : public Geom {
 public:
  Nested<TV3> val;
 NestedV3dGeom(Nested<TV3> val) : Geom(nested_v3d_kind), val(val) { }
};

class PolyGeom : public Geom {
 public:
  Nested<TV2> val;
 PolyGeom(Nested<TV2> val) : Geom(poly_kind), val(val) { }
};

class MeshGeom : public Geom {
 public:
  Mesh val;
 MeshGeom(Tuple<Ref<TriangleSoup>, Array<TV3>> val) : Geom(mesh_kind), val(const_mesh(val)) { }
 MeshGeom(Mesh val) : Geom(mesh_kind), val(val) { }
};

extern "C" int g_kind(Geom* g);

extern bool all_args_kind (std::vector<Geom*> args, int kind);
extern Geom* g_array_v2d (std::vector<Geom*> args);
extern Geom* g_array_v3d (std::vector<Geom*> args);
extern Geom* g_array_v3i (std::vector<Geom*> args);
extern Geom* g_nested_v2d (std::vector<Geom*> args);
extern Geom* g_nested_v3d (std::vector<Geom*> args);
extern Geom* g_poly (std::vector<Geom*> args);

extern "C" Geom* g_args_fab(void);
extern "C" Geom* g_args_add(Geom* g, Geom* e);
extern "C" int g_args_len(Geom* g);
extern "C" Geom* g_num_fab(int a);
extern "C" int g_num_value(Geom* g);
extern "C" Geom* g_string_fab(char* str);
extern "C" int g_string_len(Geom* g);
extern "C" char* g_string_c_str(Geom*);
extern "C" Geom* g_v2d_fab(int x, int y);
extern "C" int g_v2d_elt(Geom* g, int idx);
extern "C" int g_v2d_x(Geom* g);
extern "C" int g_v2d_y(Geom* g);
extern "C" Geom* g_v3d_fab(int x, int y, int z);
extern "C" int g_v3d_elt(Geom* g, int idx);
extern "C" int g_v3d_x(Geom* g);
extern "C" int g_v3d_y(Geom* g);
extern "C" int g_v3d_z(Geom* g);
extern "C" Geom* g_v3i_fab(int x, int y, int z);
extern "C" int g_v3i_elt(Geom* g, int idx);
extern "C" int g_v3i_x(Geom* g);
extern "C" int g_v3i_y(Geom* g);
extern "C" int g_v3i_z(Geom* g);
extern "C" Geom* g_bbox2_min(Geom* g);
extern "C" Geom* g_bbox2_max(Geom* g);
extern "C" Geom* g_bbox3_min(Geom* g);
extern "C" Geom* g_bbox3_max(Geom* g);
extern "C" Geom* g_mat_fab(int i00, int i01, int i02, int i03, int i10, int i11, int i12, int i13,
                           int i20, int i21, int i22, int i23, int i30, int i31, int i32, int i33);
extern "C" int g_mat_elt(Geom* g, int i, int j);
extern "C" Geom* g_array_v2d_fab(Geom* args);
extern "C" Geom* g_array_v2d_elt(Geom* g, int idx);
extern "C" int g_array_v2d_len(Geom* g);
extern "C" Geom* g_array_v3d_fab(Geom* args);
extern "C" Geom* g_array_v3d_elt(Geom* g, int idx);
extern "C" int g_array_v3d_len(Geom* g);
extern "C" Geom* g_array_v3i_fab(Geom* args);
extern "C" Geom* g_array_v3i_elt(Geom* g, int idx);
extern "C" int g_array_v3i_len(Geom* g);
extern "C" Geom* g_nested_v2d_fab(Geom* args);
extern "C" Geom* g_nested_v2d_elt(Geom* g, int idx);
extern "C" int g_nested_v2d_len(Geom* g);
extern "C" Geom* g_nested_v3d_fab(Geom* args);
extern "C" Geom* g_nested_v3d_elt(Geom* g, int idx);
extern "C" int g_nested_v3d_len(Geom* g);
extern "C" Geom* g_poly_fab(Geom* args);
extern "C" Geom* g_poly_elt(Geom* g, int idx);
extern "C" int g_poly_len(Geom* g);
extern "C" Geom* g_mesh_fab(Geom* points, Geom* faces);
extern "C" Geom* g_mesh_points(Geom* g);
extern "C" Geom* g_mesh_faces(Geom* g);
extern "C" Geom* g_bbox(Geom* g);
extern "C" Geom* g_dims(Geom* g);
extern "C" Geom* g_center(Geom* g);
extern "C" Geom* g_centering(Geom* g);
extern "C" Geom* g_load(Geom* s);
extern "C" Geom* g_save(Geom* s, Geom* g);
extern "C" char* g_c_str(Geom* g);
extern "C" Geom* g_to_str(Geom* g);
extern "C" Geom* g_print(Geom* g);
extern "C" Geom* g_check(Geom* g);
extern "C" Geom* g_pretty_print(Geom* g);
extern "C" Geom* g_pi(void);
extern "C" Geom* g_none2(void);
extern "C" Geom* g_all2(void);
extern "C" Geom* g_none(void);
extern "C" Geom* g_all(void);
extern "C" Geom* g_circle(Geom* a);
extern "C" Geom* g_square(Geom* a);
extern "C" Geom* g_square_lo_hi(Geom* lo, Geom* hi);
extern "C" Geom* g_letter(Geom* a);
extern "C" Geom* g_text(Geom* a);
extern "C" Geom* g_elt(Geom* g, Geom* i);
extern "C" Geom* g_mag(Geom* v, Geom* g);
extern "C" Geom* g_mag1(Geom* a, Geom* g);
extern "C" Geom* g_xmag(Geom* a, Geom* g);
extern "C" Geom* g_ymag(Geom* a, Geom* g);
extern "C" Geom* g_zmag(Geom* a, Geom* g);
extern "C" Geom* g_mov(Geom* v, Geom* g);
extern "C" Geom* g_xmov(Geom* a, Geom* g);
extern "C" Geom* g_ymov(Geom* a, Geom* g);
extern "C" Geom* g_zmov(Geom* a, Geom* g);
extern "C" Geom* g_rot(Geom* from, Geom* to, Geom* g);
extern "C" Geom* g_xrot(Geom* a, Geom* g);
extern "C" Geom* g_yrot(Geom* a, Geom* g);
extern "C" Geom* g_zrot(Geom* a, Geom* g);
extern "C" Geom* g_reflect_x(Geom* g);
extern "C" Geom* g_reflect_y(Geom* g);
extern "C" Geom* g_reflect_xy(Geom* g);
extern "C" Geom* g_reflect_z(Geom* g);
extern "C" Geom* g_reflect_yz(Geom* g);
extern "C" Geom* g_reflect_xz(Geom* g);
extern "C" Geom* g_add(Geom* a, Geom* b);
extern "C" Geom* g_mul(Geom* a, Geom* b);
extern "C" Geom* g_div(Geom* a, Geom* b);
extern "C" Geom* g_sub(Geom* a, Geom* b);
extern "C" Geom* g_intersection(Geom* a, Geom* b);
extern "C" Geom* g_union(Geom* a, Geom* b);
extern "C" Geom* g_difference(Geom* a, Geom* b);
extern "C" Geom* g_not(Geom* a);
extern "C" Geom* g_offset(Geom* a, Geom* g);
extern "C" Geom* g_hollow(Geom* a, Geom* m);
extern "C" Geom* g_simplify(Geom* g);
extern "C" Geom* g_slice(Geom* a, Geom* g);
extern "C" Geom* g_extrude(Geom* a, Geom* p);
extern "C" Geom* g_thicken(Geom* a, Geom* l);
extern "C" Geom* g_sphere(Geom* a);
extern "C" Geom* g_cube(Geom* a);
extern "C" Geom* g_cube_lo_hi(Geom* lo, Geom* hi);
extern "C" Geom* g_cone(Geom* a, Geom* p);
extern "C" Geom* g_revolve(Geom* p);
extern "C" Geom* g_hull(Geom* m);
extern "C" Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p);

extern std::vector<Geom*> g_args_val(Geom* g);
extern Geom* g_num(T a);
extern T g_num_val(Geom* g);
extern Geom* g_string(std::string s);
extern std::string g_string_val(Geom* g);

extern Geom* g_v2d(TV2 v);
extern TV2 g_v2d_val(Geom* g);

extern Geom* g_v3d(TV3 v);
extern TV3 g_v3d_val(Geom* g);

extern Geom* g_v3i(IV3 v);
extern IV3 g_v3i_val(Geom* g);

extern Geom* g_mat(Matrix<T,4> mat);
extern Matrix<T,4> g_mat_val(Geom* g);

extern Geom* g_array_v2d(Array<TV2> line);
extern Geom* g_array_v2d(RawArray<TV2> line);
extern Array<TV2> g_array_v2d_val(Geom* g);

extern Geom* g_array_v3d(Array<TV3> line);
extern Geom* g_array_v3d(RawArray<TV3> line);
extern Array<TV3> g_array_v3d_val(Geom* g);

extern Geom* g_array_v3i(Array<IV3> a);
extern Array<IV3> g_array_v3i_val(Geom* g);

extern bool is_nested_v2d(Geom* g);
extern Geom* g_nested_v2d(Nested<TV2> polyline);
extern Nested<TV2> g_nested_v2d_val(Geom* g);

extern Geom* g_nested_v3d(Nested<TV3> polyline);
extern bool is_nested_v3d(Geom* g);
extern Nested<TV3> g_nested_v3d_val(Geom* g);

extern bool is_poly(Geom* g);
extern Geom* g_poly(Nested<TV2> poly);
extern Nested<TV2> g_poly_val(Geom* g);

extern Geom* g_mesh(Mesh mesh);
extern Mesh g_mesh_val(Geom* g);
extern std::string g_to_str_val(Geom* g);


#endif
