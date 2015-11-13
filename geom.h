#ifndef __IS_GEOM__
#define __IS_GEOM__

enum GeomKind {
  args_kind,
  float_kind,
  string_kind,
  vec2_kind,
  vec3_kind,
  mat_kind,
  line2_kind,
  line3_kind,
  faces_kind,
  polyline2_kind,
  polyline3_kind,
  contour_kind,
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
 ArgsGeom(void) : Geom(float_kind) { }
};

class FloatGeom : public Geom {
 public:
  double val;
 FloatGeom(double val) : Geom(float_kind), val(val) { }
};

class StringGeom : public Geom {
 public:
  std::string val;
 StringGeom(std::string val) : Geom(string_kind), val(val) { }
};

class Vec2Geom : public Geom {
 public:
  TV2 val;
 Vec2Geom(TV2 val) : Geom(vec2_kind), val(val) { }
};

class Vec3Geom : public Geom {
 public:
  TV3 val;
 Vec3Geom(TV3 val) : Geom(vec3_kind), val(val) { }
};

class MatGeom : public Geom {
 public:
  Matrix<T,4> val;
 MatGeom(Matrix<T,4> val) : Geom(mat_kind), val(val) { }
};

class Line2Geom : public Geom {
 public:
  Array<TV2> val;
 Line2Geom(Array<TV2> val) : Geom(line2_kind), val(val) { }
};

class Line3Geom : public Geom {
 public:
  Array<TV3> val;
 Line3Geom(Array<TV3> val) : Geom(line3_kind), val(val) { }
};

class FacesGeom : public Geom {
 public:
  Array<IV3> val;
 FacesGeom(Array<IV3> val) : Geom(faces_kind), val(val) { }
};

class PolyLine2Geom : public Geom {
 public:
  Nested<TV2> val;
 PolyLine2Geom(Nested<TV2> val) : Geom(polyline2_kind), val(val) { }
};

class PolyLine3Geom : public Geom {
 public:
  Nested<TV3> val;
 PolyLine3Geom(Nested<TV3> val) : Geom(polyline3_kind), val(val) { }
};

class ContourGeom : public Geom {
 public:
  Array<TV2> val;
 ContourGeom(Array<TV2> val) : Geom(contour_kind), val(val) { }
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

extern std::vector<Geom*> g_args_val(Geom* g);
extern "C" Geom* g_args_fab(void);
extern "C" Geom* g_args_add(Geom* g);
extern "C" int g_args_len(Geom* g);

extern "C" Geom* g_num(T a);
extern T g_num_val(Geom* g);

extern Geom* g_string(std::string s);
extern std::string g_string_val(Geom* g);
extern "C" Geom* g_string_fab(char* str);
extern "C" int g_string_len(Geom* g);
extern "C" char* g_string_c_str(Geom*);

extern Geom* g_vec2(TV2 v);
extern TV2 g_vec2_val(Geom* g);
extern "C" Geom* g_vec2_fab(T x, T y);
extern "C" T g_vec2_elt(Geom* g, int idx);
extern "C" T g_vec2_x(Geom* g);
extern "C" T g_vec2_y(Geom* g);

extern Geom* g_vec3(TV3 v);
extern TV3 g_vec3_val(Geom* g);
extern "C" Geom* g_vec3_fab(T x, T y, T z);
extern "C" T g_vec3_elt(Geom* g, int idx);
extern "C" T g_vec3_x(Geom* g);
extern "C" T g_vec3_y(Geom* g);
extern "C" T g_vec3_z(Geom* g);

extern "C" Geom* g_bbox2_min(Geom* g);
extern "C" Geom* g_bbox2_max(Geom* g);
extern "C" Geom* g_bbox3_min(Geom* g);
extern "C" Geom* g_bbox3_max(Geom* g);

extern Geom* g_mat(Matrix<T,4> mat);
extern Matrix<T,4> g_mat_val(Geom* g);
extern "C" Geom* g_mat_fab(T i00, T i01, T i02, T i03, T i10, T i11, T i12, T i13,
                           T i20, T i21, T i22, T i23, T i30, T i31, T i32, T i33);
extern "C" T g_mat_fab_elt(Geom* g, int idx);

extern Geom* g_line2(Array<TV2> line);
extern Geom* g_line2(RawArray<TV2> line);
extern Array<TV2> g_line2_val(Geom* g);
extern "C" Geom* g_line2_fab(Geom* args);
extern "C" Geom* g_line2_elt(Geom* g, int idx);
extern "C" int g_line2_fab_len(Geom* g);

extern Geom* g_line3(Array<TV3> line);
extern Geom* g_line3(RawArray<TV3> line);
extern Array<TV3> g_line3_val(Geom* g);
extern "C" Geom* g_line3_fab(Geom* args);
extern "C" Geom* g_line3_elt(Geom* g, int idx);
extern "C" int g_line3_fab_len(Geom* g);

extern Geom* g_faces(Array<IV3> faces);
extern Array<IV3> g_faces_val(Geom* g);
extern "C" Geom* g_faces_fab(Geom* args);
extern "C" Geom* g_faces_elt(Geom* g, int idx);
extern "C" int g_faces_fab_len(Geom* g);

extern bool is_polyline2(Geom* g);
extern Geom* g_polyline2(Nested<TV2> polyline);
extern Nested<TV2> g_polyline2_val(Geom* g);
extern "C" Geom* g_polyline2_fab(Geom* args);
extern "C" Geom* g_polyline2_elt(Geom* g, int idx);
extern "C" int g_polyline2_fab_len(Geom* g);

extern Geom* g_polyline3(Nested<TV3> polyline);
extern bool is_polyline3(Geom* g);
extern Nested<TV3> g_polyline3_val(Geom* g);
extern "C" Geom* g_polyline3_fab(Geom* args);
extern "C" Geom* g_polyline3_elt(Geom* g, int idx);
extern "C" int g_polyline3_fab_len(Geom* g);

extern Geom* g_contour(Array<TV2> contour);
extern Geom* g_contour(RawArray<TV2> contour);
extern Array<TV2> g_contour_val(Geom* g);
extern "C" Geom* g_contour_fab(Geom* args);
extern "C" Geom* g_contour_elt(Geom* g, int idx);
extern "C" int g_contour_fab_len(Geom* g);

extern Geom* g_poly(Nested<TV2> poly);
extern bool is_poly(Geom* g);
extern Nested<TV2> g_poly_val(Geom* g);
extern "C" Geom* g_poly_fab(Geom* args);
extern "C" Geom* g_poly_elt(Geom* g, int idx);
extern "C" int g_poly_fab_len(Geom* g);

extern Geom* g_mesh(Mesh mesh);
extern Mesh g_mesh_val(Geom* g);
extern "C" Geom* g_mesh_fab(Geom* points, Geom* faces);
extern "C" Geom* g_mesh_points(Geom* g);
extern "C" Geom* g_mesh_faces(Geom* g);

extern "C" Geom* g_bbox(Geom* g);

extern "C" Geom* g_load(Geom* s);
extern "C" Geom* g_save(Geom* s, Geom* g);
extern "C" char* g_c_str(Geom* g);
extern "C" Geom* g_to_str(Geom* g);
extern std::string g_to_str_val(Geom* g);
extern "C" Geom* g_print(Geom* g);
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
extern "C" Geom* g_sub(Geom* a, Geom* b);
extern "C" Geom* g_intersection(Geom* a, Geom* b);
extern "C" Geom* g_union(Geom* a, Geom* b);
extern "C" Geom* g_difference(Geom* a, Geom* b);
extern "C" Geom* g_not(Geom* a);
extern "C" Geom* g_offset(Geom* a, Geom* g);
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
// extern "C" Geom* g_hollow(Geom* a, Geom* m);
extern "C" Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p);

#endif
