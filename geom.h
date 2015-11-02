#ifndef __IS_GEOM__
#define __IS_GEOM__

enum GeomKind {
  float_kind,
  string_kind,
  vec2_kind,
  vec3_kind,
  mat_kind,
  line2_kind,
  line3_kind,
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
  TV val;
 Vec3Geom(TV val) : Geom(vec3_kind), val(val) { }
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
  Array<TV> val;
 Line3Geom(Array<TV> val) : Geom(line3_kind), val(val) { }
};

class PolyLine2Geom : public Geom {
 public:
  Nested<TV2> val;
 PolyLine2Geom(Nested<TV2> val) : Geom(polyline2_kind), val(val) { }
};

class PolyLine3Geom : public Geom {
 public:
  Nested<TV> val;
 PolyLine3Geom(Nested<TV> val) : Geom(polyline3_kind), val(val) { }
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
 MeshGeom(Tuple<Ref<TriangleSoup>, Array<TV>> val) : Geom(mesh_kind), val(const_soup(val)) { }
 MeshGeom(Mesh val) : Geom(mesh_kind), val(val) { }
};

extern double g_val(Geom* g);
extern std::string g_string(Geom* g);
extern TV2 g_vec2(Geom* g);
extern TV g_vec3(Geom* g);
extern TV g_vec(Geom* g);
extern Matrix<T,4> g_mat(Geom* g);
extern Array<TV2> g_line2(Geom* g);
extern Array<TV> g_line3(Geom* g);
extern Array<IV> g_faces(Geom* g);
extern bool is_polyline2(Geom* g);
extern Nested<TV2> g_polyline2(Geom* g);
extern bool is_polyline3(Geom* g);
extern Nested<TV> g_polyline3(Geom* g);
extern Array<TV2> g_contour(Geom* g);
extern bool is_poly(Geom* g);
extern Nested<TV2> g_poly(Geom* g);
extern Mesh g_mesh(Geom* g);
extern Geom* g_load(Geom* s);
extern Geom* g_save(Geom* s, Geom* g);
extern Geom* g_print(Geom* g);
extern Geom* g_pretty_print(Geom* g);
extern Geom* g_float(double a);
extern Geom* g_pi(void);
extern Geom* g_none2(void);
extern Geom* g_all2(void);
extern Geom* g_none(void);
extern Geom* g_all(void);
extern Geom* g_circle(Geom* a);
extern Geom* g_square(Geom* a);
extern Geom* g_square(Geom* lo, Geom* hi);
extern Geom* g_letter(Geom* a);
extern Geom* g_text(Geom* a);
extern Geom* g_elt(Geom* g, Geom* i);
extern Geom* g_mag(Geom* v, Geom* g);
extern Geom* g_mag1(Geom* a, Geom* g);
extern Geom* g_xmag(Geom* a, Geom* g);
extern Geom* g_ymag(Geom* a, Geom* g);
extern Geom* g_zmag(Geom* a, Geom* g);
extern Geom* g_mov(Geom* v, Geom* g);
extern Geom* g_xmov(Geom* a, Geom* g);
extern Geom* g_ymov(Geom* a, Geom* g);
extern Geom* g_zmov(Geom* a, Geom* g);
extern Geom* g_rot(Geom* from, Geom* to, Geom* g);
extern Geom* g_xrot(Geom* a, Geom* g);
extern Geom* g_yrot(Geom* a, Geom* g);
extern Geom* g_zrot(Geom* a, Geom* g);
extern Geom* g_reflect_x(Geom* g);
extern Geom* g_reflect_y(Geom* g);
extern Geom* g_reflect_xy(Geom* g);
extern Geom* g_reflect_z(Geom* g);
extern Geom* g_reflect_yz(Geom* g);
extern Geom* g_reflect_xz(Geom* g);
extern Geom* g_add(Geom* a, Geom* b);
extern Geom* g_mul(Geom* a, Geom* b);
extern Geom* g_sub(Geom* a, Geom* b);
extern Geom* g_not(Geom* a);
extern Geom* g_offset(Geom* a, Geom* g);
extern Geom* g_simplify(Geom* g);
extern Geom* g_slice(Geom* a, Geom* g);
extern Geom* g_extrude(Geom* a, Geom* p);
extern Geom* g_thicken(Geom* a, Geom* l);
extern Geom* g_sphere(Geom* a);
extern Geom* g_cube(Geom* a);
extern Geom* g_cube(Geom* lo, Geom* hi);
extern Geom* g_cone(Geom* a, Geom* p);
extern Geom* g_revolve(Geom* p);
// extern Geom* g_hollow(Geom* a, Geom* m);
extern Geom* g_taper(Geom* l, Geom* r0, Geom* r1, Geom* p);

#endif