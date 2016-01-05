#ifndef __IS_GEOM_INTERFACE__
#define __IS_GEOM_INTERFACE__

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
  mesh_kind,
  expr_kind,
  octree_kind,
  meshy_kind
};

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

#endif
