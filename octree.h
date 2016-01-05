#ifndef __IS_OCTREE__
#define __IS_OCTREE__

#include "cad.h"
#include "expr.h"
#include "path.h"

class Meshy : public Geom {
 public:
  Array<TV3> vdata;
  Array<TV3> gdata;
  Array<IV3> tdata;
 Meshy(void) : Geom(meshy_kind) { }
};

const flo_t EPSILON = 1e-6;

enum voxel_kind_t {
  voxel_filled_kind,
  voxel_empty_kind,
  voxel_leaf_kind,
  voxel_parent_kind
};

class Octree : public Boxy, public Geom {
 public: 
  Octree*  parent;
  Fun*     g;
  union {
    flo_t      d[8];
    Octree*    children[8];
  } u;
  voxel_kind_t kind;
  union {
    uint32_t* tri;
    Path** contour;
  } data;
 Octree(Fun* g, Octree* parent, const Vec &lo, const Vec &hi, voxel_kind_t kind /*, flo_t *dists */) : Boxy(lo, hi), Geom(octree_kind), parent(parent), g(g), kind(kind) {
    data.tri = NULL;
    k = octree_kind;
    // for (int i = 0; i < 8; i++)
    //   d[i] = dists[i];
  }
  Octree* lookup(const Vec &p);
  Octree* leaf(const Vec &p);
  void nbrs (const Boxy &box, std::vector< Octree* >& nbrs);
  std::vector< Octree* > all_nbrs (void);
  void boundaries (const Boxy &box, std::vector< Octree* >& leaves);
  void boundary_voxels(std::vector<Octree*> &voxels);
  std::vector<Tri> boundary_voxel_triangles(flo_t r);
  inline flo_t dget(int i, int j, int k) const { return u.d[(i<<2) | (j<<1) | k]; }
  inline void dset(flo_t v, int i, int j, int k) { u.d[(i<<2) | (j<<1) | k] = v; }
  inline flo_t rad( void ) const { return hi.x - lo.x; }
  inline int child_index (const Vec &p) const {
    auto ctr = center();
    int dx = (p.x - ctr.x) >= 0.0;
    int dy = (p.y - ctr.y) >= 0.0;
    int dz = (p.z - ctr.z) >= 0.0;
    int q  = dz * 4 + dy * 2 + dx;
    return q;
  }
  long count_nodes (void) {
    long gc = 1;
    if (kind == voxel_parent_kind) {
      for (int i = 0; i < 8; i++) 
        gc += u.children[i]->count_nodes();
    }
    return gc;
  }
  long count_fun (void) {
    long gc = g->size();
    if (kind == voxel_parent_kind) {
      for (int i = 0; i < 8; i++) 
        gc += u.children[i]->count_fun();
    }
    return gc;
  }
  flo_t interpolate(const Vec &p);
  flo_t do_interpolate(const Vec &p);
  inline flo_t interpolate_one(const Vec &p) const {
    Vec diam = hi - lo;
    Vec del = (p - lo) / diam;
    flo_t d00 = dget(0,0,0) * (1 - del.x) + dget(1,0,0) * del.x;
    flo_t d10 = dget(0,1,0) * (1 - del.x) + dget(1,1,0) * del.x;
    flo_t d01 = dget(0,0,1) * (1 - del.x) + dget(1,0,1) * del.x;
    flo_t d11 = dget(0,1,1) * (1 - del.x) + dget(1,1,1) * del.x;
    flo_t d0  = d00 * (1 - del.y) + d10 * del.y;
    flo_t d1  = d01 * (1 - del.y) + d11 * del.y;
    flo_t d   = d0 * (1 - del.z) + d1 * del.z;
    return d;
  }
  Vec gradient_one(const Vec &p) const {
    const flo_t dx    = hi.x - lo.x,
                dy    = hi.y - lo.y,
                dz    = hi.z - lo.z;

    const flo_t xd    = dx ? (p.x - lo.x) / dx : 0.5,
                yd    = dy ? (p.y - lo.y) / dy : 0.5,
                zd    = dz ? (p.z - lo.z) / dz : 0.5;

    const flo_t c000  = dget(0,0,0), c001 = dget(0,0,1), c010  = dget(0,1,0), c011 = dget(0,1,1), 
                c100  = dget(1,0,0), c101 = dget(1,0,1), c110  = dget(1,1,0), c111 = dget(1,1,1);

    // Numerically find the gradients
    return vec(-(zd - 1)*((yd - 1)*(c000/dx - c100/dx) - (c010/dx - c110/dx)*yd) + (p.z - lo.z)*((yd - 1)*(c001/dx - c101/dx) - (c011/dx - c111/dx)*yd)/dz,

                (zd - 1)*(((xd - 1)*c010 - xd*c110)/dy - ((xd - 1)*c000 - xd*c100)/dy) - (((xd - 1)*c011 - xd*c111)/dy - ((xd - 1)*c001 - xd*c101)/dy)*zd,

               -((yd - 1)*((xd - 1)*c000 - xd*c100) - (p.y - lo.y)*((xd - 1)*c010 - xd*c110)/dy)/dz + ((yd - 1)*((xd - 1)*c001 - xd*c101) - (p.y - lo.y)*((xd - 1)*c011 - xd*c111)/dy)/dz);
  }
  Vec zero_crossing(const uint8_t v0, const uint8_t v1) {
    const flo_t x0 = (v0 & 4) ? hi.x : lo.x,
                y0 = (v0 & 2) ? hi.y : lo.y,
                z0 = (v0 & 1) ? hi.z : lo.z;

    const flo_t dx = ((v1 & 4) ? hi.x : lo.x) - x0,
                dy = ((v1 & 2) ? hi.y : lo.y) - y0,
                dz = ((v1 & 1) ? hi.z : lo.z) - z0;


    flo_t d = 0.5;
    flo_t step = 0.25;

    // Binary search along the edge to find the zero crossing
    // (because math is hard and for loops are easy)
    for (int iteration = 0; iteration < 16; ++iteration) {

      // Use interpolation to find the distance metric value at this edge
      const flo_t result = interpolate_one(vec(x0 + d*dx, y0 + d*dy, z0 + d*dz));

      // Change the binary search coefficient accordingly
      if (result < -EPSILON)         d += step;
      else if (result > EPSILON)     d -= step;
      else                        break;

      // And shrink the step size
      step /= 2;
    }

    return vec(x0+d*dx, y0+d*dy, z0+d*dz);
  }

  void get_neighbors_3d(const Octree* const old_nbrs[6], const Octree* new_nbrs[6], uint8_t b);
  void get_neighbors_2d(const Octree* const old_nbrs[4], const Octree* new_nbrs[4], uint8_t b);
  Vec zero_crossing(const Vec &p, const Vec &v) const {
    flo_t d = 0.5;
    flo_t step = 0.25;
    // printf("I ");
    for (int i = 0; i < 24; i++) {
      auto sp  = p + v * d;
      auto res = interpolate_one(sp);
      // printf("%f ", res);
      if (res < -EPSILON)     d += step;
      else if (res > EPSILON) d -= step;
      else { /* printf("SUCCESS\n"); */ return sp; }
      step *= 0.5;
    }
    // auto sp = p + v * d;
    printf("FAILED\n");
    auto sp = vec(INFTY,INFTY,INFTY);
    return sp;
  }
  bool has_x_nbr ( void ) const {
    for (int z = 0; z < 2; z++) 
      for (int y = 0; y < 2; y++)
        if (dget(1,y,z) < 0) return true;
    return false;
  }
  bool has_y_nbr ( void ) const {
    for (int x = 0; x < 2; x++)
      for (int z = 0; z < 2; z++)
        if (dget(x,1,z) < 0) return true;
    return false;
  }
  bool has_z_nbr ( void ) const {
    for (int x = 0; x < 2; x++)
      for (int y = 0; y < 2; y++)
        if (dget(x,y,1) < 0) return true;
    return false;
  }

  Vec boundary_point(const Vec &p) const {
    auto g = gradient_one(p);
    return zero_crossing(p, g * rad() * 2.0);
  }
  Segment<TV3> find_crossing_line ( void ) const {
    for (int x = 0; x < 2; x++) {
      for (int y = 0; y < 2; y++) {
        for (int z = 0; z < 2; z++) {
          if ((dget(x,y,z)<0) != (dget(1-x,1-y,1-z)<0)) {
            // printf("SUCCESS(%d,%d,%d): %d != %d\n", x, y, z, dget(x,y,z) < 0, dget(1-x,1-y,1-z) < 0);
            auto p0 = vec(x?lo.x:hi.x, y?lo.y:hi.y, z?lo.z:hi.z);
            auto p1 = vec(x?hi.x:lo.x, y?hi.y:lo.y, z?hi.z:lo.z);
            if (dget(x,y,z) < 0)
              return Segment<TV3>(p1, p0);
            else 
              return Segment<TV3>(p0, p1);
          }
        }
      }
    }
    printf("FAILED: "); 
    for (int x = 0; x < 2; x++) 
      for (int y = 0; y < 2; y++) 
        for (int z = 0; z < 2; z++) 
          printf("%d ", dget(x,y,z) < 0);
    printf("\n");
    return Segment<TV3>(vec(-INFTY,-INFTY,-INFTY),vec(INFTY,INFTY,INFTY));
  }
  Vec boundary_point(void) const {
    auto l = find_crossing_line();
    return zero_crossing(l.x0, l.x1 - l.x0);
  }
  Vec sample_point(void) const {
    return boundary_point(rnd_vec());
  }
  Segment<TV3> sample_gradient(void) const {
    auto p = rnd_vec();
    auto g = gradient_one(p);
    return Segment<TV3>(p, p + g);
  }
  void all_voxels(std::vector<Octree*> &voxels);
  std::vector<Tri> boundary_triangles(flo_t r);
  Vec get_normal(Vec pos);
  Octree* slice(const flo_t z);
  Octree* do_slice(Octree* parent, const flo_t z);
  Meshy* triangulate(volatile int* const halt);
  std::vector<Path*> contour(volatile int* const halt);
  double march_cube
    (double target_value, std::vector< Tri >& tris, std::vector< Tri >& normals, std::vector< Tri >& colors);
 private:
  void write_edges (std::vector<Path*>& paths);
  void find_edges(const Octree* const neighbors[4],volatile int* const halt);
  void evaluate_pixel(const Octree* const neighbors[4]);
  Path* contour_zero_crossing (const Octree* const neighbors[4], const uint8_t edge);
  uint32_t mesh_zero_crossing(const Octree* const neighbors[6], const uint8_t v0, const uint8_t v1, Meshy* mesh);
  void evaluate_voxel(const Octree* const neighbors[6], Meshy* mesh);
  void do_triangulate(const Octree* const neighbors[6], Meshy* mesh, volatile int* const halt);
};

extern Octree* octree(Fun* g, const Vec &lo, const Vec &hi, flo_t r, flo_t t);

extern void marching_cubes(std::vector<Octree*>& voxels, std::vector<Tri> &tris, std::vector<Tri> &colors, std::vector<Tri> &normals);
// extern std::vector<Tri> marching_cubes(std::vector<Octree*> &voxels);
extern std::vector<Tri> voxels_triangles(std::vector<Octree*> &voxels);
extern std::vector<Vec> voxels_points(std::vector<Octree*> &voxels, int n);
extern std::vector<Segment<TV3>> voxels_gradients(std::vector<Octree*> &voxels, int n);
extern std::vector< Segment<TV3>> voxels_boundary_lines(std::vector< Octree* > &voxels);

extern Nested< TV2 > tree_slice (Octree* tree, T z);
extern Octree* g_octree_val(Geom* g);
extern Octree* expr_to_tree(Expr* g, T rad, T thresh);
extern Mesh octree_to_mesh(Octree* g);
extern Octree* pretty_print_tree(Octree* g);

#endif
