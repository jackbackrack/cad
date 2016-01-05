#include "cad.h"
#include "octree.h"
// #include "glut-c.h"
// #include "c-stub.h"

const bool is_print = false;
extern bool is_interpolating;

inline Octree* mark_leaf(Octree* tree, bool is_empty, bool is_filled, bool is_reproducible) {
  if (is_print) printf("REPRODUCIBLE %d\n", is_reproducible);
  tree->kind = is_empty ? voxel_empty_kind : (is_filled ? voxel_filled_kind : voxel_leaf_kind);
  return tree;
}

extern bool is_epsilon_target;

Octree* do_octree(Octree* parent, Fun* g, const Vec &lo, const Vec &hi, const flo_t r, flo_t t, int d) {
  //  Octree* octree(Octree* parent, Geom* g, const Vec &p, const flo_t r, flo_t t, int d) {
  const Vec p = (lo + hi) * 0.5;
  auto bounds = box(lo, hi);
  // auto i = g->range(bounds);
  // printf("PRUNING ");
  // auto ng = g->is_prune() ? g->prune(d) : g;
  Interval i;
  auto ng = g->maybe_prune(d, i, bounds);
  // if (ng == g) printf("-->NO\n"); else printf("-->YES\n");
  // auto ng = g;
  if (is_print) for (int i = 0; i < d; i++) printf("  ");
  if (is_print) printf("VOX %f,%f,%f I %f,%f R %f ", p.x, p.y, p.z, i.lo, i.hi, r);
  if (i.lo >= 0.0) {
    if (is_print) printf("EMPTY\n");
    return new Octree(ng, parent, lo, hi, voxel_empty_kind);
  } else if (i.hi < 0.0) {
    if (is_print) printf("FILLED\n");
    return new Octree(ng, parent, lo, hi, voxel_filled_kind);
  } else {
    auto res = new Octree(ng, parent, lo, hi, voxel_parent_kind);
    bool is_empty        = true;
    bool is_filled       = true;
    int i = 0;
    // printf("GEOM %s\n", ng->to_str().c_str());
    for (int dx = -1; dx <= 1; dx += 2) {
      for (int dy = -1; dy <= 1; dy += 2) {
        for (int dz = -1; dz <= 1; dz += 2) {
          auto tp = p + vec((T)dx,(T)dy,(T)dz) * r;
          auto dv = res->u.d[i++] = ng->dist(tp);
          if (dv < 0)
            is_empty = false;
          else
            is_filled = false;
        }
      }
    }
    if (r <= t)
      return mark_leaf(res, is_empty, is_filled, false);
    if (is_interpolating) {
      bool is_reproducible = true;
      auto target = is_epsilon_target ? EPSILON : t; // TODO: > t EPSILON
      for (int dx = -1; dx <= 1; dx += 1) {
        for (int dy = -1; dy <= 1; dy += 1) {
          for (int dz = -1; dz <= 1; dz += 1) {
            if (dx == 0 || dy == 0 || dz == 0) {
              auto tp = p + vec((T)dx,(T)dy,(T)dz) * r;
              auto cd = ng->dist(tp);
              auto id = res->interpolate_one(tp);
              if (fabs(cd - id) > target || ((id < 0) != (cd < 0))) 
                is_reproducible = false;
            }
          }
        }
      }
      if (is_reproducible)
        return mark_leaf(res, is_empty, is_filled, true);
    } 
    if (is_print) printf("INTERNAL\n");
    auto r2 = r * 0.5;
    int k = 0;
    for (int dx = -1; dx <= 1; dx += 2) {
      for (int dy = -1; dy <= 1; dy += 2) {
        for (int dz = -1; dz <= 1; dz += 2) {
          const Vec c_lo(dx<0?lo.x:p.x, dy<0?lo.y:p.y, dz<0?lo.z:p.z);
          const Vec c_hi(dx<0?p.x:hi.x, dy<0?p.y:hi.y, dz<0?p.z:hi.z);
          res->u.children[k++] = do_octree(res, ng, c_lo, c_hi, r2, t, d + 1);
        }
      }
    }
    return res;
  }
}

Octree* octree(Fun* g, const Vec &lo, const Vec &hi, const flo_t r, flo_t t) {
  return do_octree(NULL, g, lo, hi, r, t, 0);
}

Octree* Octree::do_slice(Octree* parent, const flo_t z) {
  auto out = new Octree(g, parent, vec(lo.x, lo.y, z), vec(hi.x, hi.y, z), voxel_parent_kind);

  if (kind == voxel_leaf_kind) {
    out->kind = voxel_leaf_kind;
    for (int i=0; i < 8; ++i) {
      out->u.d[i] = interpolate_one(vec(i & 4 ? hi.x : lo.x, i & 2 ? hi.y : lo.y, z));
    }
  } else if (kind == voxel_filled_kind) {
    out->kind = voxel_filled_kind;
  } else if (kind == voxel_empty_kind) {
    out->kind = voxel_empty_kind;
  } else if (kind == voxel_parent_kind) {
    uint8_t upper = u.children[1] && z > u.children[1]->lo.z;
    for (int i=0; i < 8; i += 1) {
      out->u.children[i] = i&1 ? NULL : u.children[i+upper]->do_slice(out, z);
    }
  }

  return out;
}

Octree* Octree::slice(const flo_t z) {
  return do_slice(NULL, z);
}

Octree* Octree::lookup(const Vec &p) {
  if (is_inside(p))
    return leaf(p);
  else
    return NULL;
}

Octree* Octree::leaf(const Vec &p) {
  if (kind == voxel_parent_kind) 
    return u.children[child_index(p)]->leaf(p);
  else
    return this;
}

static flo_t max_double = 1e9;

flo_t Octree::interpolate(const Vec &p) {
  if (is_inside(p))
    return do_interpolate(p);
  else
    return max_double;
}

flo_t Octree::do_interpolate(const Vec &p) {
  if (kind == voxel_parent_kind) 
    return u.children[child_index(p)]->do_interpolate(p);
  else
    return interpolate_one(p);
}

void Octree::all_voxels(std::vector<Octree*> &voxels) {
  voxels.push_back(this);
  if (kind == voxel_parent_kind)
    for (int i = 0; i < 8; i++)
      u.children[i]->all_voxels(voxels);
}

void Octree::boundary_voxels(std::vector<Octree*> &voxels) {
  if (kind == voxel_leaf_kind)
    voxels.push_back(this);
  if (kind == voxel_parent_kind)
    for (int i = 0; i < 8; i++)
      u.children[i]->boundary_voxels(voxels);
}

void two_tris(const Vec &p00, const Vec &p01, const Vec &p10, const Vec &p11, std::vector<Tri> &tris) {
  tris.push_back(tri(p00, p10, p11));
  tris.push_back(tri(p00, p11, p01));
}

std::vector<Tri> voxels_triangles(std::vector<Octree*> &voxels) {
  std::vector<Tri> res;
  for (auto v : voxels) {
    const Vec l = v->lo;
    const Vec h = v->hi;
    two_tris(vec(l.x,l.y,l.z),vec(l.x,l.y,h.z),vec(l.x,h.y,l.z),vec(l.x,h.y,h.z),res);
    two_tris(vec(l.x,l.y,l.z),vec(l.x,l.y,h.z),vec(h.x,l.y,l.z),vec(h.x,l.y,h.z),res);
    two_tris(vec(l.x,l.y,l.z),vec(l.x,h.y,l.z),vec(h.x,l.y,l.z),vec(h.x,h.y,l.z),res);
    two_tris(vec(h.x,l.y,l.z),vec(h.x,l.y,h.z),vec(h.x,h.y,l.z),vec(h.x,h.y,h.z),res);
    two_tris(vec(l.x,h.y,l.z),vec(l.x,h.y,h.z),vec(h.x,h.y,l.z),vec(h.x,h.y,h.z),res);
    two_tris(vec(l.x,l.y,h.z),vec(l.x,h.y,h.z),vec(h.x,l.y,h.z),vec(h.x,h.y,h.z),res);
  }
  return res;
}

void Octree::get_neighbors_3d
    (const Octree* const old_nbrs[6], const Octree* new_nbrs[6], const uint8_t b) {
  for (int i=0; i < 6; ++i) new_nbrs[i] = NULL;

  if (!u.children[b]) return;

  for (uint8_t axis=0; axis < 6; ++axis) {

    uint8_t mask = 1 << (axis/2);
    uint8_t dir  = mask * (axis % 2);

    // If we're pointing to within our own cell, then
    // pick the interior neighbor.
    if ( (b&mask)^dir && u.children[b^mask]) {
      new_nbrs[axis] = u.children[b^mask];
    }

    // Correctly handle situations where a non-branch is
    // next to a branch, but the branch only splits on
    // the parallel axis.  Ascii art:
    //
    //   |---|---:---|
    //   | N | A : B |
    //   |---|---:---|
    //
    // N is the neighbor, the larger cell is splitting into
    // A and B on the axis along which it touches N (and no
    // others axes) so N remains a valid neighbor.
    else if (old_nbrs[axis] && old_nbrs[axis]->kind == voxel_leaf_kind) {

      bool crack = false;
      for (int split=0; split < 3; ++split) {
        if ((1 << split) & mask)    continue;
        if (u.children[1<<split] != NULL) crack = true;
      }

      new_nbrs[axis] = crack ? NULL : old_nbrs[axis];
    }

    // Otherwise, check to see that the neighbor splits on
    // the same axes (other than the one along which we're
    // joining).
    else if (old_nbrs[axis] && old_nbrs[axis]->kind == voxel_parent_kind) {
      uint8_t split_mask = (old_nbrs[axis]->u.children[1] ? 1 : 0) |
        (old_nbrs[axis]->u.children[2] ? 2 : 0) |
        (old_nbrs[axis]->u.children[4] ? 4 : 0);

      bool crack = false;
      for (int split=0; split < 3; ++split) {

        // If we split differently than our neighbor and
        // this split isn't along the relevant axis, record
        // that there's a crack here.
        bool split_A = u.children[1<<split] != NULL;
        bool split_B = old_nbrs[axis]->u.children[1<<split] != NULL;

        if ((split_A ^ split_B) && !((1 << split) & mask)) {
          crack = true;
        }
      }

      // Pick the appropriate cell in the neighbor
      if (crack) {
        new_nbrs[axis] = NULL;
      } else if (dir) {
        new_nbrs[axis] = old_nbrs[axis]->u.children[(b &~mask) & split_mask];
      } else {
        new_nbrs[axis] = old_nbrs[axis]->u.children[(b | mask) & split_mask];
      }

    }
  }

}

void Octree::get_neighbors_2d
    (const Octree* const old_nbrs[4], const Octree* new_nbrs[4], uint8_t b) {
  // Zero out the array elements
  for (int i=0; i < 4; ++i)   new_nbrs[i] = NULL;

  if (!u.children[b]) return;

  for (uint8_t axis=0; axis < 4; ++axis) {

    uint8_t mask = 1 << (axis/2 + 1);
    uint8_t dir  = mask * (axis % 2);

    // If we're pointing to within our own cell, then
    // pick the interior neighbor.
    if ( (b&mask)^dir && u.children[b^mask]) {
      new_nbrs[axis] = u.children[b^mask];
    }

    else if (old_nbrs[axis] && old_nbrs[axis]->kind != voxel_parent_kind) {
      new_nbrs[axis] = old_nbrs[axis];
    }

    // Otherwise, check to see that the neighbor splits on
    // the same axes (other than the one along which we're
    // joining).
    else if (old_nbrs[axis] && old_nbrs[axis]->kind == voxel_parent_kind) {

      uint8_t split_mask = (old_nbrs[axis]->u.children[2] ? 2 : 0) |
        (old_nbrs[axis]->u.children[4] ? 4 : 0);

      bool crack = false;
      for (int split=0; split < 3; ++split) {

        // If we split differently than our neighbor and
        // this split isn't along the relevant axis, record
        // that there's a crack here.
        bool split_A = u.children[1<<split] != NULL;
        bool split_B = old_nbrs[axis]->u.children[1<<split] != NULL;

        if ((split_A ^ split_B) && !((1 << split) & mask)) {
          crack = true;
        }
      }

      if (crack) { // Multi-scale stitching is okay in 2D
        new_nbrs[axis] = old_nbrs[axis];
      } else if (dir) {
        new_nbrs[axis] = old_nbrs[axis]->u.children[(b &~mask) & split_mask];
      } else {
        new_nbrs[axis] = old_nbrs[axis]->u.children[(b | mask) & split_mask];
      }

    }
  }

}

static int rndi () {
  return rand();
}

static int rndi (int mn, int mx) {
  return (rndi() % (mx-mn+1)) + mn;
}

std::vector<Vec> voxels_points(std::vector<Octree*> &voxels, int n) {
  flo_t tot_area = 0.0;
  for (auto voxel : voxels) 
    tot_area += sqr(voxel->rad());
  std::vector<Vec> res;
  for (auto voxel : voxels) {
    auto area = sqr(voxel->rad());
    int vn = roundf(n * area / tot_area);
    for (int i = 0; i < vn; i++) {
      // res.push_back(voxel->rnd_vec());
      res.push_back(voxel->sample_point());
    }
  }
  return res;
}

std::vector< Segment<TV3> > voxels_gradients(std::vector<Octree*> &voxels, int n) {
  flo_t tot_area = 0.0;
  for (auto voxel : voxels) 
    tot_area += sqr(voxel->rad());
  std::vector< Segment<TV3> > res;
  for (auto voxel : voxels) {
    auto area = sqr(voxel->rad());
    int vn = roundf(n * area / tot_area);
    for (int i = 0; i < vn; i++) {
      res.push_back(voxel->sample_gradient());
    }
  }
  return res;
}

// inline bool is_border (int ov, int cv) {
//   return ov != -1 && ov != cv;
// }

inline bool is_border (Octree *ov, Octree *cv) {
  return ov == NULL || (ov->kind != voxel_leaf_kind && ov->kind != voxel_filled_kind);
}

void print_voxel_box(Octree* voxel, const Boxy &box) {
  /*
  printf("CHECKING %p LO %f,%f,%f HI %f,%f,%f BOX LO %f,%f,%f HI %f,%f,%f -> %d\n",
         voxel, voxel->lo.x, voxel->lo.y, voxel->lo.z, voxel->hi.x, voxel->hi.y, voxel->hi.z, 
         box.lo.x, box.lo.y, box.lo.z, box.hi.x, box.hi.y, box.hi.z, 
         voxel->is_overlap(box));
  */
}

void Octree::boundaries(const Boxy &box, std::vector< Octree* > &leaves) {
  if (kind == voxel_leaf_kind) {
    print_voxel_box(this, box);
    if (is_overlap(box)) {
      leaves.push_back(this);
    }
  } else if (kind == voxel_parent_kind) {
    for (int i = 0; i < 8; i++) {
      if (u.children[i]->is_overlap(box)) {
        u.children[i]->boundaries(box, leaves);
      }
    }
  }
}

void Octree::nbrs (const Boxy &box, std::vector< Octree* > &nbrs) {
  Octree* p = parent;
  if (p != NULL)
    print_voxel_box(p, box);
  while (p != NULL && !p->is_overlap(box)) {
    print_voxel_box(p, box);
    p = p->parent;
  }
  if (p != NULL) {
    print_voxel_box(p, box);
    p->boundaries(box, nbrs);
  }
}

std::vector< Octree* > Octree::all_nbrs (void) {
  std::vector< Octree* > res;
  // nbrs(box(vec(lo.x+EPSILON,lo.y+EPSILON,lo.z-EPSILON),vec(hi.x-EPSILON,hi.y-EPSILON,lo.z-EPSILON)), res);
  if (has_z_nbr())
    nbrs(box(vec(lo.x+EPSILON,lo.y+EPSILON,hi.z+EPSILON),vec(hi.x-EPSILON,hi.y-EPSILON,hi.z+EPSILON)), res);
  // nbrs(box(vec(lo.x+EPSILON,lo.y-EPSILON,lo.z+EPSILON),vec(hi.x-EPSILON,lo.y-EPSILON,hi.z-EPSILON)), res);
  if (has_y_nbr())
    nbrs(box(vec(lo.x+EPSILON,hi.y+EPSILON,lo.z+EPSILON),vec(hi.x-EPSILON,hi.y+EPSILON,hi.z-EPSILON)), res);
  // nbrs(box(vec(lo.x-EPSILON,lo.y+EPSILON,lo.z+EPSILON),vec(lo.x-EPSILON,hi.y-EPSILON,hi.z-EPSILON)), res);
  if (has_x_nbr())
    nbrs(box(vec(hi.x+EPSILON,lo.y+EPSILON,lo.z+EPSILON),vec(hi.x+EPSILON,hi.y-EPSILON,hi.z-EPSILON)), res);
  return res;
}

std::vector<Tri> Octree::boundary_voxel_triangles(flo_t r) {
  std::vector<Octree*> voxels;
  boundary_voxels(voxels);
  printf("%d BOUNDARY VOXELS R %f\n", (int)voxels.size(), r);
  return voxels_triangles(voxels);
}

std::vector<Tri> Octree::boundary_triangles(flo_t r) {
  std::vector<Octree*> voxels;
  boundary_voxels(voxels);
  printf("%d BOUNDARY VOXELS R %f\n", (int)voxels.size(), r);
  std::vector<Tri> res;
  for (auto voxel : voxels) {
    Vec off = voxel->hi - voxel->lo;
    Vec c   = voxel->center();
    Vec l1  = c - off;
    Vec h1  = c + off;
    Vec l2  = voxel->lo;
    Vec h2  = voxel->hi;
    // int lc  = lookup(c);
    auto lc = voxel;
    if (is_border(lookup(vec(l1.x, c.y, c.z)), lc)) 
      two_tris(vec(l2.x,l2.y,l2.z),vec(l2.x,l2.y,h2.z),vec(l2.x,h2.y,l2.z),vec(l2.x,h2.y,h2.z),res);
    if (is_border(lookup(vec(c.x, l1.y, c.z)), lc))
      two_tris(vec(l2.x,l2.y,l2.z),vec(l2.x,l2.y,h2.z),vec(h2.x,l2.y,l2.z),vec(h2.x,l2.y,h2.z),res);
    if (is_border(lookup(vec(c.x, c.y, l1.z)), lc))
      two_tris(vec(l2.x,l2.y,l2.z),vec(l2.x,h2.y,l2.z),vec(h2.x,l2.y,l2.z),vec(h2.x,h2.y,l2.z),res);
    if (is_border(lookup(vec(h1.x, c.y, c.z)), lc))
      two_tris(vec(h2.x,l2.y,l2.z),vec(h2.x,l2.y,h2.z),vec(h2.x,h2.y,l2.z),vec(h2.x,h2.y,h2.z),res);
    if (is_border(lookup(vec(c.x, h1.y, c.z)), lc))
      two_tris(vec(l2.x,h2.y,l2.z),vec(l2.x,h2.y,h2.z),vec(h2.x,h2.y,l2.z),vec(h2.x,h2.y,h2.z),res);
    if (is_border(lookup(vec(c.x, c.y, h1.z)), lc))
      two_tris(vec(l2.x,l2.y,h2.z),vec(l2.x,h2.y,h2.z),vec(h2.x,l2.y,h2.z),vec(h2.x,h2.y,h2.z),res);
  }
  return res;
}

std::vector< Segment<TV3> > voxels_boundary_lines(std::vector< Octree* > &voxels) {
  std::vector< Segment<TV3> > res;
  std::map< Octree*, Vec > centers;
  for (auto voxel : voxels) {
    centers[voxel] = voxel->boundary_point();
    /*
    do {
      voxel->boundary_ctr = voxel->boundary_point(voxel->ctr + (rnd_vec_of(-voxel->rad, voxel->rad) * 0.1));
      // voxel->boundary_ctr = voxel->sample_point();
    } while (false && voxel->boundary_ctr.x == INFTY);
    if (voxel->boundary_ctr.x == INFTY) {
      // printf("FAILED\n");
      voxel->boundary_ctr = voxel->ctr;
    }
    */
    // printf("BC %f,%f,%f\n", voxel->boundary_ctr.x, voxel->boundary_ctr.y, voxel->boundary_ctr.z);
  }
  for (auto voxel : voxels) {
    auto nbrs = voxel->all_nbrs();
    for (auto nbr : nbrs) {
      res.push_back(segment(centers[voxel], centers[nbr]));
    }
  }
  return res;
}

static const uint8_t VERTEX_LOOP[] = {6, 4, 5, 1, 3, 2, 6};

// Based on which vertices are filled, this map tells you which
// edges to interpolate between when forming zero, one, or two
// triangles.
// (filled vertex is first in the pair)
static const int EDGE_MAP[16][2][3][2] = {
    {{{-1,-1}, {-1,-1}, {-1,-1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // ----
    {{{ 0, 2}, { 0, 1}, { 0, 3}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // ---0
    {{{ 1, 0}, { 1, 2}, { 1, 3}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // --1-
    {{{ 1, 2}, { 1, 3}, { 0, 3}}, {{ 0, 3}, { 0, 2}, { 1, 2}}}, // --10
    {{{ 2, 0}, { 2, 3}, { 2, 1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // -2--
    {{{ 0, 3}, { 2, 3}, { 2, 1}}, {{ 2, 1}, { 0, 1}, { 0, 3}}}, // -2-0
    {{{ 1, 0}, { 2, 0}, { 2, 3}}, {{ 2, 3}, { 1, 3}, { 1, 0}}}, // -21-
    {{{ 2, 3}, { 1, 3}, { 0, 3}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // -210

    {{{ 3, 0}, { 3, 1}, { 3, 2}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 3---
    {{{ 3, 2}, { 0, 2}, { 0, 1}}, {{ 0, 1}, { 3, 1}, { 3, 2}}}, // 3--0
    {{{ 1, 2}, { 3, 2}, { 3, 0}}, {{ 3, 0}, { 1, 0}, { 1, 2}}}, // 3-1-
    {{{ 1, 2}, { 3, 2}, { 0, 2}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 3-10
    {{{ 3, 0}, { 3, 1}, { 2, 1}}, {{ 2, 1}, { 2, 0}, { 3, 0}}}, // 32--
    {{{ 3, 1}, { 2, 1}, { 0, 1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 32-0
    {{{ 3, 0}, { 1, 0}, { 2, 0}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 321-
    {{{-1,-1}, {-1,-1}, {-1,-1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 3210
};

////////////////////////////////////////////////////////////////////////////////


/*  mesh_zero_crossing
 *
 *  Finds a zero crossing for a given leaf.  May store new vertices in the
 *  vdata buffer (if we haven't already solved for this crossing).
 *
 *  Returns an index into the vdata buffer.
 */
uint32_t Octree::mesh_zero_crossing(const Octree* const neighbors[6],
                                    const uint8_t v0, const uint8_t v1,
                                    Meshy* mesh) {

  if (!data.tri) {
    data.tri = (uint32_t*)calloc(64, sizeof(uint32_t));
  }

  // If we've already solved for this zero crossing then return
  // the appropriate vertex id
  if (data.tri[(v0<<3)|v1]) {
    return data.tri[(v0<<3)|v1] - 1;
  }

  // First, look up this edge in our neighbors to see if we
  // find a match.  If we do, then save the vertex ID and set
  // found to 'true'
  bool found = false;
  uint32_t id;

  for (uint8_t axis=0; axis < 3 && !found; ++axis) {

    // If the edge doesn't vary along this axis, we can compare
    // with the appropriate neighbor.
    uint8_t mask = 1 << axis;
    if (!(mask & ~(v0 ^ v1))) continue;

    // Find the appropriate neighbor
    const Octree* neighbor = neighbors[axis*2 + ((v0 & mask) ? 1 : 0)];
    // If we don't have this neighbor or this neighbor doesn't have
    // a vertex cache, then keep going.
    if (!neighbor || !neighbor->data.tri)  continue;

    // Look up this vertex in the neighbor's edge cache
    uint16_t index = ((v0 ^ mask)<<3) | (v1 ^ mask);
    // If we don't find it, keep going
    if (!neighbor->data.tri[index]) continue;

    // We found the vertex!  Save it in the id variable.
    found = true;
    id = neighbor->data.tri[index] - 1;
    data.tri[(v0 << 3)|v1] = id + 1;
    data.tri[(v1 << 3)|v0] = id + 1;
  }

  // If we didn't find a match among neighbors, solve for the
  // zero crossing via interpolation.
  if (!found) {
    auto c = zero_crossing(v0, v1);

    id = mesh->vdata.size();
    data.tri[(v0<<3)|v1] = id+1;
    data.tri[(v1<<3)|v0] = id+1;

    mesh->vdata.append(c);
    mesh->gdata.append(vec(0.0, 0.0, 0.0));
    // printf("  C [%f,%f,%f]\n", c.x, c.y, c.z);
  }

  // The gradient will be a sum from all of the cells
  // that share this vertex.
  auto g = gradient_one(mesh->vdata[id]);
  mesh->gdata[id] = mesh->gdata[id] + g;
  // printf("ID %d [%f,%f,%f] -> [%f,%f,%f]\n", id, g.x, g.y, g.z, mesh->gdata[id].x, mesh->gdata[id].y, mesh->gdata[id].z);

  return id;
}
////////////////////////////////////////////////////////////////////////////////


/*  evaluate_voxel
 *
 *  Uses marching tetrahedrons to generate a set of triangles
 *  for a given leaf cell.  Triangle vertices in idata as indexes into
 *  the vdata list, and icount is updated appropriately.
 */
void Octree::evaluate_voxel(const Octree* const neighbors[6], Meshy* mesh) {

  // Loop over the six tetrahedra that make up a voxel cell
  for (int t = 0; t < 6; ++t) {

    // Find vertex positions for this tetrahedron
    uint8_t vertices[] = {0, 7, VERTEX_LOOP[t], VERTEX_LOOP[t+1]};

    // Figure out which of the sixteen possible combinations
    // we're currently experiencing.
    uint8_t lookup = 0;
    for (int vertex = 3; vertex >= 0; --vertex) {
      lookup = (lookup << 1) + (u.d[vertices[vertex]] < 0);
    }

    if (EDGE_MAP[lookup][0][0][0] == -1)    continue;


    // Store the first triangle in the vertex array
    IV3 tri;
    for (int v=0; v < 3; ++v) {
      tri[v] = mesh_zero_crossing(neighbors,
                                  vertices[EDGE_MAP[lookup][0][v][0]],
                                  vertices[EDGE_MAP[lookup][0][v][1]],
                                  mesh);
    }
    mesh->tdata.append(tri);

    // Move on to the second triangle, aborting if there isn't one.
    if (EDGE_MAP[lookup][1][0][0] == -1)    continue;

    // Store the second triangle in the vertex array
    for (int v=0; v < 3; ++v) {
      tri[v] = mesh_zero_crossing(neighbors,
                                  vertices[EDGE_MAP[lookup][1][v][0]],
                                  vertices[EDGE_MAP[lookup][1][v][1]],
                                  mesh);
    }
    mesh->tdata.append(tri);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Octree::do_triangulate(const Octree* const neighbors[6], Meshy* mesh, volatile int* const halt) {
  if (*halt) return;

  if (kind == voxel_leaf_kind) {
    evaluate_voxel(neighbors, mesh);
  } else if (kind == voxel_parent_kind) {
    const Octree* new_neighbors[6];

    for (int i=0; i < 8; ++i) {
      get_neighbors_3d(neighbors, new_neighbors, i);
      u.children[i]->do_triangulate(new_neighbors, mesh, halt);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

Meshy* Octree::triangulate(volatile int* const halt) {
  Meshy* mesh = new Meshy();
  const Octree* neighbors[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
  do_triangulate(neighbors, mesh, halt);
  return mesh;
}

/*
    Data tables for marching squares


    Vertices have the following IDs
   1        3
    --------
    |      |
    |      |       ^ Y
    |      |       |
    --------       --> X
   0        2

    (divided by two from d[i] indices, since we're
     looking at a flat ASDF and we don't care about z)


    Edges are numbered as follows:
        1
    --------
    |      |
  2 |      | 3     ^ Y
    |      |       |
    --------       --> X
        0


*/

// For a given set of filled corners, this array defines
// the cell edges from which we draw interior edges
static const int8_t EDGE_MAP2[16][2][2] = {
    {{-1, -1}, {-1, -1}}, // ----
    {{0, 2}, {-1, -1}},  // ---0
    {{2, 1}, {-1, -1}},  // --1-
    {{0, 1}, {-1, -1}},  // --10
    {{3, 0}, {-1, -1}},  // -2--
    {{3, 2}, {-1, -1}},  // -2-0
    {{3, 0}, { 2,  1}},  // -21-
    {{3, 1}, {-1, -1}},  // -210

    {{1, 3}, {-1, -1}},  // 3---
    {{1, 3}, { 0,  2}},  // 3--0
    {{2, 3}, {-1, -1}},  // 3-1-
    {{0, 3}, {-1, -1}},  // 3-10
    {{1, 0}, {-1, -1}},  // 32--
    {{1, 2}, {-1, -1}},  // 32-0
    {{2, 0}, {-1, -1}},  // 321-
    {{-1,-1}, {-1, -1}}  // 3210
};


// Indexed by edge number, returns vertex index
static const int8_t VERTEX_MAP2[4][2] = {
    {0, 2},
    {3, 1},
    {1, 0},
    {2, 3}
};

std::vector<Path*> Octree::contour(volatile int* const halt) {
  const Octree* neighbors[4] = {NULL, NULL, NULL, NULL};
  find_edges(neighbors, halt);

  std::vector<Path*> paths;
  write_edges(paths);
  // TODO: free_data(tree);

  return paths;
}

void Octree::write_edges (std::vector<Path*>& paths) {
  // printf("WRITE-EDGES %d [%f,%f][%f,%f]\n",
  //        kind, lo.x, lo.y, hi.x, hi.y);
  if (kind == voxel_leaf_kind) {
    for (int e=0; e < 4; ++e) {
      Path* p = data.contour[e];
      if (!p) continue;

      // Grab the path and trace it out
      // printf("BACKTRACE PATH %d\n", p->id);
      Path* start = backtrace_path(p->prev, p);
      // printf("DECIMATE PATH %lx\n", start);
      start = decimate_path(start, EPSILON);

      paths.push_back(start);
      // printf("DISCONNECT PATH %d\n", start->id);
      disconnect_path(start);
    }
  } else if (kind == voxel_parent_kind) {
    // Otherwise, recurse.
    for (int i=0; i < 8; i += 2) {
      u.children[i]->write_edges(paths);
    }
  }
}


void Octree::find_edges(const Octree* const neighbors[4],volatile int* const halt) {
  if (*halt) return;

  if (kind == voxel_leaf_kind) {

    evaluate_pixel(neighbors);

  } else if (kind == voxel_parent_kind) {

    // Evaluate the children
    for (int i=0; i < 8; i+=2) {
      const Octree* new_neighbors[4];
      get_neighbors_2d(neighbors, new_neighbors, i);
      u.children[i]->find_edges(new_neighbors, halt);
    }

    // Upgrade disconnected edges from children
    // (necessary for multi-scale path merging; trust me on this)
    data.contour = (Path**)calloc(4, sizeof(Path*));

    for (int e=0; e < 4; ++e) {
      for (int i=0; i < 8; i += 2) {

        // Only pull from children that touch this edge.
        if (e == 0 &&  (i & 2)) continue;
        if (e == 1 && !(i & 2)) continue;
        if (e == 2 &&  (i & 4)) continue;
        if (e == 3 && !(i & 4)) continue;

        // If we've already filled this edge or
        // this cell doesn't exist or doesn't have data,
        // or doesn't have data on this edge, then skip it.
        if (data.contour[e] ||
            !u.children[i] ||
            !u.children[i]->data.contour ||
            !u.children[i]->data.contour[e]) {
          continue;
        }

        // If this edge has a loose end, upgrade it
        Path* p = u.children[i]->data.contour[e];

        data.contour[e] = p;
        // printf("WRITE CONTOUR %d\n", p->id);

        // Add a pointer so that this path can disconnect
        // itself from the grid when needed
        p->ptrs.push_back(&(data.contour[e]));
      }
    }
  }
}

void Octree::evaluate_pixel(const Octree* const neighbors[4]) {
  if (!data.contour) {
    data.contour = (Path**)calloc(4, sizeof(Path*));
  }

  // Figure out which of the 16 possible cases we've encountered.
  // (ignore the z levels, just check xy)
  uint8_t lookup = 0;
  for (int vertex=3; vertex >= 0; --vertex) {
    lookup = (lookup << 1) + (u.d[vertex << 1] < 0);
  }

  if (EDGE_MAP2[lookup][0][0] != -1) {
    Path *p1 = contour_zero_crossing(neighbors, EDGE_MAP2[lookup][0][0]),
         *p2 = contour_zero_crossing(neighbors, EDGE_MAP2[lookup][0][1]);
    p1->next = p2;
    p2->prev = p1;
  }

  if (EDGE_MAP2[lookup][1][0] != -1) {
    Path *p1 = contour_zero_crossing(neighbors, EDGE_MAP2[lookup][1][0]),
         *p2 = contour_zero_crossing(neighbors, EDGE_MAP2[lookup][1][1]);
    p1->next = p2;
    p2->prev = p1;
  }
}

static int id = 0;

Path* Octree::contour_zero_crossing (const Octree* const neighbors[4], const uint8_t edge) {
  Path* t = NULL;

  // Check for edge intersections among neighbors:
  // Neighbors are in the order -y, +y, -x, +x

    // Lower edge
  if (edge == 0 && neighbors[0] && neighbors[0]->data.contour) {
    t = neighbors[0]->data.contour[1];
    // Upper edge
  } else if (edge == 1 && neighbors[1] && neighbors[1]->data.contour) {
    t = neighbors[1]->data.contour[0];
    // Left edge
  } else if (edge == 2 && neighbors[2] && neighbors[2]->data.contour) {
    t = neighbors[2]->data.contour[3];
    // Right edge
  } else if (edge == 3 && neighbors[3] && neighbors[3]->data.contour) {
    t = neighbors[3]->data.contour[2];
  }

  // If we've found one, then we're done.
  if (t) {
    data.contour[edge] = t;
    t->ptrs.push_back(&(data.contour[edge]));
    return t;
  }


  // If we didn't find a match among neighbors, solve for the
  // zero crossing via interpolation.
  const uint8_t v0 = VERTEX_MAP2[edge][0]<<1,
                v1 = VERTEX_MAP2[edge][1]<<1;
  const flo_t d0 = u.d[v0],
              d1 = u.d[v1];

  // Interpolation from d0 to d1
  const flo_t d = (d0)/(d0-d1);

  // Find interpolated coordinates
  const flo_t x0 = (v0 & 4) ? hi.x : lo.x,
              y0 = (v0 & 2) ? hi.y : lo.y;

  const flo_t x1 = (v1 & 4) ? hi.x : lo.x,
              y1 = (v1 & 2) ? hi.y : lo.y;

  auto c = vec(x0*(1-d)+x1*d, y0*(1-d)+y1*d, 0.0);

  // Create a point
  // printf("CREATING PATH %d [%f,%f]\n", id, c.x, c.y);
  t = new Path(c.x, c.y, id++);
  data.contour[edge] = t;

  // Add a back-reference to the containing Octree's edge pointer
  // (so that the path can disconnect itself when needed)
  t->ptrs.push_back(&(data.contour[edge]));

  return t;
}

//// MARCHING CUBES ON VOXELS

// These tables are used so that everything can be done in little loops that you can look at all at once
// rather than in pages and pages of unrolled code.

//a2fVertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static std::vector< Vec > a2fVertexOffset;
//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
static std::vector< Vec > a2fEdgeDirection;

static bool isInitMesher = false;

static void initMesher (void) {
  std::vector< Vec > r;
  r.push_back(vec(0.0, 0.0, 0.0));
  r.push_back(vec(1.0, 0.0, 0.0));
  r.push_back(vec(1.0, 1.0, 0.0));
  r.push_back(vec(0.0, 1.0, 0.0));
  r.push_back(vec(0.0, 0.0, 1.0));
  r.push_back(vec(1.0, 0.0, 1.0));
  r.push_back(vec(1.0, 1.0, 1.0));
  r.push_back(vec(0.0, 1.0, 1.0));
  a2fVertexOffset = r;
  r.clear();
  r.push_back(vec(1.0, 0.0, 0.0));
  r.push_back(vec(0.0, 1.0, 0.0));
  r.push_back(vec(-1.0, 0.0, 0.0));
  r.push_back(vec(0.0, -1.0, 0.0));
  r.push_back(vec(1.0, 0.0, 0.0));
  r.push_back(vec(0.0, 1.0, 0.0));
  r.push_back(vec(-1.0, 0.0, 0.0));
  r.push_back(vec(0.0, -1.0, 0.0));
  r.push_back(vec(0.0, 0.0, 1.0));
  r.push_back(vec(0.0, 0.0, 1.0));
  r.push_back(vec(0.0, 0.0, 1.0));
  r.push_back(vec(0.0, 0.0, 1.0));
  a2fEdgeDirection = r;
};

//a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
static const int a2iEdgeConnection[12][2] = {
  {0,1}, {1,2}, {2,3}, {3,0},
  {4,5}, {5,6}, {6,7}, {7,4},
  {0,4}, {1,5}, {2,6}, {3,7}
};

//get_offset finds the approximate point of intersection of the surface
// between two points with the values fValue1 and fValue2
double get_offset(double fValue1, double fValue2, double fValueDesired) {
  double fDelta = fValue2 - fValue1;

  if(fDelta == 0.0) 
    return 0.5;
  else
    return (fValueDesired - fValue1)/fDelta;
}

static Vec normalize_vector(Vec &rfVectorSource) {
  double fOldLength = rfVectorSource.magnitude();

  if (fOldLength == 0.0) {
    return rfVectorSource;
  } else {
    return rfVectorSource * (1.0 / fOldLength);
  }
}

//vGetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
Vec Octree::get_normal(Vec pos) {
  auto v = vec(g->dist(vec(pos.x-0.01, pos.y, pos.z)) - g->dist(vec(pos.x+0.01, pos.y, pos.z)),
               g->dist(vec(pos.x, pos.y-0.01, pos.z)) - g->dist(vec(pos.x, pos.y+0.01, pos.z)),
               g->dist(vec(pos.x, pos.y, pos.z-0.01)) - g->dist(vec(pos.x, pos.y, pos.z+0.01)));
  return normalize_vector(v);
}

void marching_cubes(std::vector<Octree*>& voxels, 
                    std::vector<Tri> &tris, std::vector<Tri> &colors, std::vector<Tri> &normals) {
  for (auto voxel : voxels) {
    voxel->march_cube(0.0, tris, normals, colors);
  }
}


//march_cube performs the Marching Cubes algorithm on a single cube
double Octree::march_cube
    (double target_value, std::vector< Tri >& tris, std::vector< Tri >& normals, std::vector< Tri >& colors) {
  extern int aiCubeEdgeFlags[256];
  extern int a2iTriangleConnectionTable[256][16];

  // double fX = lo.x;
  // double fY = lo.y;
  // double fZ = lo.z;
  double fScale = rad();
  int iEdgeFlags;
  Vec sColor;
  double afCubeValue[8];
  Vec asEdgeVertex[12];
  Vec asEdgeNorm[12];

  if (!isInitMesher) {
    initMesher(); isInitMesher = true;
  }
  // printf("CUBE %f,%f,%f:", fX, fY, fZ);
  //Make a local copy of the values at the cube's corners
  for(int iVertex = 0; iVertex < 8; iVertex++) {
    afCubeValue[iVertex] = g->dist(lo + a2fVertexOffset[iVertex] * fScale);
    // printf(" %f", afCubeValue[iVertex]);
    // if (iVertex == 0 && fabs(afCubeValue[iVertex] - target_value) > 2*fScale)
    //   return afCubeValue[iVertex];
  }
  // printf("\n");

  //Find which vertices are inside of the surface and which are outside
  unsigned int iFlagIndex = 0;
  for(int iVertexTest = 0; iVertexTest < 8; iVertexTest++) {
    if (afCubeValue[iVertexTest] < target_value) 
      iFlagIndex |= 1<<iVertexTest;
  }

  //Find which edges are intersected by the surface
  iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

  //If the cube is entirely inside or outside of the surface, then there will be no intersections
  if(iEdgeFlags == 0) {
    // post("HOMO AT [%f,%f,%f]\n", fX, fY, fZ);
    return afCubeValue[0];
  } else {
    // post("EDGE AT [%f,%f,%f]\n", fX, fY, fZ);
  }

  //Find the point of intersection of the surface with each edge
  //Then find the normal to the surface at those points
  for(int iEdge = 0; iEdge < 12; iEdge++) {
    //if there is an intersection on this edge
    if(iEdgeFlags & (1<<iEdge)) {
      auto fOffset = get_offset(afCubeValue[ a2iEdgeConnection[iEdge][0] ], 
                                afCubeValue[ a2iEdgeConnection[iEdge][1] ], target_value);

      asEdgeVertex[iEdge] = lo + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ] + a2fEdgeDirection[iEdge] * fOffset) * fScale;
      asEdgeNorm[iEdge] = get_normal(asEdgeVertex[iEdge]);
    }
  }

  //Draw the triangles that were found.  There can be up to five per cube
  for(int iTriangle = 0; iTriangle < 5; iTriangle++) {
    if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0)
      break;

    Vec vtri[3];
    Vec vcolor[3];
    Vec vnormal[3];
    for(int iCorner = 0; iCorner < 3; iCorner++) {
      auto iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];

      auto sColor = get_color(asEdgeNorm[iVertex]);
      // glColor3f(sColor.x, sColor.y, sColor.z);
      // glNormal3f(asEdgeNorm[iVertex].x,   asEdgeNorm[iVertex].y,   asEdgeNorm[iVertex].z);
      vtri[iCorner] = asEdgeVertex[iVertex];
      vcolor[iCorner] = sColor;
      vnormal[iCorner] = asEdgeNorm[iVertex];
      // tri[iCorner] = asEdgeVertex[iVertex];
      // post("EDGE [%f,%f,%f]\n", asEdgeVertex[iVertex].x, asEdgeVertex[iVertex].y, asEdgeVertex[iVertex].z);
      // post("EDGE [%f,%f,%f]\n", tri[iCorner].x, tri[iCorner].y, tri[iCorner].z);
    }
    tris.push_back(tri(vtri[0],vtri[1],vtri[2] /*,vnormal[0] */));
    colors.push_back(tri(vcolor[0],vcolor[1],vcolor[2]));
    normals.push_back(tri(vnormal[0],vnormal[1],vnormal[2]));
    // printf("TRI %f,%f,%f\n", vtri[0].x, vtri[0].y, vtri[0].z);
  }
  return afCubeValue[0];
}

// For any edge, if one vertex is inside of the surface and the other is outside of the surface
//  then the edge intersects the surface
// For each of the 8 vertices of the cube can be two possible states : either inside or outside of the surface
// For any cube the are 2^8=256 possible sets of vertex states
// This table lists the edges intersected by the surface for all 256 possible vertex states
// There are 12 edges.  For each entry in the table, if edge #n is intersected, then bit #n is set to 1

int aiCubeEdgeFlags[256]= {
  0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 
  0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 
  0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 
  0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 
  0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0, 
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460, 
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0, 
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230, 
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190, 
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

//  For each of the possible vertex states listed in aiCubeEdgeFlags there is a specific triangulation
//  of the edge intersection points.  a2iTriangleConnectionTable lists all of them in the form of
//  0-5 edge triples with the list terminated by the invalid value -1.
//  For example: a2iTriangleConnectionTable[3] list the 2 triangles formed when corner[0] 
//  and corner[1] are inside of the surface, but the rest of the cube is not.
//
//  I found this table in an example program someone wrote long ago.  It was probably generated by hand

int a2iTriangleConnectionTable[256][16] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
  {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
  {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
  {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
  {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
  {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
  {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
  {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
  {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
  {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
  {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
  {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
  {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
  {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
  {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
  {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
  {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
  {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
  {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
  {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
  {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
  {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
  {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
  {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
  {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
  {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
  {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
  {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
  {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
  {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
  {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
  {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
  {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
  {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
  {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
  {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
  {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
  {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
  {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
  {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
  {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
  {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
  {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
  {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
  {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
  {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
  {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
  {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
  {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
  {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
  {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
  {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
  {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
  {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
  {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
  {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
  {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
  {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
  {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
  {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
  {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
  {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
  {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
  {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
  {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
  {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
  {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
  {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
  {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
  {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
  {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
  {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
  {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
  {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
  {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
  {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
  {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
  {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
  {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
  {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
  {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
  {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
  {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
  {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
  {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
  {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
  {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
  {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
  {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
  {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
  {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
  {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
  {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
  {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
  {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
  {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
  {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

Octree* g_octree_val(Geom* g) {
  ensure(g->k == octree_kind, "EXPECTING OCTREE KIND");
  return ((Octree*)g);
}

Octree* expr_to_tree(Expr* e, T radius, T threshold) {
  auto s  = e->to_str();
  auto mg = e->map(expr_x(),expr_y(),expr_z());
  auto fg = mg->fold();
  fg->build();
  // printf("GEOM = %s\n", s.c_str());
  auto tree = octree(fg, vec(-radius, -radius, -radius), vec(radius, radius, radius), radius, threshold);
  // tree = octree(NULL, e, vec(0, 0, 0), radius, threshold, 0);
  // if (is_timing)
  //   exit(-1);
  // printf("TREE %d\n", tree->k);
  return tree;
}

Mesh meshy_to_mesh (Meshy* mesh) {
  return fab_mesh(mesh->tdata, mesh->vdata);
}

Mesh octree_to_mesh(Octree* g) {
  int halt = 0;
  auto tree = (Octree*)g;
  auto mesh = tree->triangulate(&halt);
  return meshy_to_mesh(mesh);
}

extern flo_t get_radius(void), get_threshold(void);

Octree* pretty_print_tree(Octree* g) {
  auto tree = (Octree*)g;
  auto radius = get_radius();
  auto threshold = get_threshold();
  printf("TREE NODE           = %luB\n", sizeof(Octree));
  printf("TREE NODES          = %ld\n", tree->count_nodes());
  printf("MAX VOXELS          = %.2fM\n", pow(2 * radius / threshold,3) / (double)1e6);
  printf("TREE COMPRESSION    = %.2f\n", pow(2 * radius / threshold,3) / tree->count_nodes());
  printf("TREE SIZE           = %.2fMB\n", tree->count_nodes() * sizeof(Octree) / 1e6);
  printf("ROOT EXPR NODES     = %ld\n", tree->g->count());
  printf("MAX TREE EXPR NODES = %.2fM\n", (tree->count_nodes() * tree->g->count() / 1e6));
  printf("TREE EXPR NODES     = %.2fM\n", tree->count_fun() / 1e6);
  printf("EXPR COMPRESSION    = %.2f\n", (tree->count_nodes() * tree->g->count()) / (double)tree->count_fun());
  return tree;
}

Nested< TV2 > tree_slice (Octree* tree, T z) {
  // TODO: Z is ignored
  auto slice = tree->slice(0.0);
  int halt = 0;
  auto paths = slice->contour(&halt);
  Nested< TV2, false > contours;
  for (auto start : paths) {
    Array< TV2 > points;
    auto p = start;
    do {
      points.append(vec(p->x, p->y));
      p = p->next;
    } while (p && p != start);
    contours.append(points);
  }
  contours.freeze();
  return contours;
}
