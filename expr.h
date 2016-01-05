#ifndef __IS_EXPR__
#define __IS_EXPR__

#include "cad.h"
#include "geom-interface.h"
#include <geode/geometry/SimplexTree.h>

enum ExprKind {
  expr_lit_kind,
  expr_x_kind,
  expr_y_kind,
  expr_z_kind,
  expr_space_kind,
  expr_neg_kind,
  expr_abs_kind,
  expr_sqrt_kind,
  expr_sin_kind,
  expr_cos_kind,
  expr_tan_kind,
  expr_asin_kind,
  expr_acos_kind,
  expr_atan_kind,
  expr_min_kind,
  expr_max_kind,
  expr_add_kind,
  expr_sub_kind,
  expr_mul_kind,
  expr_pow_kind,
  expr_div_kind,
  expr_xform_kind,
  expr_mesh_kind
};

class Expr;

extern Expr* expr_lit(flo_t a);
extern Expr* expr_all(void);
extern Expr* expr_none(void);
extern Expr* expr_pi(void);
extern Expr* expr_x(void);
extern Expr* expr_y(void);
extern Expr* expr_z(void);
extern Expr* expr_neg(Expr* g);
extern Expr* expr_not(Expr* g);
extern Expr* expr_min (Expr* a, Expr* b);
extern Expr* expr_or (Expr* a, Expr* b);
extern Expr* expr_max (Expr* a, Expr* b);
extern Expr* expr_and (Expr* a, Expr* b);
extern Expr* expr_add (Expr* a, Expr* b);
extern Expr* expr_sub (Expr* a, Expr* b);
extern Expr* expr_mul (Expr* a, Expr* b);
extern Expr* expr_sqr (Expr* a);
extern Expr* expr_div (Expr* a, Expr* b);
extern Expr* expr_abs (Expr* a);
extern Expr* expr_sqrt (Expr* a);
extern Expr* expr_sin (Expr* g);
extern Expr* expr_cos (Expr* g);
extern Expr* expr_tan (Expr* a);
extern Expr* expr_pow (Expr* a, Expr* b);
extern Expr* expr_asin (Expr* a);
extern Expr* expr_acos (Expr* a);
extern Expr* expr_atan (Expr* a);
extern Expr* expr_xform(Expr* nx, Expr* ny, Expr* nz, Expr* g);
extern Expr* expr_mov(Expr* dx, Expr* dy, Expr* dz, Expr* g);
extern Expr* expr_xmov(Expr* dx, Expr* g);
extern Expr* expr_ymov(Expr* dy, Expr* g);
extern Expr* expr_zmov(Expr* dz, Expr* g);
extern Expr* expr_mag(Expr* dx, Expr* dy, Expr* dz, Expr* g);
extern Expr* expr_mag1(Expr* dxyz, Expr* g);
extern Expr* expr_xmag(Expr* dx, Expr* g);
extern Expr* expr_ymag(Expr* dy, Expr* g);
extern Expr* expr_zmag(Expr* dz, Expr* g);
extern Expr* expr_yrevolve(Expr* g);
extern Expr* expr_xrevolve(Expr* g);
extern Expr* expr_shear_x_z(Expr* z0, Expr* z1, Expr* dx0, Expr* dx1, Expr* g);
extern Expr* expr_taper(Expr* dz, Expr* s0, Expr* s1, Expr* g);
extern Expr* expr_xrot(Expr* a, Expr* g);
extern Expr* expr_reflect_x(Expr* g);
extern Expr* expr_reflect_xy(Expr* g);
extern Expr* expr_yrot(Expr* a, Expr* g);
extern Expr* expr_reflect_y(Expr* g);
extern Expr* expr_zrot(Expr* a, Expr* g);
extern Expr* expr_reflect_z(Expr* g);
extern Expr* expr_reflect_xz(Expr* g);
extern Expr* expr_reflect_yz(Expr* g);
extern Expr* expr_half(Expr* nx, Expr* ny, Expr* nz, Expr* d);
extern Expr* expr_edge(Expr* x0, Expr* y0, Expr* x1, Expr* y1);
extern Expr* expr_triangle(Expr* x0, Expr* y0, Expr* x1, Expr* y1, Expr* x2, Expr* y2);
extern Expr* expr_square(Expr* d);
extern Expr* expr_rem(Expr* a, Expr* b);
extern Expr* expr_hollow(Expr* a, Expr* r);
extern Expr* expr_circle(Expr* d);
extern Expr* expr_sphere(Expr* d);
extern Expr* expr_extrude(Expr* d, Expr* g);
extern Expr* expr_xbox(Expr* a, Expr* b);
extern Expr* expr_align_xmin(Expr* a, Expr* b);
extern Expr* expr_align_xmax(Expr* a, Expr* b);
extern Expr* expr_align_ymin(Expr* a, Expr* b);
extern Expr* expr_align_ymax(Expr* a, Expr* b);
extern Expr* expr_align_zmin(Expr* a, Expr* b);
extern Expr* expr_align_zmax(Expr* a, Expr* b);
extern Expr* expr_ybox(Expr* a, Expr* b);
extern Expr* expr_zbox(Expr* a, Expr* b);
extern Expr* expr_cylinder(Expr* xyrad, Expr* zrad);
extern Expr* expr_capsule(Expr* xyrad, Expr* zrad);
// extern Meshy* tree_save(Meshy* g, std::string s);
extern Expr* expr_cube(Expr* r);
extern Expr* expr_space(Expr* r);
extern Expr* expr_rect2(Expr* dx, Expr* dy);
extern Expr* expr_rect3(Expr* dx, Expr* dy, Expr* dz);
extern Expr* expr_cone(Expr* xyd, Expr* zd);
extern Expr* expr_pyramid(Expr* xyd, Expr* zd);
extern Expr* expr_blend(Expr* g1, Expr* g2, Expr* amount);
extern Expr* pretty_print_expr(Expr* e);
extern Expr* polygonize (Nested<TV2> p);
extern Expr* expr_mesh(Mesh mesh);

extern Geom* parse_geom(std::string s);
extern Expr* parse_expr(std::string s);

class Fun {
 public:
  virtual flo_t dist (const Vec &p) { return 0.0; }
  virtual Fun* maybe_prune (int d, Interval &i, const Boxy &bounds) {
    printf("DEFAULT MAYBE_PRUNE\n");
    return this; }
  virtual long size (void) { return 0; }
  virtual long count (void) { return 0; }
  virtual std::string to_str (void) { return ""; }
};

class Expr : public Fun, public Geom {
 public:
  ExprKind ek; // kind
  bool     is_dead;
  flo_t    d; // dist
  Interval i; // interval
  int      c; // count
  int      o; // offset
  Expr *xmin, *xmax;
  Expr *ymin, *ymax;
  Expr *zmin, *zmax;
  // Boxy      bbox;
  // bool     is_bbox;
  std::vector< Expr* > elts;
  void copy_bbox(Expr* g) { 
    xmin = g->xmin; xmax = g->xmax; 
    ymin = g->ymin; ymax = g->ymax; 
    zmin = g->zmin; zmax = g->zmax; 
  }
  virtual long size(void) { return elts.size(); }
  virtual long count(void) { return c = 0; }
  virtual bool is_same(Expr* o) { return o->ek == ek && o->c == c; }
  Interval range (const Boxy &bounds) {
    for (auto g : elts)
      g->do_range(bounds);
    return i;
  }
  bool is_prune (void) {
    for (auto g : elts)
      if (g->do_is_prune())
        return true;
    return false;
  }
  virtual flo_t dist (const Vec &p) { 
    // printf("DIST %d ELTS\n", (int)elts.size());
    for (auto g : elts)
      g->do_dist(p);
    return d;
  }
  virtual Fun* maybe_prune (int d, Interval &i, const Boxy &bounds) {
    i = range(bounds);
    auto ng = is_prune() ? prune(d) : this;
    return ng;
  }
  Expr* prune (int d) {
    // for (int i = 0; i < (d+1); i++) printf("  ");
    // int before = elts.size();
    auto pg = do_prune();
    pg->build();
    // int after = pg->elts.size();
    // printf("PRUNING... %d -> %d - %.2f%%\n",
    //        before, after, 100 * (before - after)/(double)before);
    return pg;
  }
  void build (void) {
    count();
    // printf("GEOM %s\n", to_str().c_str());
    std::map< int,std::vector<Expr*> > seen;
    cse(seen);
    std::set< Expr* > is_visited;
    elts.clear(); sched(elts, is_visited);
    for (size_t i = 0; i < elts.size(); i++) 
      elts[i]->o = i;
    // for (size_t i = 0; i < elts.size(); i++) 
    //   printf("%2zu: %s\n", i, elts[i]->sched_str().c_str());
    // printf("BUILD %d ELTS\n", (int)elts.size());
  }
  Expr* lookup(std::vector<Expr*> &seen) {
    for (auto g : seen) {
      if (is_same(g)) {
        // printf("SAME %s == %s\n", this->to_str().c_str(), g->to_str().c_str());
        return g;
      }
    }
    return NULL;
  }
  Expr* find_update(std::map< int, std::vector<Expr*> > &seen) { 
    // printf("FIND UPDATE K %d C %d S %s\n", k, c, to_str().c_str());
    auto elts = seen.find(k);
    if (elts == seen.end()) {
      std::vector< Expr* > res;
      res.push_back(this);
      seen[k] = res;
      // printf("  K NOT FOUND -- ADDING\n");
      return this;
    } else {
      auto olds = seen[k];
      auto old  = lookup(olds);
      // printf("  K FOUND -- ");
      if (old == NULL) {
        // printf("GEOM NOT FOUND -- ADDING\n");
        olds.push_back(this);
        seen[k] = olds;
        return this;
      } else {
        // printf("GEOM FOUND\n");
        return old;
      }
    }
  }
  virtual void cse(std::map< int, std::vector<Expr*> > &seen) { }
  virtual Expr* get_a(void) { return NULL; }
  virtual Expr* get_b(void) { return NULL; }
  virtual Expr* get_q(void) { return NULL; }
  virtual void sched (std::vector< Expr* > &elts, std::set< Expr* > &is_visited) { }
  virtual void do_dist (const Vec &p) { d = 0.0; }
  virtual void do_range (const Boxy &r) { i = interval(0.0, 0.0); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return this; }
  virtual Expr* fold (void) { return this; }
  virtual Expr* do_prune (void) { return this; }
  virtual bool do_is_prune (void) { return false; }
  virtual std::string op_str (void) { return ""; }
  virtual std::string sched_str (void) { return ""; }
  virtual bool is_lit ( void ) { return false; }
  virtual flo_t lit_val ( void ) { return 0.0; }
  Expr (ExprKind k) : Geom(expr_kind),ek(k),is_dead(true),xmin(NULL),xmax(NULL),ymin(NULL),ymax(NULL),zmin(NULL),zmax(NULL)
    // xmin(neg_inf()),xmax(pos_inf()),ymin(neg_inf()),ymax(pos_inf()),zmin(neg_inf()),zmax(pos_inf())
    { }
};

class LeafExpr : public Expr {
 public:
  virtual long count (void) { return c = 1; }
  virtual void sched (std::vector< Expr* > &elts, std::set< Expr* > &is_visited) { 
    if (is_visited.find(this) == is_visited.end()) {
      is_visited.insert(this);
      elts.push_back(this); 
    }
  }
  virtual std::string op_str (void) { return to_str(); }
  virtual std::string sched_str (void) { return to_str(); }
 LeafExpr(ExprKind k) : Expr(k) { }
};

class ExprLit : public LeafExpr {
 public:
  flo_t val;
  virtual void do_dist (const Vec &p) { d = val; }
  virtual void do_range (const Boxy &r) { i = interval(val, val); }
  virtual bool is_same(Expr* o) { return o->is_lit() && o->lit_val() == val; }
  virtual std::string to_str (void) { return std::to_string(val); }
  virtual bool is_lit ( void ) { return true; }
  virtual flo_t lit_val ( void ) { return val; }
 ExprLit(flo_t val) : LeafExpr(expr_lit_kind), val(val) { }
};

class ExprX : public LeafExpr {
 public:
  virtual void do_dist (const Vec &p) { d = p.x; }
  virtual void do_range (const Boxy &r) { i = interval(r.lo.x, r.hi.x); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return nx; }
  virtual std::string to_str (void) { return "X()"; }
 ExprX(void) : LeafExpr(expr_x_kind) { }
};

class ExprY : public LeafExpr {
 public:
  virtual void do_dist (const Vec &p) { d = p.y; }
  virtual void do_range (const Boxy &r) { i = interval(r.lo.y, r.hi.y); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return ny; }
  virtual std::string to_str (void) { return "Y()"; }
 ExprY(void) : LeafExpr(expr_y_kind) { }
};

class ExprZ : public LeafExpr {
 public:
  virtual void do_dist (const Vec &p) { d = p.z; }
  virtual void do_range (const Boxy &r) { i = interval(r.lo.z, r.hi.z); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return nz; }
  virtual std::string to_str (void) { return "Z()"; }
 ExprZ(void) : LeafExpr(expr_z_kind) { }
};

static const T sqrt2 = sqrt(2.0);

class ExprSpace : public LeafExpr {
 public:
  virtual void do_dist (const Vec &p) { d = INFTY; }
  virtual void do_range (const Boxy &r) { i = interval(INFTY, INFTY); }
  virtual std::string to_str (void) { return "Space()"; }
 ExprSpace(void) : LeafExpr(expr_space_kind) { }
};

class UnaryExpr : public Expr {
 public:
  Expr* a;
  virtual Expr* get_a(void) { return a; }
  virtual void sched (std::vector< Expr* > &elts, std::set< Expr* > &is_visited) { 
    if (is_visited.find(this) == is_visited.end()) {
      is_visited.insert(this);
      a->sched(elts,is_visited); elts.push_back(this); 
    }
  }
  virtual bool is_same(Expr* o) { 
    return o->ek == ek && o->c == c && a->is_same(o->get_a());
  }
  virtual void cse(std::map< int, std::vector<Expr*> > &seen) { 
    // printf("CSE %s\n", to_str().c_str());
    auto na = a->find_update(seen);
    if (na == a) 
      a->cse(seen);
    else
      a = na;
  }
  virtual std::string to_str (void) { return op_str() + "(" + a->to_str() + ")"; }
  virtual std::string sched_str (void) { return op_str() + "(" + std::to_string(a->o) + ")"; }
  virtual long count (void) { return c = a->count() + 1; }
 UnaryExpr(ExprKind k, Expr* a) : Expr(k), a(a) { }
};

class ExprNeg : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = -(a->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval(-a->i.hi, -a->i.lo);
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(-fa->lit_val());
    else
      return expr_neg(fa);
  }
  virtual Expr* do_prune (void) { return expr_neg(a->do_prune()); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_neg(a->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "neg"; }
 ExprNeg(Expr* a) : UnaryExpr(expr_neg_kind, a) { }
};

class ExprAbs : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = fabs(a->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval(a->i.lo < 0.0 ? 0.0 : fmin(fabs(a->i.lo), fabs(a->i.hi)),
                 fmax(fabs(a->i.lo), fabs(a->i.hi)));
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_abs(a->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(fabs(fa->lit_val()));
    else
      return expr_abs(fa);
  }
  virtual Expr* do_prune (void) { return expr_abs(a->do_prune()); }
  virtual std::string op_str (void) { return "abs"; }
 ExprAbs(Expr* a) : UnaryExpr(expr_abs_kind, a) { }
};

class ExprSqrt : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = sqrt(a->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval(a->i.lo <= 0.0 ? 0.0 : sqrt(a->i.lo),
                 a->i.hi <= 0.0 ? 0.0 : sqrt(a->i.hi));
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_sqrt(a->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(sqrt(fa->lit_val()));
    else
      return expr_sqrt(fa);
  }
  virtual Expr* do_prune (void) { return expr_sqrt(a->do_prune()); }
  virtual std::string op_str (void) { return "sqrt"; }
 ExprSqrt(Expr* a) : UnaryExpr(expr_sqrt_kind, a) { }
};

class ExprSin : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = sin(a->d); }
  virtual void do_range (const Boxy &r) { 
    if (a->i.lo == a->i.hi) {
      i = interval(sin(a->i.lo), a->i.lo);
    } else {
      i = interval(-1.0, 1.0);
    }
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_sin(a->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(sin(fa->lit_val()));
    else
      return expr_sin(fa);
  }
  virtual Expr* do_prune (void) { return expr_sin(a->do_prune()); }
  virtual std::string op_str (void) { return "sin"; }
 ExprSin(Expr* a) : UnaryExpr(expr_sin_kind, a) { }
};


class ExprCos : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = cos(a->d); }
  virtual void do_range (const Boxy &r) { 
    if (a->i.lo == a->i.hi) {
      i = interval(cos(a->i.lo), a->i.lo);
    } else {
      i = interval(-1.0, 1.0);
    }
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(cos(fa->lit_val()));
    else
      return expr_cos(fa);
  }
  virtual Expr* do_prune (void) { return expr_cos(a->do_prune()); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_cos(a->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "cos"; }
 ExprCos(Expr* a) : UnaryExpr(expr_cos_kind, a) { }
};


class ExprTan : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = tan(a->d); }
  virtual void do_range (const Boxy &r) { 
    if (a->i.lo == a->i.hi) {
      i = interval(tan(a->i.lo), a->i.lo);
    } else {
      i = interval(-INFTY, INFTY);
    }
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_tan(a->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(tan(fa->lit_val()));
    else
      return expr_tan(fa);
  }
  virtual Expr* do_prune (void) { return expr_tan(a->do_prune()); }
  virtual std::string op_str (void) { return "tan"; }
 ExprTan(Expr* a) : UnaryExpr(expr_tan_kind, a) { }
};

class ExprAsin : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = asin(a->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval((a->i.lo <= -1) ? -M_PI/2 : ((a->i.lo >= 1) ? M_PI/2 : asin(a->i.lo)),
                 (a->i.hi <= -1) ? -M_PI/2 : ((a->i.hi >= 1) ? M_PI/2 : asin(a->i.hi)));
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_asin(a->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(asin(fa->lit_val()));
    else
      return expr_asin(fa);
  }
  virtual Expr* do_prune (void) { return expr_asin(a->do_prune()); }
  virtual std::string op_str (void) { return "asin"; }
 ExprAsin(Expr* a) : UnaryExpr(expr_asin_kind, a) { }
};


class ExprAcos : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = acos(a->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval((a->i.lo <= -1) ? M_PI : ((a->i.lo >= 1) ? 0 : acos(a->i.lo)),
                 (a->i.hi <= -1) ? M_PI : ((a->i.hi >= 1) ? 0 : acos(a->i.hi)));
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(acos(fa->lit_val()));
    else
      return expr_acos(fa);
  }
  virtual Expr* do_prune (void) { return expr_acos(a->do_prune()); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_acos(a->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "acos"; }
 ExprAcos(Expr* a) : UnaryExpr(expr_acos_kind, a) { }
};


class ExprAtan : public UnaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = atan(a->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval(atan(a->i.lo), atan(a->i.hi));
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_atan(a->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    if (fa->is_lit())
      return expr_lit(atan(fa->lit_val()));
    else
      return expr_atan(fa);
  }
  virtual Expr* do_prune (void) { return expr_atan(a->do_prune()); }
  virtual std::string op_str (void) { return "atan"; }
 ExprAtan(Expr* a) : UnaryExpr(expr_atan_kind, a) { }
};

class BinaryExpr : public Expr {
 public:
  Expr* a;
  Expr* b;
  virtual Expr* get_a(void) { return a; }
  virtual Expr* get_b(void) { return b; }
  virtual void sched (std::vector< Expr* > &elts, std::set< Expr* > &is_visited) { 
    if (is_visited.find(this) == is_visited.end()) {
      is_visited.insert(this);
      a->sched(elts,is_visited); b->sched(elts,is_visited); elts.push_back(this); 
    }
  }
  virtual long count (void) { return c = a->count() + b->count() + 1; }
  virtual bool is_same(Expr* o) { 
    return o->ek == ek && o->c == c && a->is_same(o->get_a()) && b->is_same(o->get_b());
  }
  virtual void cse(std::map< int, std::vector<Expr*> > &seen) { 
    // printf("CSE %s\n", to_str().c_str());
    auto na = a->find_update(seen);
    if (na == a)
      a->cse(seen);
    else
      a = na;
    auto nb = b->find_update(seen);
    if (nb == b)
      b->cse(seen);
    else
      b = nb;
  }
  virtual std::string to_str (void) { return op_str() + "(" + a->to_str() + "," + b->to_str() + ")"; }
  virtual std::string sched_str (void) { return op_str() + "(" + std::to_string(a->o) + "," + std::to_string(b->o) + ")"; }
 BinaryExpr(ExprKind k, Expr* a, Expr* b) : Expr(k), a(a), b(b) { }
};

class ExprMin : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = fmin(a->d, b->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval(fmin(a->i.lo, b->i.lo), fmin(a->i.hi, b->i.hi));
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(fmin(fa->lit_val(), fb->lit_val()));
    else if (fa->is_lit() && fa->lit_val() == INFTY)
      return fb;
    else if (fb->is_lit() && fb->lit_val() == INFTY)
      return fa;
    else
      return expr_min(fa, fb);
  }
  virtual bool do_is_prune (void) { 
    return (a->i >= b->i || b->i >= a->i);
  }
  virtual Expr* do_prune (void) { 
    if (a->i >= b->i)
      return b->do_prune();
    else if (b->i >= a->i)
      return a->do_prune();
    else
      return expr_min(a->do_prune(), b->do_prune()); 
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_min(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "min"; }
 ExprMin(Expr* a, Expr *b) : BinaryExpr(expr_min_kind, a, b) { }
};

class ExprMax : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = fmax(a->d, b->d); }
  virtual void do_range (const Boxy &r) { 
    i = interval(fmax(a->i.lo, b->i.lo), fmax(a->i.hi, b->i.hi));
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(fmax(fa->lit_val(), fb->lit_val()));
    else if (fa->is_lit() && fa->lit_val() == -INFTY)
      return fb;
    else if (fb->is_lit() && fb->lit_val() == -INFTY)
      return fa;
    else
      return expr_max(fa, fb);
  }
  virtual Expr* do_prune (void) { 
    if (a->i <= b->i)
      return b->do_prune();
    else if (b->i <= a->i)
      return a->do_prune();
    else
      return expr_max(a->do_prune(), b->do_prune()); 
  }
  virtual bool do_is_prune (void) { 
    return (a->i <= b->i || b->i <= a->i);
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) { return expr_max(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "max"; }
 ExprMax(Expr* a, Expr *b) : BinaryExpr(expr_max_kind, a, b) { }
};

class ExprAdd : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = a->d + b->d; }
  virtual void do_range (const Boxy &r) { 
    i = interval(a->i.lo + b->i.lo, a->i.hi + b->i.hi);
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(fa->lit_val() + fb->lit_val());
    else if (fa->is_lit() && fa->lit_val() == 0.0)
      return fb;
    else if (fb->is_lit() && fb->lit_val() == 0.0)
      return fa;
    else
      return expr_add(fa, fb);
  }
  virtual Expr* do_prune (void) { return expr_add(a->do_prune(), b->do_prune()); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) {
    return expr_add(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "add"; }
 ExprAdd(Expr* a, Expr* b) : BinaryExpr(expr_add_kind, a, b) { }
};

class ExprSub : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = a->d - b->d; }
  virtual void do_range (const Boxy &r) { 
    i = interval(a->i.lo - b->i.hi, a->i.hi - b->i.lo);
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) {
    return expr_sub(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(fa->lit_val() - fb->lit_val());
    else if (fa->is_lit() && fa->lit_val() == 0.0)
      return expr_neg(fb);
    else if (fb->is_lit() && fb->lit_val() == 0.0)
      return fa;
    else
      return expr_sub(fa, fb);
  }
  virtual Expr* do_prune (void) { return expr_sub(a->do_prune(), b->do_prune()); }
  virtual std::string op_str (void) { return "sub"; }
 ExprSub(Expr* a, Expr* b) : BinaryExpr(expr_sub_kind, a, b) { }
};

class ExprMul : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = a->d * b->d; }
  virtual void do_range (const Boxy &r) { 
    flo_t ll = a->i.lo * b->i.lo;
    flo_t lh = a->i.lo * b->i.hi;
    flo_t hl = a->i.hi * b->i.lo;
    flo_t hh = a->i.hi * b->i.hi;
    i = interval(fmin(fmin(ll,lh),fmin(hl,hh)),fmax(fmax(ll,lh),fmax(hl,hh)));
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(fa->lit_val() * fb->lit_val());
    else if ((fa->is_lit() && fa->lit_val() == 0.0) || (fb->is_lit() && fb->lit_val() == 0.0))
      return expr_lit(0.0);
    else if (fb->is_lit() && fb->lit_val() == 1.0)
      return fa;
    else if (fa->is_lit() && fa->lit_val() == 1.0)
      return fb;
    else
      return expr_mul(fa, fb);
  }
  virtual Expr* do_prune (void) { return expr_mul(a->do_prune(), b->do_prune()); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) {
    return expr_mul(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "mul"; }
 ExprMul(Expr* a, Expr* b) : BinaryExpr(expr_mul_kind, a, b) { }
};

class ExprDiv : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = a->d / b->d; }
  virtual void do_range (const Boxy &r) { 
    if (b->i.lo <= 0.0 && b->i.hi >= 0.0) {
      i = interval(-INFTY, INFTY);
      return;
    }
    auto cr = interval(1.0/b->i.lo, 1.0/b->i.hi);
    flo_t ll = a->i.lo * cr.lo;
    flo_t lh = a->i.lo * cr.hi;
    flo_t hl = a->i.hi * cr.lo;
    flo_t hh = a->i.hi * cr.hi;
    i = interval(fmin(fmin(ll,lh),fmin(hl,hh)),fmax(fmax(ll,lh),fmax(hl,hh)));
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) {
    return expr_div(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(fa->lit_val() / fb->lit_val());
    else if (fa->is_lit() && fa->lit_val() == 0.0)
      return expr_lit(0.0);
    else if (fb->is_lit() && fb->lit_val() == 1.0)
      return fa;
    else
      return expr_div(fa, fb);
  }
  virtual Expr* do_prune (void) { return expr_div(a->do_prune(), b->do_prune()); }
  virtual std::string op_str (void) { return "div"; }
 ExprDiv(Expr* a, Expr* b) : BinaryExpr(expr_div_kind, a, b) { }
};

class ExprPow : public BinaryExpr {
 public:
  virtual void do_dist (const Vec &p) { d = pow(a->d, b->d); }
  virtual void do_range (const Boxy &r) { 
    int p = b->i.lo;
    if (p % 2) {
      i = interval(pow(a->i.lo, p), pow(a->i.hi, p));
    } else {
      flo_t al = fabs(a->i.lo);
      flo_t ah = fabs(a->i.hi);
      i = interval((a->i.lo <= 0 && a->i.hi >= 0) ? 0 : pow(al < ah ? al : ah, p),
                   pow((al < ah) ? ah : al, p));
    }
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    if (fa->is_lit() && fb->is_lit())
      return expr_lit(pow(fa->lit_val(), fb->lit_val()));
    else if (fb->is_lit()) {
      if (fa->lit_val() == 0.0)
        return expr_lit(1.0);
      else if (fa->lit_val() == 1.0) 
        return fa;
      else
        return expr_pow(fa, fb);
    } else
      return expr_pow(fa, fb);
  }
  virtual Expr* do_prune (void) { return expr_pow(a->do_prune(), b->do_prune()); }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) {
    return expr_pow(a->map(nx, ny, nz), b->map(nx, ny, nz)); }
  virtual std::string op_str (void) { return "pow"; }
 ExprPow(Expr* a, Expr* b) : BinaryExpr(expr_pow_kind, a, b) { }
};

class TernaryExpr : public Expr {
 public:
  Expr* a;
  Expr* b;
  Expr* q;
  virtual Expr* get_a(void) { return a; }
  virtual Expr* get_b(void) { return b; }
  virtual Expr* get_q(void) { return q; }
  virtual void sched (std::vector< Expr* > &elts, std::set< Expr* > &is_visited) { 
    if (is_visited.find(this) == is_visited.end()) {
      is_visited.insert(this);
      a->sched(elts,is_visited); b->sched(elts,is_visited); q->sched(elts,is_visited); elts.push_back(this); 
    }
  }
  virtual long count (void) { return c = a->count() + b->count() + q->count() + 1; }
  virtual bool is_same(Expr* o) { 
    return o->ek == ek && o->c == c && a->is_same(o->get_a()) && b->is_same(o->get_b()) && q->is_same(o->get_q());
  }
  virtual void cse(std::map< int, std::vector<Expr*> > &seen) { 
    // printf("CSE %s\n", to_str().c_str());
    auto na = a->find_update(seen);
    if (na == a)
      a->cse(seen);
    else
      a = na;
    auto nb = b->find_update(seen);
    if (nb == b)
      b->cse(seen);
    else
      b = nb;
    auto nq = q->find_update(seen);
    if (nq == q)
      q->cse(seen);
    else
      q = nq;
  }
  virtual std::string to_str (void) { return op_str() + "(" + a->to_str() + "," + b->to_str() + "," + q->to_str() + ")"; }
  virtual std::string sched_str (void) {
    return op_str() + "(" + std::to_string(a->o) + "," + std::to_string(b->o) + "," + std::to_string(q->o) + ")"; }
 TernaryExpr(ExprKind k, Expr* a, Expr* b, Expr* q) : Expr(k), a(a), b(b), q(q) { }
};

class ExprMesh : public TernaryExpr {
 private:
  Mesh mesh;
  Ref<SimplexTree<TV3,2>> tree;
 public:
  Vec transform_point (const Vec &p) const {
    for (size_t i = 0; i < elts.size()-1; i++)
      elts[i]->do_dist(p);
    return vec(a->d, b->d, q->d);
  }    
  void calc_dist (const Vec &p) {
    // if (p.x != a->d || p.y != b->d || p.z != q->d)
    //    printf("[%f,%f,%f] -> [%f,%f,%f] (%s, %s, %s)\n", p.x, p.y, p.z, a->d, b->d, q->d,
    //           a->to_str().c_str(), b->to_str().c_str(), q->to_str().c_str());
    auto np = transform_point(p);
    d = tree->distance(np);
    // d = tree->distance(p);
    try {
      d = tree->inside(np) ? -d : d;
    } catch (const ArithmeticError&) {
      d = d;
    }
  }
  virtual void do_dist (const Vec &p) {
    calc_dist(p);
  }
  virtual Expr* map (Expr* nx, Expr* ny, Expr* nz) {
    return new ExprMesh(mesh, a->map(nx, ny, nz), b->map(nx, ny, nz), q->map(nx, ny, nz));
  }
  virtual Expr* fold (void) { 
    auto fa = a->fold();
    auto fb = b->fold();
    auto fq = q->fold();
    auto res = new ExprMesh(mesh, fa, fb, fq);
    res->build();
    return res;
  }
  virtual void do_range (const Boxy &r) {
    auto p = r.center();
    calc_dist(p);
    Boxy nr(transform_point(r.lo), transform_point(r.hi));
    // printf("DISTANCE %f POINT [%f,%f,%f] -> [%f,%f,%f] BOX [%f,%f,%f]->[%f,%f,%f]\n",
    //        d, p.x, p.y, p.z, a->d, b->d, q->d, r.lo.x, r.lo.y, r.lo.z, r.hi.x, r.hi.y, r.hi.z);
    auto mdm = 0.5 * magnitude(nr.hi - nr.lo);
    // auto dms = r.dims();
    // auto mdm = 0.6 * sqrt2 * max(dms.x, max(dms.y, dms.z));
    if (d == INFTY)
      i = interval(-d, d);
    else 
      i = interval(d - mdm, d + mdm);
  }
  virtual std::string op_str (void) { return "Mesh"; }
  virtual std::string to_str (void) { return "Mesh(" + a->to_str() + "," + b->to_str() + "," + q->to_str() + ")"; }
 ExprMesh(Mesh mesh, Expr* a, Expr* b, Expr* q) : TernaryExpr(expr_mesh_kind, a, b, q), mesh(mesh), tree(new_<SimplexTree<TV3,2>>(mesh.soup, mesh.points, 4)) {  }
};

class ExprXform : public Expr {
 public:
  Expr* nx;
  Expr* ny;
  Expr* nz;
  Expr* g;
  virtual void do_dist (const Vec &p) { d = 0.0; }
  virtual Interval range (const Boxy &r) { return interval(0.0, 0.0); }
  virtual Expr* map (Expr* nnx, Expr* nny, Expr* nnz) { 
    return g->map(nx->map(nnx, nny, nnz), ny->map(nnx, nny, nnz), nz->map(nnx, nny, nnz)); 
  }
  virtual std::string op_str (void) { return "xform"; }
  virtual std::string to_str (void) { return "Xform(" + nx->to_str() + "," + ny->to_str() + "," + nz->to_str() + "," + g->to_str() + ")"; }
 ExprXform(Expr* nx, Expr* ny, Expr* nz, Expr* g) : Expr(expr_xform_kind), nx(nx), ny(ny), nz(nz), g(g) { }
};

#endif
