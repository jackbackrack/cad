#include "cad.h"
#include "expr.h"
#include "octree.h"
#include <fstream>

Expr* expr_lit(flo_t a) { return new ExprLit(a); }
Expr* expr_all(void) { 
  auto res = expr_lit(-INFTY); 
  res->xmin = expr_lit(-INFTY);
  res->ymin = expr_lit(-INFTY);
  res->zmin = expr_lit(-INFTY);
  res->xmax = expr_lit( INFTY);
  res->ymax = expr_lit( INFTY);
  res->zmax = expr_lit( INFTY);
  return res;
}
Expr* expr_none(void) { 
  auto res = expr_lit(INFTY); 
  res->xmin = expr_lit( INFTY);
  res->ymin = expr_lit( INFTY);
  res->zmin = expr_lit( INFTY);
  res->xmax = expr_lit(-INFTY);
  res->ymax = expr_lit(-INFTY);
  res->zmax = expr_lit(-INFTY);
  return res;
}
Expr* expr_pi(void) { return expr_lit(M_PI); }
Expr* expr_x(void) { return new ExprX(); }
Expr* expr_y(void) { return new ExprY(); }
Expr* expr_z(void) { return new ExprZ(); }
Expr* expr_neg(Expr* g) { return new ExprNeg(g); }
Expr* expr_not(Expr* g) { return expr_neg(g); }
Expr* expr_min (Expr* a, Expr* b) { 
  auto res = new ExprMin(a, b); 
  // printf("CHECKING A %p INF XMIN %p YMIN %p ZMIN %p XMAX %p YMAX %p ZMAX %p \n", 
  //        a, a->xmin, a->ymin, a->zmin, a->xmax, a->ymax, a->zmax);
  // printf("CHECKING B %p INF XMIN %p YMIN %p ZMIN %p XMAX %p YMAX %p ZMAX %p \n", 
  //        b, b->xmin, b->ymin, b->zmin, b->xmax, b->ymax, b->zmax);
  // if (b->zmax == NULL) {
  //   printf("BAD ZMAX\n");
  // }
  if (a->xmin != NULL && b->xmin != NULL)
  res->xmin = expr_min(a->xmin, b->xmin);
  if (a->ymin != NULL && b->ymin != NULL)
  res->ymin = expr_min(a->ymin, b->ymin);
  if (a->zmin != NULL && b->zmin != NULL)
  res->zmin = expr_min(a->zmin, b->zmin);
  if (a->xmax != NULL && b->xmax != NULL)
  res->xmax = expr_max(a->xmax, b->xmax);
  if (a->ymax != NULL && b->ymax != NULL)
  res->ymax = expr_max(a->ymax, b->ymax);
  if (a->zmax != NULL && b->zmax != NULL)
  res->zmax = expr_max(a->zmax, b->zmax);
  return res;
}
Expr* expr_or (Expr* a, Expr* b) { return expr_min(a, b); }
Expr* expr_max (Expr* a, Expr* b) { return new ExprMax(a, b); }
Expr* expr_and (Expr* a, Expr* b) { return expr_max(a, b); }
Expr* expr_add (Expr* a, Expr* b) { return new ExprAdd(a, b); }
Expr* expr_sub (Expr* a, Expr* b) { return new ExprSub(a, b); }
Expr* expr_mul (Expr* a, Expr* b) { return new ExprMul(a, b); }
Expr* expr_sqr (Expr* a) { return expr_mul(a, a); }
Expr* expr_div (Expr* a, Expr* b) { return new ExprDiv(a, b); }
Expr* expr_abs (Expr* a) { return expr_mul(a, a); }
Expr* expr_sqrt (Expr* a) { return new ExprSqrt(a); }
Expr* expr_sin (Expr* g) { return new ExprSin(g); }
Expr* expr_cos (Expr* g) { return new ExprCos(g); }
Expr* expr_tan (Expr* a) { return new ExprTan(a); }
Expr* expr_pow (Expr* a, Expr* b) { return new ExprPow(a, b); }
Expr* expr_asin (Expr* a) { return new ExprAsin(a); }
Expr* expr_acos (Expr* a) { return new ExprAcos(a); }
Expr* expr_atan (Expr* a) { return new ExprAtan(a); }

Expr* expr_xform(Expr* nx, Expr* ny, Expr* nz, Expr* g) { 
  return new ExprXform(nx, ny, nz, g);
}

Expr* expr_mov(Expr* dx, Expr* dy, Expr* dz, Expr* g) {
  auto tg = new ExprXform(expr_sub(expr_x(), dx), expr_sub(expr_y(), dy), expr_sub(expr_z(), dz), g);
  tg->copy_bbox(g);
  tg->xmin = expr_add(g->xmin, dx);
  tg->xmax = expr_add(g->xmax, dx);
  tg->ymin = expr_add(g->ymin, dy);
  tg->ymax = expr_add(g->ymax, dy);
  tg->zmin = expr_add(g->zmin, dz);
  tg->zmax = expr_add(g->zmax, dz);
  return tg;
}

Expr* expr_xmov(Expr* dx, Expr* g) {
  return expr_mov(dx,expr_lit(0.0),expr_lit(0.0),g);
}

Expr* expr_ymov(Expr* dy, Expr* g) {
  return expr_mov(expr_lit(0.0),dy,expr_lit(0.0),g);
}

Expr* expr_zmov(Expr* dz, Expr* g) {
  return expr_mov(expr_lit(0.0),expr_lit(0.0),dz,g);
}

Expr* expr_mag(Expr* dx, Expr* dy, Expr* dz, Expr* g) {
  auto tg = new ExprXform(expr_div(expr_x(), dx), expr_div(expr_y(), dy), expr_div(expr_z(), dz), g);
  tg->xmin = expr_mul(dx, g->xmin);
  tg->xmax = expr_mul(dx, g->xmax);
  tg->ymin = expr_mul(dy, g->ymin);
  tg->ymax = expr_mul(dy, g->ymax);
  tg->zmin = expr_mul(dz, g->zmin);
  tg->zmax = expr_mul(dz, g->zmax);
  return tg;
}

Expr* expr_mag1(Expr* dxyz, Expr* g) {
  auto tg = new ExprXform(expr_div(expr_x(), dxyz), expr_div(expr_y(), dxyz), expr_div(expr_z(), dxyz), g);
  tg->xmin = expr_mul(dxyz, g->xmin);
  tg->xmax = expr_mul(dxyz, g->xmax);
  tg->ymin = expr_mul(dxyz, g->ymin);
  tg->ymax = expr_mul(dxyz, g->ymax);
  tg->zmin = expr_mul(dxyz, g->zmin);
  tg->zmax = expr_mul(dxyz, g->zmax);
  return tg;
}

Expr* expr_xmag(Expr* dx, Expr* g) {
  auto tg = expr_mag(dx,expr_lit(1.0),expr_lit(1.0),g);
  tg->copy_bbox(g);
  tg->xmin = expr_mul(dx, g->xmin);
  tg->xmax = expr_mul(dx, g->xmax);
  return tg;
}

Expr* expr_ymag(Expr* dy, Expr* g) {
  auto tg = expr_mag(expr_lit(1.0),dy,expr_lit(1.0),g);
  tg->copy_bbox(g);
  tg->ymin = expr_mul(dy, g->ymin);
  tg->ymax = expr_mul(dy, g->ymax);
  return tg;
}

Expr* expr_zmag(Expr* dz, Expr* g) {
  auto tg = expr_mag(expr_lit(1.0),expr_lit(1.0),dz,g);
  tg->copy_bbox(g);
  tg->zmin = expr_mul(dz, g->zmin);
  tg->zmax = expr_mul(dz, g->zmax);
  return tg;
}

Expr* expr_yrevolve(Expr* g) {
  return new ExprXform(expr_sqrt(expr_add(expr_sqr(expr_x()),expr_sqr(expr_z()))),expr_y(),expr_z(),g);
}

Expr* expr_xrevolve(Expr* g) {
  return new ExprXform(expr_x(),expr_sqrt(expr_add(expr_sqr(expr_y()),expr_sqr(expr_z()))),expr_z(),g);
}

Expr* expr_shear_x_z(Expr* z0, Expr* z1, Expr* dx0, Expr* dx1, Expr* g) {
  return new ExprXform(expr_sub(expr_x(),
                                expr_sub(dx0, expr_div(expr_mul(expr_sub(dx1,dx0),expr_sub(expr_z(),z0)),
                                                       expr_sub(z1,z0)))),
                       expr_y(),
                       expr_z(),
                       g);
}

Expr* expr_taper(Expr* dz, Expr* s0, Expr* s1, Expr* g) {
  auto zrad = expr_mul(expr_lit(0.5), dz);
  auto z0   = expr_neg(zrad);
  auto z1   = zrad;
  return new ExprXform(expr_mul(expr_x(), expr_div(dz,expr_add(expr_mul(s1,expr_sub(expr_z(),z0)), expr_mul(s0,expr_sub(z1,expr_z()))))),
                       expr_mul(expr_y(), expr_div(dz,expr_add(expr_mul(s1,expr_sub(expr_z(),z0)), expr_mul(s0,expr_sub(z1,expr_z()))))),
                       expr_z(),
                       g);
}

Expr* expr_xrot(Expr* a, Expr* g) {
  auto sa = expr_sin(a);
  auto ca = expr_cos(a);
  return new ExprXform(expr_x(), expr_add(expr_mul(ca, expr_y()), expr_mul(sa, expr_z())), expr_add(expr_mul(expr_neg(sa), expr_y()), expr_mul(ca, expr_z())), g);
}

Expr* expr_reflect_x(Expr* g) {
  auto tg = new ExprXform(expr_neg(expr_x()), expr_y(), expr_z(), g);
  tg->copy_bbox(g);
  tg->xmin = expr_neg(g->xmax);
  tg->xmax = expr_neg(g->xmin);
  return tg;
}

Expr* expr_reflect_xy(Expr* g) {
  auto tg = new ExprXform(expr_y(), expr_x(), expr_z(), g);
  tg->copy_bbox(g);
  tg->xmin = tg->ymin;
  tg->xmax = tg->ymax;
  tg->ymin = tg->xmin;
  tg->ymax = tg->xmax;
  return tg;
}

Expr* expr_yrot(Expr* a, Expr* g) {
  auto sa = expr_sin(a);
  auto ca = expr_cos(a);
  return new ExprXform(expr_add(expr_mul(ca, expr_x()), expr_mul(sa, expr_z())), expr_y(), expr_add(expr_mul(expr_neg(sa), expr_x()), expr_mul(ca, expr_z())), g);
}

Expr* expr_reflect_y(Expr* g) {
  auto tg = new ExprXform(expr_x(), expr_neg(expr_y()), expr_z(), g);
  tg->copy_bbox(g);
  tg->ymin = expr_neg(g->ymax);
  tg->ymax = expr_neg(g->ymin);
  return tg;
}

Expr* expr_zrot(Expr* a, Expr* g) {
  auto sa = expr_sin(a);
  auto ca = expr_cos(a);
  return new ExprXform(expr_add(expr_mul(ca, expr_x()), expr_mul(sa, expr_y())), expr_add(expr_mul(expr_neg(sa), expr_x()), expr_mul(ca, expr_y())), expr_z(), g);
}

Expr* expr_reflect_z(Expr* g) {
  auto tg = new ExprXform(expr_x(), expr_y(), expr_neg(expr_z()), g);
  tg->copy_bbox(g);
  tg->zmin = expr_neg(g->zmax);
  tg->zmax = expr_neg(g->zmin);
  return tg;
}

Expr* expr_reflect_xz(Expr* g) {
  auto tg = new ExprXform(expr_z(), expr_y(), expr_x(), g);
  tg->copy_bbox(g);
  tg->xmin = g->zmin;
  tg->xmax = g->zmax;
  tg->zmin = g->xmin;
  tg->zmax = g->xmax;
  return tg;
}

Expr* expr_reflect_yz(Expr* g) {
  auto tg = new ExprXform(expr_x(), expr_z(), expr_y(), g);
  tg->copy_bbox(g);
  tg->ymin = g->zmin;
  tg->ymax = g->zmax;
  tg->zmin = g->ymin;
  tg->zmax = g->ymax;
  return tg;
}

Expr* expr_half(Expr* nx, Expr* ny, Expr* nz, Expr* d) {
  return expr_sub(d, expr_add(expr_mul(expr_x(), nx), expr_add(expr_mul(expr_y(), ny), expr_mul(expr_z(), nz))));
}

Expr* expr_edge(Expr* x0, Expr* y0, Expr* x1, Expr* y1) {
  auto dx = expr_sub(x1,x0);
  auto dy = expr_sub(y1,y0);
  auto res = expr_sub(expr_mul(dy, expr_sub(expr_x(),x0)), expr_mul(dx, expr_sub(expr_y(),y0)));
  // printf("EDGE %s\n", res->to_str().c_str());
  return res;
}

Expr* expr_triangle(Expr* x0, Expr* y0, Expr* x1, Expr* y1, Expr* x2, Expr* y2) {
  auto e0 = expr_edge(x1,y1,x0,y0);
  auto e1 = expr_edge(x2,y2,x1,y1);
  auto e2 = expr_edge(x0,y0,x2,y2);
  auto g  = expr_max(e0, expr_max(e1, e2));
  g->xmin = expr_min(x0, expr_min(x1, x2));
  g->xmax = expr_max(x0, expr_max(x1, x2));
  g->ymin = expr_min(y0, expr_min(y1, y2));
  g->ymax = expr_max(y0, expr_max(y1, y2));
  return g;
}

Expr* expr_square(Expr* d) {
  auto r = expr_mul(expr_lit(0.5), d);
  auto g = expr_and(expr_and(expr_sub(expr_neg(r),expr_x()), expr_sub(expr_x(),r)), 
                    expr_and(expr_sub(expr_neg(r),expr_y()), expr_sub(expr_y(),r)));
  g->xmin = expr_neg(r);
  g->xmax = r;
  g->ymin = expr_neg(r);
  g->ymax = r;
  return g;
}

Expr* expr_rem(Expr* a, Expr* b) {
  return expr_and(a, expr_neg(b));
}

Expr* expr_hollow(Expr* d, Expr* a) {
  return expr_rem(a, expr_add(a, d));
}

Expr* expr_circle(Expr* d) {
  auto r = expr_mul(expr_lit(0.5), d);
  auto g = expr_sub(expr_sqrt(expr_add(expr_sqr(expr_x()), expr_sqr(expr_y()))),r);
  g->xmin = expr_neg(r);
  g->xmax = r;
  g->ymin = expr_neg(r);
  g->ymax = r;
  return g;
}

Expr* expr_sphere(Expr* d) {
  auto r = expr_mul(expr_lit(0.5), d);
  auto g = expr_sub(expr_sqrt(expr_add(expr_sqr(expr_x()), expr_add(expr_sqr(expr_y()), expr_sqr(expr_z())))),r);
  g->xmin = expr_neg(r);
  g->xmax = r;
  g->ymin = expr_neg(r);
  g->ymax = r;
  g->zmin = expr_neg(r);
  g->zmax = r;
  return g;
}

Expr* expr_extrude(Expr* d, Expr* g) {
  auto r  = expr_mul(expr_lit(0.5), d);
  auto lo = expr_neg(r);
  auto hi = r;
  auto tg = expr_and(expr_and(expr_sub(lo,expr_z()), expr_sub(expr_z(),r)), g);
  tg->xmin = g->xmin;
  tg->xmax = g->xmax;
  tg->ymin = g->ymin;
  tg->ymax = g->ymax;
  tg->zmin = lo;
  tg->zmax = hi;
  // printf("CHECKING EXT TG %p INF XMIN %p YMIN %p ZMIN %p XMAX %p YMAX %p ZMAX %p \n", 
  //        tg, tg->xmin, tg->ymin, tg->zmin, tg->xmax, tg->ymax, tg->zmax);
  return tg;
}

Expr* expr_xbox(Expr* a, Expr* b) {
  Expr* wa = expr_sub(a->xmax,a->xmin);
  Expr* wb = expr_sub(b->xmax,b->xmin);
  Expr* r = expr_mul(expr_lit(0.5),expr_add(wa,wb));
  auto tg = expr_or(expr_xmov(expr_sub(expr_neg(r),a->xmin),a),expr_xmov(expr_sub(expr_add(expr_neg(r),wa),b->xmin),b));
  tg->xmin = expr_neg(r);
  tg->xmax = r;
  tg->ymin = expr_min(a->ymin, b->ymin);
  tg->ymax = expr_max(a->ymax, b->ymax);
  if (a->zmin != NULL && b->zmin != NULL) {
    tg->zmin = expr_min(a->zmin, b->zmin);
    tg->zmax = expr_max(a->zmax, b->zmax);
  }
  return tg;
}

Expr* expr_align_xmin(Expr* a, Expr* b) {
  auto del = expr_sub(a->xmin,b->xmin);
  auto tg  = expr_xmov(del,b);
  tg->copy_bbox(b);
  tg->xmin = expr_add(b->xmin, del);
  tg->xmax = expr_add(b->xmax, del);
  return tg;
}

Expr* expr_align_xmax(Expr* a, Expr* b) {
  auto del = expr_sub(a->xmax,b->xmax);
  auto tg  = expr_xmov(del,b);
  tg->copy_bbox(b);
  tg->xmin = expr_add(b->xmin, del);
  tg->xmax = expr_add(b->xmax, del);
  return tg;
}

Expr* expr_align_ymin(Expr* a, Expr* b) {
  auto del = expr_sub(a->ymin,b->ymin);
  auto tg  = expr_ymov(del,b);
  tg->copy_bbox(b);
  tg->ymin = expr_add(b->ymin, del);
  tg->ymax = expr_add(b->ymax, del);
  return tg;
}

Expr* expr_align_ymax(Expr* a, Expr* b) {
  auto del = expr_sub(a->ymax,b->ymax);
  auto tg  = expr_ymov(del,b);
  tg->copy_bbox(b);
  tg->ymin = expr_add(b->ymin, del);
  tg->ymax = expr_add(b->ymax, del);
  return tg;
}

Expr* expr_align_zmin(Expr* a, Expr* b) {
  auto del = expr_sub(a->zmin,b->zmin);
  auto tg  = expr_zmov(del,b);
  tg->copy_bbox(b);
  tg->zmin = expr_add(b->zmin, del);
  tg->zmax = expr_add(b->zmax, del);
  return tg;
}

Expr* expr_align_zmax(Expr* a, Expr* b) {
  auto del = expr_sub(a->zmax,b->zmax);
  auto tg  = expr_zmov(del,b);
  tg->copy_bbox(b);
  tg->zmin = expr_add(b->zmin, del);
  tg->zmax = expr_add(b->zmax, del);
  return tg;
}

Expr* expr_ybox(Expr* a, Expr* b) {
  Expr* ha = expr_sub(a->ymax,a->ymin);
  Expr* hb = expr_sub(b->ymax,b->ymin);
  Expr* r = expr_mul(expr_lit(0.5),expr_add(ha,hb));
  auto tg = expr_or(expr_ymov(expr_sub(expr_neg(r),a->ymin),a),expr_ymov(expr_sub(expr_add(expr_neg(r),ha),b->ymin),b));
  tg->ymin = expr_neg(r);
  tg->ymax = r;
  tg->xmin = expr_min(a->xmin, b->xmin);
  tg->xmax = expr_max(a->xmax, b->xmax);
  if (a->zmin != NULL && b->zmin != NULL) {
    tg->zmin = expr_min(a->zmin, b->zmin);
    tg->zmax = expr_max(a->zmax, b->zmax);
  }
  return tg;
}

Expr* expr_zbox(Expr* a, Expr* b) {
  Expr* za = expr_sub(a->zmax,a->zmin);
  Expr* zb = expr_sub(b->zmax,b->zmin);
  Expr* r  = expr_mul(expr_lit(0.5),expr_add(za,zb));
  auto tg  = expr_or(expr_zmov(expr_sub(expr_neg(r),a->zmin),a),
                     expr_zmov(expr_sub(expr_add(expr_neg(r),za),b->zmin),b));
  tg->zmin = expr_neg(r);
  tg->zmax = r;
  tg->xmin = expr_min(a->xmin, b->xmin);
  tg->xmax = expr_max(a->xmax, b->xmax);
  tg->ymin = expr_min(a->ymin, b->ymin);
  tg->ymax = expr_max(a->ymax, b->ymax);
  return tg;
}

Expr* expr_cylinder(Expr* xyrad, Expr* zrad) {
  return expr_extrude(zrad, expr_circle(xyrad));
}

Expr* expr_capsule(Expr* xyrad, Expr* zrad) {
  return expr_or(expr_cylinder(xyrad, zrad),
                 expr_or(expr_zmov(expr_neg(zrad), expr_sphere(xyrad)), expr_zmov(zrad, expr_sphere(xyrad))));
}

Expr* expr_cube(Expr* d) {
  return expr_extrude(d, expr_square(d));
}

Expr* expr_space(Expr* d) {
  auto r = expr_mul(expr_lit(0.5), d);
  auto g = new ExprSpace();
  g->xmin = expr_neg(r);
  g->xmax = r;
  g->ymin = expr_neg(r);
  g->ymax = r;
  g->zmin = expr_neg(r);
  g->zmax = r;
  return g;
}

Expr* expr_rect2(Expr* xd, Expr* yd) {
  auto xrad = expr_mul(expr_lit(0.5), xd);
  auto yrad = expr_mul(expr_lit(0.5), yd);
  auto g = expr_and(expr_and(expr_sub(expr_neg(xrad),expr_x()), expr_sub(expr_x(),xrad)), 
                    expr_and(expr_sub(expr_neg(yrad),expr_y()), expr_sub(expr_y(),yrad)));
  g->xmin = expr_neg(xrad);
  g->xmax = xrad;
  g->ymin = expr_neg(yrad);
  g->ymax = yrad;
  return g;
}

Expr* expr_rect3(Expr* xd, Expr* yd, Expr* zd) {
  return expr_extrude(zd, expr_rect2(xd, yd));
}

Expr* expr_cone(Expr* xyd, Expr* zd) {
  auto cyl = expr_cylinder(xyd, zd);
  return expr_taper(zd, expr_lit(1.0), expr_lit(0.0), cyl);
}

Expr* expr_pyramid(Expr* xyd, Expr* zd) {
  auto c = expr_rect3(xyd, xyd, zd);
  return expr_taper(zd, expr_lit(1.0), expr_lit(0.0), c);
}

Expr* expr_blend(Expr* g1, Expr* g2, Expr* amount) {
  auto j = expr_add(g1, g2);
  auto f = expr_sub(expr_add(expr_sqrt(expr_abs(g1)), expr_sqrt(expr_abs(g2))), amount);
  return expr_add(j, f);
}

Expr* polygonize (Nested<TV2> p) {
  Mesh mesh = triangulate(p);
  // printf("%d VERTICES %d INDICES\n", (int)mesh.vertices.size(), (int)mesh.indices.size());
  Expr* res = expr_none();
  Boxy box(vec( INFTY,  INFTY,  INFTY), vec(-INFTY, -INFTY, -INFTY));
  for (auto tri : mesh.soup->elements) {
    Vec p0 = mesh.points[tri.x];
    Vec p1 = mesh.points[tri.y];
    Vec p2 = mesh.points[tri.z];
    box = add(add(add(box, p0), p1), p2);
    res = expr_or(res, expr_triangle(expr_lit(p0.x), expr_lit(p0.y), expr_lit(p1.x), expr_lit(p1.y), expr_lit(p2.x), expr_lit(p2.y)));
  }
  res->xmin = expr_lit(box.lo.x);
  res->ymin = expr_lit(box.lo.y);
  res->zmin = expr_lit(box.lo.z);
  res->xmax = expr_lit(box.hi.x);
  res->ymax = expr_lit(box.hi.y);
  res->zmax = expr_lit(box.hi.z);
  return res;
}

Expr* expr_mesh(Mesh mesh) {
  return new ExprMesh(mesh, expr_x(), expr_y(), expr_z());
}

Expr* pretty_print_expr(Expr* e) {
  auto s = e->to_str();
  printf("%s\n", s.c_str());
  return e;
}
