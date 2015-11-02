#include "cad.h"
#include <geode/geometry/polygon.h>
#include <geode/geometry/offset_mesh.h>
#include <geode/exact/mesh_csg.h>
// #include <geode/exact/delaunay.h>
#include <geode/exact/polygon_csg.h>
#include <geode/mesh/SegmentSoup.h>
#include <geode/mesh/TriangleTopology.h>
#include <geode/mesh/improve_mesh.h>

void error (std::string msg) {
  fprintf(stderr, "MSG %s\n", msg.c_str());
  exit(-1);
}

void error (std::string msg, std::string arg) {
  fprintf(stderr, "MSG %s %s\n", msg.c_str(), arg.c_str());
  exit(-1);
}

void ensure (bool val, std::string msg) {
  if (!val) error(msg);
}

Mesh fab_mesh (Array<IV> faces, Array<TV> points) {
  return tuple(new_<const TriangleSoup>(faces),points);
}

typedef real T;
typedef Vector<T,2> TV2;
typedef Vector<T,3> TV;
typedef Vector<int,3> IV;
typedef Vector<int,2> IV2;

Matrix<T,4> to_matrix44(Matrix<T,3> M) {
  return Matrix<T,4>(M.x[0][0],M.x[1][0],M.x[2][0],0,M.x[0][1],M.x[1][1],M.x[2][1],0,M.x[0][2],M.x[1][2],M.x[2][2],0,0,0,0,1);
}
Matrix<T,4> rotation_matrix(TV from, TV to) {
  auto m = Rotation<Vector<T,3> >::from_rotated_vector(from, to).matrix();
  return to_matrix44(m);
}
Matrix<T,4> rotation_matrix(TV2 from, TV2 to) {
  auto a = atan2(to.y,to.x) - atan2(from.y,from.x);
  T c=cos(a),s=sin(a);
  return Matrix<T,4>(c,s,0,0,-s,c,0,0,0,0,1,0,0,0,0,1);
}
Matrix<T,4> scale_matrix(TV v) {
  return Matrix<T,4>(v.x,0,0,0,0,v.y,0,0,0,0,v.z,0,0,0,0,1);
}
Matrix<T,4> translation_matrix(TV v) {
  return Matrix<T,4>(1,0,0,0,0,1,0,0,0,0,1,0,v.x,v.y,v.z,1);
}
Matrix<T,4> translation_matrix(TV2 v) {
  return Matrix<T,4>(1,0,0,0,0,1,0,0,0,0,1,0,v.x,v.y,0,1);
}

double is_clockwise (Array<TV2> contour) {
  auto area = polygon_area(RawArray<TV2>(contour));
  return area < 0;
}

// Combines multiple triangle soups into one
static Mesh concat_soups(Mesh soup0, Mesh soup1) {
  Array<IV> triangles;
  Array<TV> vertices;

  for(const auto& s : { soup0, soup1 }) {
    // 0th vertex in each soup will be first vertex after all previous ones
    // Save as a vector of ints to we can directly add it to each triangle
    const auto index_offset = IV::repeat(vertices.size());
    vertices.extend(s.y);
    for(const auto& t : s.x->elements) {
      triangles.append(t+index_offset); // Add new triangle with indices mapped to combined vertices
    }
  }

  return fab_mesh(triangles,vertices);
}

Ref<const TriangleSoup> const_soup(Ref<TriangleSoup> val) {
  return new_<const TriangleSoup>(val->elements);
}

Mesh const_soup(Tuple<Ref<TriangleSoup>, Array<TV>> val) {
  return tuple(const_soup(val.x), val.y);
}

Mesh simplify_mesh(Mesh mesh) {
  printf("SIMPLIFYING\n");
  Array<IV> faces;
  for (auto face : mesh.x->elements)
    faces.append(face);
  auto points = mesh.y.copy();
  Array<int> mapping;
  for (int i = 0; i < points.size(); i++)
    mapping.append(i);
  for (;;) {
    bool is_changed = false;
    for (auto face : faces) {
      if (face.x != face.y && points[face.x] == points[face.y]) {
        printf("MAPPING %d TO %d\n", face.x, face.y);
        mapping[face.x] = face.y; is_changed = true; break;
      }
      if (face.y != face.z && points[face.y] == points[face.z]) {
        printf("MAPPING %d TO %d\n", face.y, face.z);
        mapping[face.y] = face.z; is_changed = true; break;
      }
      if (face.x != face.z && points[face.x] == points[face.z]) {
        printf("MAPPING %d TO %d\n", face.x, face.z);
        mapping[face.x] = face.z; is_changed = true; break;
      }
    }
    if (!is_changed) break;
    for (int i = 0; i < faces.size(); i++) {
      auto face = faces[i];
      if (face.x < 0 || face.x >= points.size())
        printf("BAD FACE X %d\n", face.x);
      if (face.y < 0 || face.y >= points.size())
        printf("BAD FACE X %d\n", face.y);
      if (face.z < 0 || face.z >= points.size())
        printf("BAD FACE X %d\n", face.z);
      faces[i] = vec(mapping[face.x], mapping[face.y], mapping[face.z]);
    }
  }
  Array<IV> new_faces;
  for (auto face : mesh.x->elements) {
    if (face.x != face.y && face.x != face.z && face.y != face.z) {
      new_faces.append(face);
    } else
      printf("REMOVED FACE %d,%d,%d\n", face.x, face.y, face.z);
  }
  printf("%d FACES NOW %d REMOVED %d FACES\n",
         faces.size(), new_faces.size(), faces.size() - new_faces.size() );
  return fab_mesh(new_faces, points);
}

/*
Mesh report_simplify_mesh(Mesh mesh) {
  Array<IV> new_faces;
  auto points = mesh.y;
  for (auto face : mesh.x->elements) {
    if (points[face.x] == points[face.y])
      printf("zero edge %d -> %d\n", face.x, face.y);
    if (points[face.y] == points[face.z])
      printf("zero edge %d -> %d\n", face.y, face.z);
    if (points[face.x] == points[face.z])
      printf("zero edge %d -> %d\n", face.x, face.z);
    new_faces.append(face);
  }
  printf("%d FACES NOW %d REMOVED %d FACES\n",
         mesh.x->elements.size(), new_faces.size(), mesh.x->elements.size() - new_faces.size() );
  return fab_mesh(new_faces, points);
}

Mesh simplify_mesh(Mesh mesh) {
  // printf("STARTING SIMPLIFICATION\n");
  // report_simplify_mesh(mesh);
  Array<TV> pos(mesh.y);
  Field<TV,VertexId> field(pos.copy());
  auto topo = new_<MutableTriangleTopology>();
  topo->add_vertices(mesh.y.size());
  for (auto face : mesh.x->elements)
    topo->add_face(vec((VertexId)face.x, (VertexId)face.y, (VertexId)face.z));
  //  ImproveOptions options(1.1,0.01,0.01);
  ImproveOptions options(1.0000001,0.0000001,0.0000001);
  improve_mesh_inplace(topo, field, options);
  // printf("BEFORE GC %d\n", field.size());
  auto updates = topo->collect_garbage();
  // printf("AFTER GC %d UPDATES %d\n", field.size(), updates.x.size());
  int i = 0, tot = 0;
  for (auto update : updates.x) {
    if (update >= 0) tot += 1;
    // printf("  %d -> %d\n", i, update);
    i += 1;
  }
  Array<TV> new_points(tot);
  i = 0;
  for (auto update : updates.x) {
    if (update >= 0) new_points[update] = pos[i];
    i += 1;
  }
  auto new_soup = topo->face_soup().x;
  printf("BEFORE %d,%d AFTER %d,%d\n",
         pos.size(), mesh.x->elements.size(), new_points.size(), new_soup->elements.size());
  return tuple(const_soup(new_soup), new_points);
}
*/

// invert triangle soup so normals point inwards
Mesh invert_soup(Mesh soup) {
  Array<IV> triangles;

  for(const auto& t : soup.x->elements) {
    triangles.append(vec(t[0], t[2], t[1]));
  }

  return fab_mesh(triangles, soup.y);
}

Nested<TV2> invert_poly(Nested<TV2> poly) {
  Nested<TV2,false> pres;
  for (auto contour : poly) {
    Array<TV2> cres;
    for (int i = contour.size()-1; i >= 0; i--)
      cres.append(contour[i]);
    pres.append(cres);
  }
  pres.freeze();
  return pres;
}

Array<TV2> simplify_contour(RawArray<TV2> contour) {
  // TODO: LIMITS
  return polygon_simplify(contour, 179.9, 0.0001);
}

Nested<TV2> simplify_poly(Nested<TV2> poly) {
  Nested<TV2,false> res;
  for (auto contour : poly) 
    res.append(simplify_contour(contour));
  res.freeze();
  return res;
}

Nested<TV2> add(Nested<TV2> c0, Nested<TV2> c1) {
  auto res = polygon_union(c0, c1);
  // printf("UNION POLY\n");
  if (true)
    return simplify_poly(res);
  else
    return res;
}

Nested<TV2> mul(Nested<TV2> c0, Nested<TV2> c1) {
  return polygon_intersection(c0, c1);
}

Nested<TV2> sub(Nested<TV2> c0, Nested<TV2> c1) {
  return polygon_union(c0, invert_poly(c1));
}

Nested<TV2> offset(T a, Nested<TV2> c) {
  // auto merge = concat_polygons(c0, invert_polygon(c1));
  // return split_polygons(merge.x, merge.y, 0);
  return c;
}

Mesh add(Mesh soup0, Mesh soup1) {
  auto merge = concat_soups(soup0, soup1);
  // printf("--- START UNIONING ---\n");
  auto res   = split_soup(merge.x, merge.y, 0);
  // printf("--- DONE  UNIONING ---\n");
  res = simplify_mesh(res);
  return res;
}

void pretty_print_soup(Mesh soup) {
  int i = 0;
  for (auto pt : soup.y) {
    printf("PT[%2d] %g,%g,%g\n", i, pt.x, pt.y, pt.z);
    i += 1;
  }
  for (auto tri : soup.x->elements)
    printf("TRI %2d,%2d,%2d\n", tri.x, tri.y, tri.z);
}

void pretty_print_line2(Array<TV2> line) {
  int i = 0;
  for (auto pt : line) {
    printf("PT[%2d] %g,%g\n", i, pt.x, pt.y);
    i += 1;
  }
}

void pretty_print_line3(Array<TV> line) {
  int i = 0;
  for (auto pt : line) {
    printf("PT[%2d] %g,%g,%g\n", i, pt.x, pt.y, pt.z);
    i += 1;
  }
}

void pretty_print_contour(Array<TV2> contour) {
  pretty_print_line2(contour);
}

void pretty_print_poly(Nested<TV2> poly) {
  int i = 0;
  for (auto elt : poly) {
    Array<TV2> contour; for (auto e : elt) contour.append(e);
    printf("CONTOUR %d\n", i);
    pretty_print_line2(contour);
    i += 1;
  }
}

void pretty_print_matrix(Matrix<T,4> M) {
  printf("%g, %g, %g, %g\n", M.x[0][0],M.x[1][0],M.x[2][0],M.x[3][0]);
  printf("%g, %g, %g, %g\n", M.x[0][1],M.x[1][1],M.x[2][1],M.x[3][1]);
  printf("%g, %g, %g, %g\n", M.x[0][2],M.x[1][2],M.x[2][2],M.x[3][2]);
  printf("%g, %g, %g, %g\n", M.x[0][3],M.x[1][3],M.x[2][3],M.x[3][3]);
}

void pretty_print_polyline3(Nested<TV> polyline) {
  int i = 0;
  for (auto elt : polyline) {
    Array<TV> line; for (auto e : elt) line.append(e);
    printf("LINE %d\n", i);
    pretty_print_line3(line);
    i += 1;
  }
}

void pretty_print_polyline2(Nested<TV2> polyline) {
  int i = 0;
  for (auto elt : polyline) {
    Array<TV2> line; for (auto e : elt) line.append(e);
    printf("LINE %d\n", i);
    pretty_print_line2(line);
    i += 1;
  }
}

void do_print_line2(std::string name, Array<TV2> line) {
  int i = 0;
  printf("%s(", name.c_str());
  for (auto pt : line) {
    if (i > 0) printf(", ");
    printf("vec(%g,%g)", pt.x, pt.y);
    i += 1;
  }
  printf(")");
}

void do_print_line3(std::string name, Array<TV> line) {
  int i = 0;
  printf("%s(", name.c_str());
  for (auto pt : line) {
    if (i > 0) printf(", ");
    printf("vec(%g,%g,%g)", pt.x, pt.y, pt.z);
    i += 1;
  }
  printf(")");
}

void do_print_faces(std::string name, Array<const IV> line) {
  int i = 0;
  printf("%s(", name.c_str());
  for (auto pt : line) {
    if (i > 0) printf(", ");
    printf("vec(%d,%d,%d)", pt.x, pt.y, pt.z);
    i += 1;
  }
  printf(")");
}

void print_soup(Mesh soup) {
  printf("mesh(");
  do_print_line3("points", soup.y);
  printf(", ");
  do_print_faces("faces", soup.x->elements);
  printf(")\n");
}

void print_contour(Array<TV2> contour) {
  do_print_line2("contour", contour);
}

void print_line3(Array<TV> contour) {
  do_print_line3("line3", contour);
}

void print_line2(Array<TV2> contour) {
  do_print_line2("line2", contour);
}

void print_poly(Nested<TV2> poly) {
  int i = 0;
  printf("poly(");
  for (auto elt : poly) {
    Array<TV2> contour; for (auto e : elt) contour.append(e);
    Array<TV2> line; for (auto e : elt) line.append(e);
    if (i > 0) printf(", ");
    print_contour(line);
    i += 1;
  }
  printf(")\n");
}

void print_polyline3(Nested<TV> polyline) {
  int i = 0;
  printf("polyline3(");
  for (auto elt : polyline) {
    Array<TV> line; for (auto e : elt) line.append(e);
    if (i > 0) printf(", ");
    print_line3(line);
    i += 1;
  }
  printf(")\n");
}

void print_polyline2(Nested<TV2> polyline) {
  int i = 0;
  printf("polyline2(");
  for (auto elt : polyline) {
    Array<TV2> line; for (auto e : elt) line.append(e);
    if (i > 0) printf(", ");
    print_line2(line);
    i += 1;
  }
  printf(")\n");
}

void print_matrix(Matrix<T,4> M) {
  printf("mat(%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g)\n",
         M.x[0][0],M.x[1][0],M.x[2][0],M.x[3][0],
         M.x[0][1],M.x[1][1],M.x[2][1],M.x[3][1],
         M.x[0][2],M.x[1][2],M.x[2][2],M.x[3][2],
         M.x[0][3],M.x[1][3],M.x[2][3],M.x[3][3]);
}

Mesh mul(Mesh soup0, Mesh soup1) {
  auto merge = concat_soups(soup0, soup1);
  auto res = split_soup(merge.x, merge.y, 1);
  return res;
}

Mesh mesh_from(int start, Mesh soup) {
  Array<TV> pts;
  for (int i = start; i < soup.y.size(); i++) {
    pts.append(soup.y[i]);
  }
  /*
  for (int i = 0; i < pts.size(); i++) {
    for (int j = i + 1; j < pts.size(); j++) {
      if (pts[i] == pts[j]) {
        printf("%d (%f,%f,%f) === %d (%f,%f,%f)\n", i, pts[i].x, pts[i].y, pts[i].z, j, pts[j].x, pts[j].y, pts[j].z);
      }
    }
  }
  */
  Array<IV> faces;
  for (auto tri : soup.x->elements) {
    bool is_all = tri.x >= start && tri.y >= start && tri.z >= start;
    if (is_all)
      faces.append(vec(tri.x - start, tri.y - start, tri.z - start));
  }
  return tuple(new_<const TriangleSoup>(faces), pts);
}

Nested<TV2> slice(T z, Mesh soup) {
  T t = 1e8;
  auto smesh = const_soup(cube_mesh(vec(-t, -t, -t), vec( t,  t,  z)));
  // print_soup(smesh);
  auto res = mul(smesh, soup);
  int start = smesh.y.size() + soup.y.size();
  /*
  printf("---->>>>\n");
  print_soup(res);
  printf("========\n");
  printf("%d OLD %d NEW VERTICES\n", start, res.y.size() - start);
  for (auto tri : res.x->elements) {
    bool is_all = tri.x >= start && tri.y >= start && tri.z >= start;
    if (is_all)
      printf("SLICE %2d,%2d,%2d\n", tri.x, tri.y, tri.z);
  }
  */
  auto new_mesh = mesh_from(start, res);
  auto boundary = new_mesh.x->boundary_mesh();
  auto polys = boundary->polygons();
  /*
  int i = 0;
  for (auto p : new_mesh.y) {
    printf("BP[%2d] %f,%f,%f\n", i, p.x, p.y, p.z);
    i += 1;
  }
  for (auto elt : boundary->elements) {
    printf("ELT %d,%d\n", elt.x, elt.y);
  }
  printf("POLYGONS\n");
  for (auto poly : polys.x) {
    for (auto p : poly) 
      printf("%d ", p);
    printf("\n");
  }
  */

  Nested<TV2,false> pres;
  for (auto poly : polys.x) {
    Array<TV2> contour;
    for (auto p : poly) {
      auto pt = new_mesh.y[p];
      contour.append(vec(pt.x, pt.y));
    }
    pres.append(contour);
  }
  pres.freeze();
  return pres;
}

Mesh sub(Mesh soup0, Mesh soup1) {
  auto merge = concat_soups(soup0, invert_soup(soup1));
  return split_soup(merge.x, merge.y, 0);
}

template<class ET>
Nested<ET> contour_to_poly(Array<ET> contour) {
  Nested<ET,false> poly;
  poly.append(contour);
  poly.freeze();
  return poly;
}

template<class ET>
Array<ET> poly_to_contour(Nested<ET> poly, int i) {
  Array<ET> contour;
  printf("POLY TO CONTOUR[%d] %d\n", i, poly[i].size());
  for (auto e : poly[i])
    contour.append(e);
  return contour;
}

Mesh all_soup(void) {
  T x = 1e8;
  return const_soup(cube_mesh(vec(-x, -x, -x), vec(x, x, x)));
}

Mesh none_soup(void) {
  Array<IV> faces;
  Array<TV> points;
  return fab_mesh(faces, points);
}

Nested<TV2> square_poly(TV2 min, TV2 max) {
  Array<Vector<real, 2>> pts;
  pts.append(vec(max.x,min.y));
  pts.append(vec(max.x,max.y));
  pts.append(vec(min.x,max.y));
  pts.append(vec(min.x,min.y));
  return contour_to_poly(pts);
}

Nested<TV2> square_poly(T rad) {
  return square_poly(vec(-rad,-rad), vec(rad,rad));
}

Nested<TV2> all_poly(void) {
  T x = 10.0;
  return square_poly(vec(-x, -x), vec(x, x));
}

Nested<TV2> none_poly(void) {
  Nested<TV2> res;
  return res;
}

Nested<TV2> circle_poly(T rad, int n) {
  Array<Vector<real, 2>> pts;
  for (int i = n-1; i >= 0; i--) {
    T a = (2 * M_PI * i) / n;
    // printf("A = %f\n", a);
    pts.append(vec(rad * sin(a),rad * cos(a)));
  }
  return contour_to_poly(pts);
}

Nested<TV2> star_poly(T rad_min, T rad_max, int n) {
  Array<Vector<real, 2>> pts;
  for (int i = n-1; i >= 0; i--) {
    T a_min = (2 * M_PI * (2 * i)) / (2 * n);
    pts.append(vec(rad_min * sin(a_min),rad_min * cos(a_min)));
    T a_max = (2 * M_PI * ((2 * i)+1)) / (2 * n);
    // printf("A = [%f, %f]\n", a_min, a_max);
    pts.append(vec(rad_max * sin(a_max),rad_max * cos(a_max)));
  }
  return contour_to_poly(pts);
}

struct Meshy {
public:
  Array<TV> points;
  Array<int> indices;
};

#include <GLUT/glut.h> 

typedef void (*GluTessCallbackType)();

static void begin_callback(GLenum type, void* mesh) {
}

static void edge_flag_callback(GLboolean flag, void* mesh) {
}

static void vertex_callback (const int index, void* mesh) {
  reinterpret_cast<Meshy*>(mesh)->indices.append(index);
}

static void end_callback (void* mesh) {
}

static void error_callback (GLenum errno, void* mesh) {
}

static int combine_callback (GLdouble coords[3], unsigned int vertexData[4], 
                             GLfloat weight[4], unsigned int* outData, void* mesh_) {
  auto mesh = reinterpret_cast<Meshy*>(mesh_);
  mesh->points.append(vec(coords[0], coords[1], coords[2]));
  return mesh->points.size() - 1;
}

static GLUtesselator* triangulator_new (void) {
  auto tess = gluNewTess();
  gluTessCallback(tess, GLU_TESS_BEGIN_DATA,
                  reinterpret_cast<GluTessCallbackType>(begin_callback));
  gluTessCallback(tess, GLU_TESS_EDGE_FLAG_DATA,
                  reinterpret_cast<GluTessCallbackType>(edge_flag_callback));
  gluTessCallback(tess, GLU_TESS_VERTEX_DATA,
                  reinterpret_cast<GluTessCallbackType>(vertex_callback));
  gluTessCallback(tess, GLU_TESS_END_DATA,
                  reinterpret_cast<GluTessCallbackType>(end_callback));
  gluTessCallback(tess, GLU_TESS_COMBINE_DATA,
                  reinterpret_cast<GluTessCallbackType>(combine_callback));
  gluTessCallback(tess, GLU_TESS_ERROR_DATA,
                  reinterpret_cast<GluTessCallbackType>(error_callback));
  return tess;
}

Mesh triangulate (Nested<TV> poly) { 
  auto tess = triangulator_new();
  Meshy mesh;
  gluTessBeginPolygon(tess, reinterpret_cast<void*>(&mesh));
  for (auto c : poly) {
    gluTessBeginContour(tess);
    for (auto p : c) {
      mesh.points.append(p);
      gluTessVertex(tess, reinterpret_cast<double*>(&p),
                    reinterpret_cast<void*>(mesh.points.size()-1));
    }
    gluTessEndContour(tess);
  }
  gluTessEndPolygon(tess);
  gluDeleteTess(tess);
  Array<IV> faces;
  for (int i = 0; i < mesh.indices.size()/3; i++) {
    faces.append(vec(mesh.indices[i*3], mesh.indices[i*3 + 2], mesh.indices[i*3 + 1]));
  }
  return fab_mesh(faces, mesh.points);
}

/*
Tuple<Ref<TriangleSoup>,Array<TV>> triangulate(Array<TV2> poly) {
  Array<IV2> edges;
  int n = poly.size();
  for (int i = 0; i < n; i++)
    edges.append(vec(i, (i + 1)%n));
  RawArray<const TV2> raw_poly(poly);
  RawArray<const IV2> raw_edges(edges);
  auto tri_topo = delaunay_points(raw_poly, raw_edges);
  return triangle_soup(tri_topo);
}
*/

TV mul(Matrix<T,4> m, TV pt) {
  return m.homogeneous_times(pt);
}

TV2 mul(Matrix<T,4> m, TV2 pt) {
  auto res = m.homogeneous_times(vec(pt.x, pt.y, 1.0));
  return vec(res.x, res.y);
}

Mesh xxx(Matrix<T,4> m, Mesh soup) {
  return tuple(soup.x, xxx(m, soup.y));
}

Mesh cone_mesh(T len, Array<TV2> poly) {
  Array<TV> bot;
  Array<TV> vh;

  T zmin = - len / 2.0;
  T zmax =   len / 2.0;
  auto c = vec(0.0, 0.0);
  int n_boundary = poly.size();
  for (auto elt : poly) {
    bot.append(vec(elt.x, elt.y, zmin));
    c = c + elt;
  }
  c = c / (double)n_boundary;
  auto lid = triangulate(contour_to_poly(bot));

  vh = lid.y;
  vh.append(vec(c.x, c.y, zmax));
    
  Array<IV> faces;
  for (auto face : lid.x->elements)
    faces.append(face);
                 
  for (int i = 0; i < n_boundary; i++) {
    int ni = (i + 1)%n_boundary;
    faces.append(vec(i, n_boundary, ni));
  }
  return fab_mesh(faces, vh);
}

Mesh cone_mesh(T len, Nested<TV2> contours) {
  ensure(contours.size() == 1, "CONE ONLY WORKS ON SINGLE CONTOUR POLYGONS");
  return cone_mesh(len, poly_to_contour(contours));
}

Mesh taper_mesh(T len, T r0, T r1, Array<TV2> poly) {
  Array<TV> top, bot;
  Array<TV> vh;

  T zmin = - len / 2.0;
  T zmax =   len / 2.0;
  int n_boundary = poly.size();
  // printf("N_BOUNDARY = %d\n", n_boundary);
  for (auto elt : poly)
    top.append(vec(r0 * elt.x, r0 * elt.y, zmin));
  auto top_mesh = triangulate(contour_to_poly(top));
  int n_mesh = top_mesh.y.size();
  // printf("N_MESH = %d\n", n_mesh);

  for (auto elt : poly) 
    bot.append(vec(r1 * elt.x, r1 * elt.y, zmax));
  auto bot_mesh = triangulate(contour_to_poly(bot));

  // auto lids = concat_soups(top_mesh, invert_soup(bot_mesh));
  auto lids = concat_soups(top_mesh, invert_soup(bot_mesh));

  // for (auto pt : vh)
  //   printf("PT [%f,%f,%f]\n", pt.x, pt.y, pt.z);

  Array<IV> faces;
  for (auto face : lids.x->elements)
    faces.append(face);
                 
  for (int i = 0; i < n_boundary; i++) {
    int ni = (i + 1)%n_boundary;
    // faces.append(vec(i, ni, n + ni));
    // printf("FACE [%d,%d,%d]\n", i,           ni,          n_mesh + ni);
    // faces.append(vec(n + ni, n + i, i));
    // printf("FACE [%d,%d,%d]\n", n_mesh + ni, n_mesh + i,  i);
    faces.append(vec(i,           ni,         n_mesh + ni));
    faces.append(vec(n_mesh + ni, n_mesh + i, i));
  }
  return fab_mesh(faces, lids.y);
}

Mesh taper_mesh(T len, T r0, T r1, Nested<TV2> contours) {
  if (contours.size() == 1) {
    if (true || contours.size() == 1)
      return taper_mesh(len, r0, r1, poly_to_contour(contours, 0));
    else
      return taper_mesh(len, r0, r1, poly_to_contour(contours, 1));
  } else {
    auto res = none_soup();
    for (int i = 0; i < contours.size(); i++) {
      auto contour = poly_to_contour(contours, i);
      if (!is_clockwise(contour)) {
        printf("OUTER %d\n", i);
        auto mesh = taper_mesh(len, r0, r1, contour);
        // pretty_print_soup(mesh);
        res = add(res, mesh);
      }
    }
    for (int i = 0; i < contours.size(); i++) {
      auto contour = poly_to_contour(contours, i);
      if (is_clockwise(contour)) {
        printf("INNER %d\n", i);  
        auto mesh = xxx(scale_matrix(vec(1.0,1.0,2.0)), taper_mesh(len, r0, r1, contour));
        // pretty_print_soup(mesh);
        res = add(res, mesh);
      }
    }
    return res;
  }
}

Mesh extrude(T len, Nested<TV2> poly) {
  return taper_mesh(len, 1.0, 1.0, poly);
}

Mesh revolve(int n, Array<TV2> poly) {
  Array<TV> pts;

  for (int i = 0; i < n; i++) {
    T a = (2 * M_PI * i) / n;
    for (auto elt : poly) {
      auto rad = elt.x;
      pts.append(vec(rad * sin(a), elt.y, rad * cos(a)));
    }
  }

  int m = poly.size();

  Array<IV> faces;

  for (int j = 0; j < n; j++) {
    int nj = (j + 1)%n; 
    for (int i = 0; i < m; i++) {
      int ni = (i + 1)%m;
      faces.append(vec((j * m) + i,   (nj * m) + ni, (j * m) + ni));
      faces.append(vec((nj * m) + ni, (j * m) + i, (nj * m) + i));
    }
  }
  return fab_mesh(faces, pts);
}

Mesh revolve(int n, Nested<TV2> contours) {
  if (contours.size() == 1)
    return revolve(n, poly_to_contour(contours, 0));
  else {
    auto res = none_soup();
    for (int i = 0; i < contours.size(); i++) {
      auto contour = poly_to_contour(contours, i);
      if (!is_clockwise(contour))
        res = add(res, revolve(n, contour));
    }
    for (int i = 0; i < contours.size(); i++) {
      auto contour = poly_to_contour(contours, i);
      if (is_clockwise(contour)) 
        res = add(res, revolve(n, contour));
    }
    return res;
    // return revolve(n, poly_to_contour(contours));
  }
}

Nested<TV2> fat_square_edge(int n, T rad, TV2 from, TV2 to) {
  auto q0 = vec(min(from.x, to.x), min(from.y, to.y));
  auto q1 = vec(max(from.x, to.x), max(from.y, to.y));
  auto v0 = q0 - vec(rad, rad);
  auto v1 = q1 + vec(rad, rad);
  if (from.x == to.x || from.y == to.y) {
    // printf("SQUARE FAT\n");
    return square_poly(v0, v1);
  } else {
    // printf("GEN FAT\n");
    Array< TV2 > points;
    points.append(v0);
    points.append(vec(v0.x, v1.y));
    points.append(v1);
    points.append(vec(v1.x, v0.y));
    return contour_to_poly(points);
  }
}

Nested<TV2> thicken(int n, T rad, Nested<TV2> line) {
  Nested<TV2> res = none_poly();
  int j = 0;
  // printf("THICKEN START\n");
  for (auto contour : line) {
    // printf("CONTOUR %d/%d\n", j + 1, line.size());
    for (int i = 1; i < contour.size(); i++) {
      // printf("EDGE %d/%d\n", i, contour.size());
      auto edge = fat_square_edge(n, rad, contour[i-1], contour[i]);
      // printf("UNIONING RES\n");
      // print_polyline2(res);
      // printf("TO EDGE\n");
      // print_polyline2(edge);
      res = add(res, edge);
    }
    j += 1;
  }
  // printf("THICKEN DONE\n");
  return res;
}

Mesh fat_edge(int n, T rad, TV from, TV to) {
  auto v = to - from;
  auto res = extrude(magnitude(v), circle_poly(rad, n));
  auto c = (to + from) * 0.5;
  return xxx(translation_matrix(c), xxx(rotation_matrix(vec(0.0, 0.0, 2*rad), v), res));
}

Nested<TV2> fat_edge(int n, T rad, TV2 from, TV2 to) {
  auto v = to - from;
  auto res = square_poly(vec(-rad, -magnitude(v)/2), vec(rad, magnitude(v)/2));
  auto c = (to + from) * 0.5;
  return xxx(translation_matrix(c), xxx(rotation_matrix(vec(0.0, magnitude(v)), v), res));
}

Nested<TV2> fat_dot(int n, T rad, TV2 pt) {
  auto circ = circle_poly(rad, n);
  auto res = xxx(translation_matrix(pt), circ);
  return res;
}

Nested<TV2> offset_poly(int n, T rad, Nested<TV2> poly) {
  Nested<TV2> res = poly;
  for (auto contour : poly) {
    for (int i = 0; i < contour.size(); i++) {
      res = add(res, add(fat_dot(n, rad, contour[i]), fat_edge(n, rad, contour[i], contour[(i+1)%contour.size()])));
    }
  }
  return res;
}

Nested<TV2> offset_polyline(int n, T rad, Nested<TV2> polyline) {
  Nested<TV2> res = none_poly();
  for (auto line : polyline) {
    for (auto pt : line) 
      res = add(res, fat_dot(n, rad, pt));
    for (int i = 0; i < (line.size()-1); i++) 
      res = add(res, fat_edge(n, rad, line[i], line[(i+1)%line.size()]));
  }
  return res;
}

Mesh thicken(int n, T rad, Mesh mesh) {
  Mesh res = const_soup(sphere_mesh(n, mesh.y[0], rad));
  // printf("TH %f,%f,%f\n", line[0].x, line[0].y, line[0].z);

  for (int i = 1; i < mesh.y.size(); i++) {
    // printf("TH %f,%f,%f\n", line[i].x, line[i].y, line[i].z);
    res = add(res, const_soup(sphere_mesh(n, mesh.y[i], rad)));
  }
  auto edges = mesh.x->triangle_edges();
  for (auto edge : edges)
    printf("EDGE %d to %d\n", edge.x, edge.y);
  for (auto edge : edges) {
    res = add(res, fat_edge(16, rad, mesh.y[edge.x], mesh.y[edge.y]));
  }
  
  return res;
}

Mesh thicken(int n, T rad, Nested<TV> polyline) {
  Mesh res = none_soup();
  for (auto line : polyline) {
    // printf("TH %f,%f,%f\n", line[0].x, line[0].y, line[0].z);
    for (int i = 0; i < line.size(); i++) {
      // printf("TH %f,%f,%f\n", line[i].x, line[i].y, line[i].z);
      res = add(res, const_soup(sphere_mesh(n, line[i], rad)));
    }
    for (int i = 1; i < line.size(); i++) {
      res = add(res, fat_edge(12, rad, line[i-1], line[i]));
    }
  }
  
  return res;
}

Mesh offset_mesh(int n, T rad, Mesh mesh) {
  /*
  Array<TV> pos(mesh.y);
  auto topo = new_<TriangleTopology>(mesh.x->elements);
  auto new_mesh = rough_offset_mesh(topo, RawField<const TV,VertexId>(pos), rad);
  auto new_soup = new_mesh.x->face_soup().x;
  Array<TV> new_points;
  for (auto point : new_mesh.y.flat)
    new_points.append(point);
  return tuple(const_soup(new_soup), new_points);
  */
  return mesh;
}

/*
Mesh offset_mesh(int n, T rad, Mesh mesh) {
  Mesh res = none_soup();

  // printf("%d POINTS\n", mesh.y.size());
  // for (auto pt : mesh.y) {
  //   printf("TH %f,%f,%f\n", pt.x, pt.y, pt.z);
  //   res = add(res, const_soup(sphere_mesh(n, pt, rad)));
  // }
  // res = simplify_mesh(res);
  auto edges = mesh.x->triangle_edges();
  int i = 0;
  for (auto edge : edges) {
    printf("EDGE %d/%d %d to %d\n", i, edges.size(), edge.x, edge.y);
    res = add(res, fat_edge(8, rad, mesh.y[edge.x], mesh.y[edge.y]));
    res = simplify_mesh(res);
    i += 1;
  }
  return res;
}
*/

static std::string letter_codes[256];
static Nested<TV2> letter_outlines[256];
static double x_specs[256];
static double y_specs[256];

void fab_stroke (std::string &s, int i, Array<TV2>& stroke) {
  auto cx = s[i];
  if (cx != ':' && cx != ';') {
    stroke.append(vec(x_specs[(int)cx],y_specs[(int)s[i + 1]]));
    fab_stroke(s, i + 2, stroke);
  }
}

Nested<TV2> fab_letter_outline (char c, std::string &s) {
  size_t k = 0;
  Nested<TV2,false> segments;
  while (k < s.size()) {
    Array<TV2> stroke;
    fab_stroke(s, k, stroke);
    k = k + stroke.size() * 2 + 1;
    segments.append(stroke);
  }
  segments.freeze();
  return segments;
}
      
Nested<TV2> stroke_char (char letter) {
  return letter_outlines[(int)letter];
}

Nested<TV2> stroke_text (std::string txt) {
  Nested<TV2,false> res;
  auto smat = scale_matrix(vec(0.8, 0.8, 1.0));
  for (size_t i = 0; i < txt.size(); i++) {
    auto dx = (i + 0.5) - (txt.size() * 0.5);
    auto tmat = translation_matrix(vec(dx,0.0,0.0));
    res.extend(xxx(tmat * smat, letter_outlines[(int)txt[i]]));
  }
  return res;
}


static bool is_init_letters = false;

static void init_letters ( void ) {
  if (!is_init_letters) {
    for (int i = 0; i < 256; i++)
      letter_codes[i] = "";
    letter_codes[(int)'0'] = "LBLTRTRBLB;CTCB;";
    letter_codes[(int)'1'] = "CTCB;";
    letter_codes[(int)'2'] = "LTRTRCLCLBRB;"; // straight up
    letter_codes[(int)'3'] = "LTRTRBLB;LCRC;";
    letter_codes[(int)'4'] = "LTLCRC;RTRB;";
    letter_codes[(int)'5'] = "LBRBRCLCLTRT;";
    letter_codes[(int)'6'] = "LCRCRBLBLTRT;";
    letter_codes[(int)'7'] = "LTRTRB;";
    letter_codes[(int)'8'] = "LBLTRTRBLB;LCRC;";
    letter_codes[(int)'9'] = "LBRBRTLTLCRC;";
    letter_codes[(int)'A'] = "LBLTRTRB;LCRC;";
    letter_codes[(int)'B'] = "LBLTRTRBLB;LCRC;CTCC;";
    letter_codes[(int)'C'] = "RBLBLTRT;";
    letter_codes[(int)'D'] = "LBLTRTRBLB;LCCCCT;"; // extra square up top
    letter_codes[(int)'E'] = "RTLTLBRB;LCRC;";
    letter_codes[(int)'F'] = "LBLTRT;LCRC;";
    letter_codes[(int)'G'] = "RTLTLBRBRCCC;";
    letter_codes[(int)'H'] = "LTLB;LCRC;RTRB;";
    letter_codes[(int)'I'] = "LTRT;CTCB;LBRB;";
    letter_codes[(int)'J'] = "LCLBRBRT;";
    letter_codes[(int)'K'] = "LTLB;LCRC;RTRB;CTCC;"; // like H with extra top vert
    letter_codes[(int)'L'] = "LTLBRB;";
    letter_codes[(int)'M'] = "LBLTRTRB;CTCB;";
    letter_codes[(int)'N'] = "LBLTRTRB;";
    letter_codes[(int)'O'] = "LBLTRTRBLB;";
    letter_codes[(int)'P'] = "LBLBLTRTRCLC;";
    letter_codes[(int)'Q'] = "LBLTRTRBLB;CBCC;";
    letter_codes[(int)'R'] = "LBLTRTRB;LCRC;CTCC;"; // like A with extra top vert
    letter_codes[(int)'S'] = "LBRBRCLCLTRT;";
    letter_codes[(int)'T'] = "LTRT;CTCB;";
    letter_codes[(int)'U'] = "LTLBRBRT;";
    letter_codes[(int)'V'] = "LTLBRBRT;LCCCCB;"; 
    letter_codes[(int)'W'] = "LTLBRBRT;CTCB;";
    letter_codes[(int)'X'] = "LTLB;LCRC;RTRB;CBCT;"; // like H with extra vert
    letter_codes[(int)'Y'] = "LTLCRC;RTRBLB;"; // y WITH TAIL
    letter_codes[(int)'Z'] = "LTRTRCLCLBRB;CTCB;";
    x_specs[(int)'L'] = -0.5;
    x_specs[(int)'C'] = 0.0;
    x_specs[(int)'R'] = 0.5;
    y_specs[(int)'B'] = -0.5;
    y_specs[(int)'C'] = 0.0;
    y_specs[(int)'T'] = 0.5;
    for (int i = 0; i < 256; i++) {
      auto code = letter_codes[i];
      if (code != "") {
        letter_outlines[i] = invert_poly(fab_letter_outline((char)i, code));
      }
    }
    is_init_letters = true;
  }
}

void init_cad ( void ) {
  init_letters();
}

// static void ensure (bool val, std::string msg, std::string arg) {
//   if (!val) error(msg, arg);
// }
