#include "cad.h"
#include "ZippedView.h"
#include <geode/array/sort.h>
#include <geode/mesh/improve_mesh.h>
// #include <geode/mesh/decimate.h>
#include <fstream>

bool has_dups(const RawArray<const Vec3i> raw_elements) {
  Array<Vec3i> elements = raw_elements.copy();

  std::sort(elements.begin(), elements.end(), [](const Vec3i lhs, const Vec3i rhs) {
    return lex_less(lhs.sorted(), rhs.sorted());
  });

  for(const int i : range(elements.size()-1)) {
    if(elements[i].sorted() == elements[i+1].sorted()) {
      return true;
    }
  }
  return false;
}


double rndd () {
  return (double)((double)rand() / (double)RAND_MAX);
}

double rndd (double mn, double mx) {
  return rndd() * (mx-mn) + mn;
}

Vec rnd_vec_of (const Vec &lo, const Vec &hi) { 
  return vec(rndd(lo.x,hi.x), rndd(lo.y,hi.y), rndd(lo.z,hi.z));
}

Vec rnd_vec_of (const flo_t lo, const flo_t hi) { 
  return vec(rndd(lo,hi), rndd(lo,hi), rndd(lo,hi));
}

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

Mesh fab_mesh (Array<IV3> faces, Array<TV3> points) {
  // return tuple(new_<const TriangleSoup>(faces),points);
  return Mesh(new_<const TriangleSoup>(faces),points);
}
Mesh fab_mesh (Ref<const TriangleSoup> soup, Array<TV3> points) {
  return Mesh(soup,points);
}

typedef real T;
typedef Vector<T,2> TV2;
typedef Vector<T,3> TV;
typedef Vector<int,3> IV3;
typedef Vector<int,2> IV2;

Matrix<T,4> to_matrix44(Matrix<T,3> M) {
  return Matrix<T,4>(M.x[0][0],M.x[1][0],M.x[2][0],0,M.x[0][1],M.x[1][1],M.x[2][1],0,M.x[0][2],M.x[1][2],M.x[2][2],0,0,0,0,1);
}
Matrix<T,4> rotation_matrix(TV angles) {
  auto m = Rotation<Vector<T,3> >::from_euler_angles(angles).matrix();
  return to_matrix44(m);
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
static Mesh concat_meshes(Mesh mesh0, Mesh mesh1) {
  Array<IV3> triangles;
  Array<TV3> vertices;

  for(const auto& s : { mesh0, mesh1 }) {
    // 0th vertex in each soup will be first vertex after all previous ones
    // Save as a vector of ints to we can directly add it to each triangle
    const auto index_offset = IV3::repeat(vertices.size());
    vertices.extend(s.points);
    for(const auto& t : s.soup->elements) {
      triangles.append(t+index_offset); // Add new triangle with indices mapped to combined vertices
    }
  }

  return fab_mesh(triangles,vertices);
}

Ref<const TriangleSoup> const_soup(Ref<TriangleSoup> val) {
  return new_<const TriangleSoup>(val->elements);
}

Mesh const_mesh(Tuple<Ref<TriangleSoup>, Array<TV3>> val) {
  return Mesh(const_soup(val.x), val.y);
}

void report_simplify_mesh(Mesh mesh) {
  bool is_printed_headline = false;
  bool is_found_problem = false;
  int fi = 0;
  for (auto face : mesh.soup->elements) {
    Triangle<TV3> tri(mesh.points[face[0]],
                      mesh.points[face[1]],
                      mesh.points[face[2]]);
    if (tri.area() == 0.0) {
      is_found_problem = true;
      if (!is_printed_headline) {
        is_printed_headline = true;
        printf("START REPORTING MESH ISSUES\n");
      }
      printf("ZERO AREA TRI %d: ", fi);
      printf("[%f,%f,%f] [%f,%f,%f] [%f,%f,%f]: ",
             tri.x0.x, tri.x0.y, tri.x0.z, tri.x1.x, tri.x1.y, tri.x1.z, tri.x2.x, tri.x2.y, tri.x2.z);
      if (tri.x0 != tri.x1 && tri.x0 != tri.x2 && tri.x1 != tri.x2) 
        printf("CUP\n");
      else
        printf("SPLINTER\n");
    }
    fi += 1;
  }
  Hashtable<TV3,int> same_vertices;
  for (int i = 0; i < mesh.points.size(); i++) {
    auto p = mesh.points[i];
    auto j = same_vertices.get_or_insert(p, i);
    if (i != j) {
      is_found_problem = true;
      if (!is_printed_headline) {
        is_printed_headline = true;
        printf("START REPORTING MESH ISSUES\n");
      }
      printf("REDUNDANT POINT [%f,%f,%f] %d -> %d \n", p.x, p.y, p.z, i, j);
    }
  }
  if (is_found_problem)
    printf("END REPORTING MESH ISSUES\n");
}

Mesh dither_mesh(Mesh mesh, double delta) {
  Array<TV3> points;
  for (auto point : mesh.points)
    points.append(point + vec(rndd(-delta, delta), rndd(-delta, delta), rndd(-delta, delta)));
  return fab_mesh(mesh.soup, points);
}

Mesh unify_mesh_vertices (Mesh mesh) {
  // printf("GCING\n");
  printf("START COMPRESS VERTICES\n");
  // pretty_print_mesh(mesh);
  // printf("---\n");
  Hashtable<TV3,int> same_vertices;
  Hashtable<int,int> vertex_map;
  for (int i = 0; i < mesh.points.size(); i++) {
    auto p = mesh.points[i];
    auto j = same_vertices.get_or_insert(p, i);
    if (i != j)
      printf("MAPPING [%f,%f,%f] %d -> %d [%f,%f,%f] \n", p.x, p.y, p.z, i, j, mesh.points[j].x, mesh.points[j].y, mesh.points[j].z);
    vertex_map[i] = j;
  }
  Array<IV3> new_faces;
  for (auto face : mesh.soup->elements) 
    new_faces.append(vec(vertex_map[face.x], vertex_map[face.y], vertex_map[face.z]));
  printf("END COMPRESS VERTICES\n");
  return fab_mesh(new_faces, mesh.points);
}

/*
Mesh quick_simplify_mesh(Mesh mesh) {
  printf("SIMPLIFYING\n");
  Array<IV3> faces;
  for (auto face : mesh.soup->elements)
    faces.append(face);
  auto points = mesh.points.copy();
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
      auto new_face = vec(mapping[face.x], mapping[face.y], mapping[face.z]);
      if (new_face.x != faces[i].x || new_face.y != faces[i].y || new_face.z != faces[i].z) {
        printf("MAPPING [%d,%d,%d] => [%d,%d,%d]\n",
               faces[i].x, faces[i].y, faces[i].z, new_face.x, new_face.y, new_face.z);
        faces[i] = new_face;
      }
    }
  }
  Array<IV3> new_faces;
  for (auto face : faces) {
    if (face.x != face.y && face.x != face.z && face.y != face.z) {
      new_faces.append(face);
    } else
      printf("REMOVED FACE %d,%d,%d\n", face.x, face.y, face.z);
  }
  printf("%d FACES NOW %d REMOVED %d FACES\n",
         faces.size(), new_faces.size(), faces.size() - new_faces.size() );
  return gc_mesh(fab_mesh(new_faces, points));
}
*/

int max_segment (Segment<TV3>& s0, Segment<TV3>& s1, Segment<TV3>& s2) {
  // printf("S0 %f S1 %f S2 %f\n", s0.length(), s1.length(), s2.length());
  auto ml = max(s0.length(), max(s1.length(), s2.length()));
  if (ml == s0.length())
    return 0;
  else if (ml == s1.length())
    return 1;
  else if (ml == s2.length())
    return 2;
  else  {
    error("COULD NOT FIND MAX SEGMENT");
    return -1;
  }
}

VertexId common_vertex (VertexId e0v0, VertexId e0v1, VertexId e1v0, VertexId e1v1) {
  // FIND VERTEX COMMON TO BOTH EDGES
  // TODO: PROBABLY IMPLICIT IN REPRESENTATION
  if (e0v0 == e1v0 || e0v0 == e1v1) return e0v0;
  else if (e0v1 == e1v0 || e0v1 == e1v1) return e0v1;
  else {
    error("COULD NOT FIND COMMON VERTEX"); return (VertexId)0;
  }
}

Mesh gc_mesh(Mesh mesh) {
  // printf("GCING\n");
  // pretty_print_mesh(mesh);
  Array<bool> is_points;
  for (int i = 0; i < mesh.points.size(); i++)
    is_points.append(false);
  Array<IV3> new_faces;
  for (auto face : mesh.soup->elements) 
    is_points[face.x] = is_points[face.y] = is_points[face.z] = true;
  int delta = 0;
  Array<int> mapping;
  Array<TV3> new_points;
  for (int i = 0; i < mesh.points.size(); i++) {
    mapping.append(i - delta);
    if (is_points[i]) {
      // printf("MAPPING %d TO %d\n", i, i - delta);
      new_points.append(mesh.points[i]);
    } else {
      // printf("REMOVING %d DELTA %d\n", i, delta);
      delta += 1;
    }
  }
  for (auto face : mesh.soup->elements) 
    new_faces.append(vec(mapping[face.x], mapping[face.y], mapping[face.z]));
  // printf("%d POINTS NOW %d REMOVED %d POINTS DELTA %d\n",
  //        mesh.points.size(), new_points.size(), mesh.points.size() - new_points.size(), delta);
  return fab_mesh(new_faces, new_points);
}

Nested<TV3> nested2_to_nested3 (Nested<TV2> contours2) {
  Nested<TV3, false> contours3;
  for (auto c : contours2) {
    Array<TV3> contour3;
    for (auto p : c) 
      contour3.append(vec(p.x, p.y, 0.0));
    contours3.append(contour3);
  }
  contours3.freeze();
  return contours3;
}

Mesh topo_to_mesh (Ref<MutableTriangleTopology> topo, const Field<TV3,VertexId>& field) {
  auto updates = topo->collect_garbage();
  int i = 0, tot = 0;
  for (auto update : updates.x) {
    if (update >= 0) tot += 1;
    // printf("  %d -> %d\n", i, update);
    i += 1;
  }
  // printf("AFTER GC %d UPDATES %d - %d\n", field.size(), updates.x.size(), tot);
  Array<TV3> new_points(tot);
  i = 0;
  for (auto update : updates.x) {
    if (update >= 0) {
      new_points[update] = field[VertexId(i)];
      // printf("MAPPING %d -> %d\n", i, update);
    }
    i += 1;
  }
  auto new_soup = topo->face_soup().x;
  // printf("SIMPLIFYING: BEFORE %d,%d AFTER %d,%d\n",
  //        pos.size(), mesh.x->elements.size(), new_points.size(), new_soup->elements.size());
  // write_mesh("tst0.stl", mesh.soup, mesh.points);
  auto new_mesh = Mesh(const_soup(new_soup), new_points);
  auto gcd_mesh = gc_mesh(new_mesh);
  // auto gcd_mesh = final_mesh;
  report_simplify_mesh(gcd_mesh);
  return gcd_mesh;
}

bool report_unsafe_collapse (Ref<MutableTriangleTopology> topo, HalfedgeId h) {
  const auto o = topo->reverse(h);
  const auto v0 = topo->src(h),
             v1 = topo->dst(h);

  // If v0 and v1 are on different boundaries, we can't do this
  if (topo->is_boundary(v0) && topo->is_boundary(v1) &&
      !topo->is_boundary(h) && !topo->is_boundary(o)) {
    printf("If v0 and v1 are on different boundaries, we can't do this\n");
    return false;
  }

  // Can't snip off an isolated vl or vr
  if ((topo->is_boundary(topo->reverse(topo->next(h))) && topo->is_boundary(topo->reverse(topo->prev(h)))) ||
      (topo->is_boundary(topo->reverse(topo->next(o))) && topo->is_boundary(topo->reverse(topo->prev(o))))) {
    printf("Can't snip off an isolated vl or vr\n");
    return false;
  }

  // Look up left and right vertices
  const auto vl = topo->is_boundary(h) ? VertexId() : topo->dst(topo->next(h)),
             vr = topo->is_boundary(o) ? VertexId() : topo->dst(topo->next(o));

  // This only happens in temporarily invalid situations, such as if
  // split_along_edge is called and not cleaned up.  No good can come of it.
  if (vl==vr) {
    printf("INVALID EDGE\n");
    return false;
  }

  // One-rings of v0 and v1 cannot intersect, otherwise we'll collapse a
  // triangle-shaped tunnel
  Hashtable<VertexId> covered;
  for (auto oh: topo->outgoing(v0)) {
    auto v = topo->dst(oh);
    if (v != vl && v != vr)
      covered.set(v);
  }
  for (auto oh: topo->outgoing(v1)) {
    auto v = topo->dst(oh);
    if (covered.contains(v)) {
      printf("ONE-RINGS of V0 and V1 INTERSECT\n");
      return false;
    }
  }

  return true;
}

Mesh quick_cleanup_mesh(Mesh mesh) {
  printf("--- START CLEANUP\n");
  // pretty_print_mesh(mesh);
  // write_mesh("tst1.stl", mesh.soup, mesh.points);
  // printf("---\n");
  // report_simplify_mesh(mesh);
  // Array<TV3> pos(mesh.points.copy());
  auto topo = new_<MutableTriangleTopology>(mesh.soup);
  printf("ADDING FIELD\n");
  FieldId<TV3,VertexId> pos_id = topo->add_field(Field<TV,VertexId>(mesh.points.copy()), vertex_position_id);
  auto &field = topo->field(pos_id);

  // topo->dump_internals();
  // pretty_print_mesh(topo_to_mesh(topo, field));
  // std::vector<FaceId> split_faces;
  bool is_changed = false;
  printf("CLEANUP ZERO GEOMETRY\n");
  do {
     do {
      is_changed = false;
      for (auto face : topo->faces()) {
        auto tri = topo->triangle(field, face);
        if (tri.area() == 0.0 && tri.x0 != tri.x1 && tri.x0 != tri.x2 && tri.x1 != tri.x2) {
          auto e0  = topo->halfedge(face, 0);
          auto e1  = topo->halfedge(face, 1);
          auto e2  = topo->halfedge(face, 2);
          RawField<TV3,VertexId> rawfield(field);
          auto s0  = topo->segment(rawfield, e0);
          auto s1  = topo->segment(rawfield, e1);
          auto s2  = topo->segment(rawfield, e2);
          // auto tri = topo->triangle(field, face);
          // printf("ZERO AREA TRI %d [%f,%f,%f] [%f,%f,%f] [%f,%f,%f]\n",
          //        (int)face, tri.x0.x, tri.x0.y, tri.x0.z, tri.x1.x, tri.x1.y, tri.x1.z, tri.x2.x, tri.x2.y, tri.x2.z);
          // printf("  E0 %d S0 [%f,%f,%f] [%f,%f,%f]\n", (int)e0, s0.x0.x, s0.x0.y, s0.x0.z, s0.x1.x, s0.x1.y, s0.x1.z);
          // printf("  E1 %d S1 [%f,%f,%f] [%f,%f,%f]\n", (int)e1, s1.x0.x, s1.x0.y, s1.x0.z, s1.x1.x, s1.x1.y, s1.x1.z);
          // printf("  E2 %d S2 [%f,%f,%f] [%f,%f,%f]\n", (int)e2, s2.x0.x, s2.x0.y, s2.x0.z, s2.x1.x, s2.x1.y, s2.x1.z);
          int msi  = max_segment(s0, s1, s2);
          int o1i  = (msi + 1)%3;
          int o2i  = (msi + 2)%3;
          auto oe1 = topo->halfedge(face,o1i);
          auto oe2 = topo->halfedge(face,o2i);
          auto cvi = common_vertex(topo->src(oe1), topo->dst(oe1), topo->src(oe2), topo->dst(oe2));
          // printf("BEFORE n-VERTS %d ALLOC VERTS %d\n", topo->n_vertices(), topo->allocated_vertices());
          auto nvi = topo->split_edge(topo->halfedge(face,msi));
          auto ov  = field[cvi];
          field[nvi] = ov;
          auto nv  = field[nvi];
          // field.append(ov);
          printf("SPLITTING EDGE %d (%d,%d,%d) ON FACE %d VERTEX %d -> %d [%f,%f,%f]->[%f,%f,%f] (%d) ALLOC %d\n",
                 (int)topo->halfedge(face,msi), int(e0), int(e1), int(e2),
                 (int)face, (int)cvi, (int)nvi, 
                 ov.x, ov.y, ov.z, nv.x, nv.y, nv.z, (int)field.size(), topo->allocated_vertices());
          topo->collect_garbage();
          // pretty_print_mesh(topo_to_mesh(topo, field));
          // topo->dump_internals();
          // printf("AFTER n-VERTS %d ALLOC VERTS %d\n", topo->n_vertices(), topo->allocated_vertices());
          break;
        }
      }
    } while (is_changed == true);
    // printf("COLLAPSE ZERO LENGTH EDGES\n");
    // do {
    is_changed = false;
    for (auto e : topo->halfedges()) {
      auto src = topo->src(e);
      auto dst = topo->dst(e);
      if (field[src] == field[dst]) {
        if (topo->is_collapse_safe(e)) {
          // printf("BEFORE n-VERTS %d ALLOC VERTS %d\n", topo->n_vertices(), topo->allocated_vertices());
          printf("COLLAPSING E%d V%d->V%d\n", (int)e, (int)src, (int)dst);
          topo->collapse(e);
          topo->collect_garbage();
          // pretty_print_mesh(topo_to_mesh(topo, field));
          // topo->dump_internals();
          // printf("AFTER n-VERTS %d ALLOC VERTS %d\n", topo->n_vertices(), topo->allocated_vertices());
          is_changed = true;
          break;
        } else {
          printf("UNABLE TO COLLAPSE %d %d->%d: ", (int)e, (int)src, (int)dst);
          report_unsafe_collapse(topo, e);
        }
      }
    }
  } while (is_changed == true);

  // printf("ERASING BOUNDARY EDGES\n");
  // do {
  //   is_changed = false;
  //   for (auto e : topo->boundary_edges()) {
  //     printf("ERASING BOUNDARY EDGE %d\n", (int)e);
  //     topo->erase(e, true);
  //     is_changed = true;
  //     break;
  //   }
  // } while (is_changed == true);

  printf("COLLECTING GARBAGE\n");
  // topo->dump_internals();
  auto new_mesh = topo_to_mesh(topo, field);
  // printf("DECIMATING MESH\n");
  // RawField<TV3,VertexId> rawfield(field);
  // do_decimate_inplace(topo, rawfield, 0.01, 0.01, 0, 0);
  // auto dec_mesh = decimate_mesh(new_mesh);
  // write_mesh("tst1.stl", new_mesh.soup, new_mesh.points);
  // auto dec_mesh = simplify_mesh(new_mesh);
  auto dec_mesh = new_mesh;
  auto final_mesh = dec_mesh;
  /*
  auto gc_new_mesh  = gc_mesh(dec_mesh);
  write_mesh("tst2.stl", gc_new_mesh.soup, gc_new_mesh.points);
  auto unified_mesh = unify_mesh_vertices(new_mesh);
  write_mesh("tst3.stl", unified_mesh.soup, unified_mesh.points);
  pretty_print_mesh(unified_mesh);
  auto final_mesh = unified_mesh;
  */
  report_simplify_mesh(final_mesh);
  printf("--- GC CLEANUP\n");
  auto gcd_mesh = gc_mesh(final_mesh);
  report_simplify_mesh(gcd_mesh);
  printf("--- END CLEANUP\n");
  return gcd_mesh;
  // return unified_mesh;
  // auto new_soup = topo->face_soup().x;
  // return gc_mesh(Mesh(const_soup(new_soup), pos));
}

Mesh simplify_mesh(Mesh mesh) {
  // printf("STARTING SIMPLIFICATION\n");
  // report_simplify_mesh(mesh);
  Array<TV3> pos(mesh.points);
  Field<TV3,VertexId> field(pos.copy());
  auto topo = new_<MutableTriangleTopology>();
  topo->add_vertices(mesh.points.size());
  for (auto face : mesh.soup->elements)
    topo->add_face(vec((VertexId)face.x, (VertexId)face.y, (VertexId)face.z));
  // ImproveOptions options(1.1,0.1,0.1);
  // ImproveOptions options(1.1,0.01,0.01);
  // ImproveOptions options(1.0000001,0.0000001,0.0000001);
  ImproveOptions options(1 + 1e-6,1e-6,1e-6);
  improve_mesh_inplace(topo, field, options);
  auto final_mesh = topo_to_mesh(topo, field);
  // printf("BEFORE GC %d\n", field.size());
  return final_mesh;
}

Mesh cleanup_mesh (Mesh mesh) {
  if (true) {
    return simplify_mesh(mesh);
  } else {
    return quick_cleanup_mesh(mesh);
  }
}

// invert triangle soup so normals point inwards
Mesh invert_mesh(Mesh mesh) {
  Array<IV3> triangles;

  for(const auto& t : mesh.soup->elements) {
    triangles.append(vec(t[0], t[2], t[1]));
  }

  return fab_mesh(triangles, mesh.points);
}

Array<TV2> invert_contour(Array<TV2> contour) {
  Array<TV2> res;
  for (int i = contour.size()-1; i >= 0; i--)
    res.append(contour[i]);
  return res;
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

Nested<TV2> mul_poly(Matrix<T,4> m, Nested<TV2> poly, bool is_invert) {
  auto res = mul(m, poly);
  return is_invert ? invert_poly(res) : res;
}

Array<TV2> mul_contour(Matrix<T,4> m, Array<TV2> contour, bool is_invert) {
  auto res = mul(m, contour);
  return is_invert ? invert_contour(res) : res;
}

Array<TV2> cleanup_contour(RawArray<TV2> contour) {
  // TODO: LIMITS
  return polygon_simplify(contour, 179.9, 0.0001);
}

Nested<TV2> cleanup_poly(Nested<TV2> poly) {
  Nested<TV2,false> res;
  for (auto contour : poly) 
    res.append(cleanup_contour(contour));
  res.freeze();
  return res;
}

Nested<TV2> maybe_cleanup_poly(Nested<TV2> poly) {
  if (true)
    return cleanup_poly(poly);
  else
    return poly;
}
Nested<TV2> union_add(Nested<TV2> c0, Nested<TV2> c1) {
  auto res = polygon_union(c0, c1);
  // printf("UNION POLY\n");
  return maybe_cleanup_poly(res);
}
Nested<TV2> union_all(Nested<TV2> c0) {
  auto res = polygon_union(c0);
  // printf("UNION POLY\n");
  return maybe_cleanup_poly(res);
}

Nested<TV2> intersection(Nested<TV2> c0, Nested<TV2> c1) {
  return maybe_cleanup_poly(polygon_intersection(c0, c1));
}

Nested<TV2> difference(Nested<TV2> c0, Nested<TV2> c1) {
  return maybe_cleanup_poly(polygon_union(c0, invert_poly(c1)));
}

void save_svg_header(std::fstream &fs) {
  fs << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
  fs << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
  fs << "<svg width=\"10cm\" height=\"10cm\" viewBox=\"0 0 100000 100000\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
}

void save_poly_star(std::string filename, std::string kind, Nested<TV2> poly) {
  std::fstream fs;
  fs.open (filename, std::fstream::out);
  save_svg_header(fs);
  for (auto c : poly) {
    fs << "<" << kind << " fill=\"none\" stroke=\"black\" stroke-width=\"1\" points=\"";
    for (auto p : c) {
      fs << " " << ((1000.0 * p.x) + 50000) << "," << ((-1000.0 * p.y) + 50000);
    }
    fs << "\"  />\n";
  }
  fs << "</svg>\n";
  fs.close();
}

void save_polygon(std::string filename, Nested<TV2> polygon) {
  save_poly_star(filename, "polygon", polygon);
}

void save_polyline(std::string filename, Nested<TV2> polyline) {
  save_poly_star(filename, "polyline", polyline);
}


Mesh maybe_cleanup_mesh (Mesh mesh, bool is_cleanup) {
  // report_cleanup_mesh(mesh);
  if (is_cleanup)
    return cleanup_mesh(mesh);
  else
    return mesh;
}

static void simplify_duplicate_faces(Array<Vec3i>& elements, Array<int>& depth_weights) {

  // Sort all elements so that duplicate faces will be adjacent
  auto zipped = zipped_view(elements, depth_weights);
  using Zipped = decltype(zipped);
  std::sort(zipped.begin(), zipped.end(), [](const Zipped::Val& lhs, const Zipped::Val& rhs) {
    return lex_less(std::get<0>(lhs).sorted(), std::get<0>(rhs).sorted());
  });

  // Scan over elements only keeping the first in groups of duplicates and updating depth_weights accordingly
  // Erase any elements where net depth_weights is 0
  if(!elements.empty()) {
    // Check if sorting f flips normal of triangle 
    const auto permutation_parity = [](const Vec3i f) {
      const int f0 = f.argmin();
      return (f[(f0+1) % 3] < f[(f0+2) % 3]);
    };

    int prev = 0; // Used as output iterator

    for(const int next : range(1,zipped.size())) {
      // Is this a new triangle?
      if(elements[prev].sorted() == elements[next].sorted()) {
        // Duplicate triangle. Need to check if normal is same or opposite of prev and update weight accordingly
        depth_weights[prev] += depth_weights[next] * ((permutation_parity(elements[prev]) == permutation_parity(elements[next])) ? 1 : -1);
      }
      else {
        // Found a new triangle
        // Unless previous triangle had a nonzero weight, we overwrite it
        if(depth_weights[prev] != 0) {
          prev += 1;
        }
        elements[prev] = elements[next];
        depth_weights[prev] = depth_weights[next];
      }
    }

    // If final unique triangle had zero net weight, erase it
    if(depth_weights[prev] == 0) {
      prev -= 1;
    }

    // Resize to match compacted output
    elements.resize(prev + 1);
    depth_weights.resize(prev + 1);
  }
}

// TODO: This function should be merged into default behavior for geode::split_soup
Tuple<Ref<const TriangleSoup>, Array<TV>> safe_split_soup(const TriangleSoup& faces, Array<const TV> X, const int depth) {
  // Perturbation allows us to ignore most degeneracies, but duplicates of the same face are still an issue
  // If faces are duplicated, crossing needs to account for total depth of all duplicates
  // If we aren't filtering by depth we leave duplicated faces. This ensures we don't cull duplicates and is safe sense depths aren't computed
  if(depth == all_depths) {
    return split_soup(faces, X, depth);
  }

  // Copy elements so we can find duplicates
  Array<Vec3i> elements = faces.elements.copy();
  auto depth_weights = Array<int>{elements.size()};
  depth_weights.fill(1);
  simplify_duplicate_faces(elements, depth_weights);

  assert(!has_dups(elements));
  return split_soup(new_<TriangleSoup>(elements), X, depth_weights, depth);
}

Mesh split_mesh (Mesh mesh, int depth, bool is_cleanup) {
  auto split = safe_split_soup(mesh.soup, mesh.points, depth);
  return maybe_cleanup_mesh(Mesh(split.x, split.y), is_cleanup);
}

Nested<TV2> offset(T a, Nested<TV2> c) {
  // auto merge = concat_polygons(c0, invert_polygon(c1));
  // return split_polygons(merge.x, merge.y, 0);
  return c;
}

Mesh union_add(Mesh mesh0, Mesh mesh1, bool is_cleanup) {
  auto merge = concat_meshes(mesh0, mesh1);
  // printf("--- START UNIONING ---\n");
  auto res   = split_mesh(merge, 0, is_cleanup);
  // printf("--- DONE  UNIONING ---\n");
  return res;
}

void pretty_print_num(T n) {
  printf("%g\n", n);
}
void pretty_print_v2d(TV2 pt) {
  printf("%g,%g\n", pt.x, pt.y);
}
void pretty_print_v3d(TV3 pt) {
  printf("%g,%g,%g\n", pt.x, pt.y, pt.z);
}
void pretty_print_v3i(IV3 pt) {
  printf("%d,%d,%d\n", pt.x, pt.y, pt.z);
}

void pretty_print_mesh(Mesh mesh) {
  int i = 0;
  printf("AREA %f VOLUME %f\n", mesh.soup->surface_area(RawArray<const TV3>(mesh.points)), mesh.soup->volume(RawArray<const TV3>(mesh.points)));
  for (auto pt : mesh.points) {
    printf("PT[%2d] %f,%f,%f\n", i, pt.x, pt.y, pt.z);
    i += 1;
  }
  i = 0;
  for (auto tri : mesh.soup->elements) {
    // T a = Triangle::area(mesh.points[tri.x], mesh.points[tri.y], mesh.points[tri.z]);
    auto p0 = mesh.points[tri.x];
    auto p1 = mesh.points[tri.y];
    auto p2 = mesh.points[tri.z];
    T a = 0.5 * cross(p1 - p0, p2 - p0).magnitude();
    printf("TRI[%d] %2d,%2d,%2d [%f,%f,%f] [%f,%f,%f] [%f,%f,%f] AREA %f\n",
           i, tri.x, tri.y, tri.z, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, a);
    i += 1;
  }
}

void pretty_print_matrix(Matrix<T,4> M) {
  printf("%g, %g, %g, %g\n", M.x[0][0],M.x[1][0],M.x[2][0],M.x[3][0]);
  printf("%g, %g, %g, %g\n", M.x[0][1],M.x[1][1],M.x[2][1],M.x[3][1]);
  printf("%g, %g, %g, %g\n", M.x[0][2],M.x[1][2],M.x[2][2],M.x[3][2]);
  printf("%g, %g, %g, %g\n", M.x[0][3],M.x[1][3],M.x[2][3],M.x[3][3]);
}

void pretty_print_array_v2d(Array<TV2> line) {
  int i = 0;
  for (auto pt : line) {
    printf("PT[%2d] %g,%g\n", i, pt.x, pt.y);
    i += 1;
  }
}

void pretty_print_array_v3d(Array<TV3> line) {
  int i = 0;
  for (auto pt : line) {
    printf("PT[%2d] %g,%g,%g\n", i, pt.x, pt.y, pt.z);
    i += 1;
  }
}

void pretty_print_array_v3i(Array<IV3> line) {
  int i = 0;
  for (auto pt : line) {
    printf("PT[%2d] %d,%d,%d\n", i, pt.x, pt.y, pt.z);
    i += 1;
  }
}

void pretty_print_poly(Nested<TV2> poly) {
  int i = 0;
  for (auto elt : poly) {
    Array<TV2> contour; for (auto e : elt) contour.append(e);
    printf("CONTOUR %d\n", i);
    pretty_print_array_v2d(contour);
    i += 1;
  }
}

void pretty_print_nested_v3d(Nested<TV3> polyline) {
  int i = 0;
  for (auto elt : polyline) {
    Array<TV3> line; for (auto e : elt) line.append(e);
    printf("LINE %d\n", i);
    pretty_print_array_v3d(line);
    i += 1;
  }
}

void pretty_print_nested_v2d(Nested<TV2> polyline) {
  int i = 0;
  for (auto elt : polyline) {
    Array<TV2> line; for (auto e : elt) line.append(e);
    printf("LINE %d\n", i);
    pretty_print_array_v2d(line);
    i += 1;
  }
}

std::string num_to_str (T num) {
  std::stringstream ss;
  ss << num;
  return ss.str();
}

std::string v2d_to_str (TV2 pt) {
  std::stringstream ss;
  ss << "[" << pt.x << "," << pt.y << "]";
  return ss.str();
}

std::string v3d_to_str (TV3 pt) {
  std::stringstream ss;
  ss << "[" << pt.x << "," << pt.y << "," << pt.z << "]";
  return ss.str();
}

std::string v3i_to_str (IV3 pt) {
  std::stringstream ss;
  ss << "[" << pt.x << "," << pt.y << "," << pt.z << "]";
  return ss.str();
}

std::string array_v2d_to_str(Array<TV2> line) {
  std::stringstream ss;
  int i = 0;
  ss << "[";
  for (auto pt : line) {
    if (i > 0) ss << ", ";
    ss << v2d_to_str(pt);
    i += 1;
  }
  ss << "]";
  return ss.str();
}

std::string array_v3d_to_str(Array<TV3> line) {
  std::stringstream ss;
  int i = 0;
  ss << "[";
  for (auto pt : line) {
    if (i > 0) ss << ", ";
    ss << v3d_to_str(pt);
    i += 1;
  }
  ss << "]";
  return ss.str();
}

std::string array_v3i_to_str(Array<const IV3> line) {
  std::stringstream ss;
  int i = 0;
  ss << "[";
  for (auto pt : line) {
    if (i > 0) ss << ", ";
    ss << v3i_to_str(pt);
    i += 1;
  }
  ss << "]";
  return ss.str();
}

std::string mesh_to_str(Mesh mesh) {
  std::stringstream ss;
  ss << "mesh(";
  ss << array_v3d_to_str(mesh.points);
  ss << ", ";
  ss << array_v3i_to_str(mesh.soup->elements);
  ss << ")";
  return ss.str();
}

std::string poly_to_str(Nested<TV2> poly) {
  std::stringstream ss;
  int i = 0;
  ss << "poly(";
  for (auto elt : poly) {
    Array<TV2> contour; for (auto e : elt) contour.append(e);
    Array<TV2> line; for (auto e : elt) line.append(e);
    if (i > 0) ss << ", ";
    ss << array_v2d_to_str(line);
    i += 1;
  }
  ss << ")";
  return ss.str();
}

std::string nested_v3d_to_str(Nested<TV3> polyline) {
  std::stringstream ss;
  int i = 0;
  ss << "[";
  for (auto elt : polyline) {
    Array<TV3> line; for (auto e : elt) line.append(e);
    if (i > 0) ss << ", ";
    ss << array_v3d_to_str(line);
    i += 1;
  }
  ss << "]";
  return ss.str();
}

std::string nested_v2d_to_str(Nested<TV2> polyline) {
  std::stringstream ss;
  int i = 0;
  ss << "[";
  for (auto elt : polyline) {
    Array<TV2> line; for (auto e : elt) line.append(e);
    if (i > 0) ss << ", ";
    ss << array_v2d_to_str(line);
    i += 1;
  }
  ss << "]";
  return ss.str();
}

std::string matrix_to_str(Matrix<T,4> M) {
  std::stringstream ss;
  ss << "mat(";
  ss << M.x[0][0] << "," << M.x[1][0] << "," << M.x[2][0] << "," << M.x[3][0] << ",";
  ss << M.x[0][1] << "," << M.x[1][1] << "," << M.x[2][1] << "," << M.x[3][1] << ",";
  ss << M.x[0][2] << "," << M.x[1][2] << "," << M.x[2][2] << "," << M.x[3][2] << ",";
  ss << M.x[0][3] << "," << M.x[1][3] << "," << M.x[2][3] << "," << M.x[3][3] << ")";
  ss << ")";
  return ss.str();
}

Mesh intersection(Mesh mesh0, Mesh mesh1, bool is_cleanup) {
  auto merge = concat_meshes(mesh0, mesh1);
  return split_mesh(merge, 1, is_cleanup);
}

Mesh mesh_from(int start, Mesh mesh) {
  Array<TV3> pts;
  for (int i = start; i < mesh.points.size(); i++) {
    pts.append(mesh.points[i]);
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
  Array<IV3> faces;
  for (auto tri : mesh.soup->elements) {
    bool is_all = tri.x >= start && tri.y >= start && tri.z >= start;
    if (is_all)
      faces.append(vec(tri.x - start, tri.y - start, tri.z - start));
  }
  return fab_mesh(faces, pts);
}

Nested<TV2> slice(T z, Mesh mesh) {
  T t = 1e8; // TODO: USE INF
  auto smesh = const_mesh(cube_mesh(vec(-t, -t, -t), vec( t,  t,  z)));
  // printf("-- SLICE %f MESH --\n", z);
  // pretty_print_mesh(smesh);
  // printf("-- ABOUT TO SLICE %f --\n", z);
  auto res = intersection(smesh, mesh, false);
  // printf("-- SLICE %f INTERSECTION MESH --\n", z);
  // pretty_print_mesh(res);
  int start = smesh.points.size() + mesh.points.size();

  // printf("---->>>>\n");
  // print_mesh(res);
  // printf("========\n");
  // printf("%d OLD %d NEW VERTICES\n", start, res.points.size() - start);
  // for (auto tri : res.soup->elements) {
  //   bool is_all = tri.x >= start && tri.y >= start && tri.z >= start;
  //   if (is_all)
  //     printf("SLICE %2d,%2d,%2d\n", tri.x, tri.y, tri.z);
  // }

  auto new_mesh = mesh_from(start, res);
  auto boundary = new_mesh.soup->boundary_mesh();
  auto polys = boundary->polygons();

  // int i = 0;
  // for (auto p : new_mesh.points) {
  //   printf("BP[%2d] %f,%f,%f\n", i, p.x, p.y, p.z);
  //   i += 1;
  // }
  // for (auto elt : boundary->elements) {
  //   printf("ELT %d,%d\n", elt.x, elt.y);
  // }
  // printf("POLYGONS\n");
  // for (auto poly : polys.x) {
  //   for (auto p : poly) 
  //     printf("%d ", p);
  //   printf("\n");
  // }

  Nested<TV2,false> pres;
  for (auto poly : polys.x) {
    Array<TV2> contour;
    for (auto p : poly) {
      auto pt = new_mesh.points[p];
      contour.append(vec(pt.x, pt.y));
    }
    pres.append(contour);
  }
  pres.freeze();
  return cleanup_poly(pres);
}

Mesh difference(Mesh mesh0, Mesh mesh1, bool is_cleanup) {
  auto merge = concat_meshes(mesh0, invert_mesh(mesh1));
  return split_mesh(merge, 0, is_cleanup);
}

template<class ET>
Nested<ET> array_to_nested(Array<ET> contour) {
  Nested<ET,false> poly;
  poly.append(contour);
  poly.freeze();
  return poly;
}

template<class ET>
Array<ET> nested_elt(Nested<ET> poly, int i) {
  Array<ET> contour;
  // printf("POLY TO CONTOUR[%d] %d\n", i, poly[i].size());
  for (auto e : poly[i])
    contour.append(e);
  return contour;
}

Mesh all_mesh(void) {
  T x = 1e8;
  return const_mesh(cube_mesh(vec(-x, -x, -x), vec(x, x, x)));
}

Mesh none_mesh(void) {
  Array<IV3> faces;
  Array<TV3> points;
  return fab_mesh(faces, points);
}

Nested<TV2> square_poly(TV2 min, TV2 max) {
  Array<Vector<real, 2>> pts;
  pts.append(vec(max.x,min.y));
  pts.append(vec(max.x,max.y));
  pts.append(vec(min.x,max.y));
  pts.append(vec(min.x,min.y));
  return array_to_nested(pts);
}

Nested<TV2> square_poly(T diam) {
  auto rad = 0.5 * diam;
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

Nested<TV2> circle_poly(T diam, int n) {
  auto rad = 0.5 * diam;
  Array<Vector<real, 2>> pts;
  for (int i = n-1; i >= 0; i--) {
    T a = (2 * M_PI * i) / n;
    // printf("A = %f\n", a);
    pts.append(vec(rad * sin(a),rad * cos(a)));
  }
  return array_to_nested(pts);
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
  return array_to_nested(pts);
}

struct Meshy {
public:
  Array<TV3> points;
  Array<int> indices;
};

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

Mesh triangulate (Nested<TV3> poly) { 
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
  Array<IV3> faces;
  for (int i = 0; i < mesh.indices.size()/3; i++) {
    faces.append(vec(mesh.indices[i*3], mesh.indices[i*3 + 2], mesh.indices[i*3 + 1]));
  }
  return fab_mesh(faces, mesh.points);
}

Mesh triangulate (Nested<TV2> poly2) { 
  return triangulate(nested2_to_nested3(poly2));
}

/*
Tuple<Ref<TriangleSoup>,Array<TV3>> triangulate(Array<TV2> poly) {
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

Mesh mul(Matrix<T,4> m, Mesh mesh, bool is_invert) {
  auto res = Mesh(mesh.soup, mul(m, mesh.points));
  return is_invert ? invert_mesh(res) : res;
}

Mesh cone_mesh(T len, Array<TV2> poly) {
  Array<TV3> bot;
  Array<TV3> vh;

  T zmin = - len / 2.0;
  T zmax =   len / 2.0;
  auto c = vec(0.0, 0.0);
  int n_boundary = poly.size();
  for (auto elt : poly) {
    bot.append(vec(elt.x, elt.y, zmin));
    c = c + elt;
  }
  c = c / (double)n_boundary;
  auto lid = triangulate(array_to_nested(bot));

  vh = lid.points;
  vh.append(vec(c.x, c.y, zmax));
    
  Array<IV3> faces;
  for (auto face : lid.soup->elements)
    faces.append(face);
                 
  for (int i = 0; i < n_boundary; i++) {
    int ni = (i + 1)%n_boundary;
    faces.append(vec(i, ni, n_boundary));
  }
  return fab_mesh(faces, vh);
}

Mesh cone_mesh(T len, Nested<TV2> contours) {
  ensure(contours.size() == 1, "CONE ONLY WORKS ON SINGLE CONTOUR POLYGONS");
  return cone_mesh(len, nested_elt(contours, 0));
}

Mesh mesh_between_poly(Array<TV3> bot, Array<TV3> top) {
  Array<TV3> vh;

  // printf("BOT\n");
  // for (auto e : bot) 
  //   printf("  %f,%f,%f\n", e.x, e.y, e.z);
  // printf("TOP\n");
  // for (auto e : top) 
  //   printf("  %f,%f,%f\n", e.x, e.y, e.z);
  int n_boundary = bot.size();
  // printf("N_BOUNDARY = %d\n", n_boundary);
  auto top_mesh = triangulate(array_to_nested(top));
  int n_mesh = top_mesh.points.size();
  // printf("N_MESH = %d\n", n_mesh);

  auto bot_mesh = triangulate(array_to_nested(bot));

  auto lids = concat_meshes(top_mesh, invert_mesh(bot_mesh));

  // for (auto pt : vh)
  //   printf("PT [%f,%f,%f]\n", pt.x, pt.y, pt.z);

  Array<IV3> faces;
  for (auto face : lids.soup->elements)
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
  return fab_mesh(faces, lids.points);
}

Mesh taper_mesh(T len, T r0, T r1, Array<TV2> poly) {
  Array<TV3> top, bot;
  T zmin = - len / 2.0;
  T zmax =   len / 2.0;
  for (auto elt : poly)
    top.append(vec(r0 * elt.x, r0 * elt.y, zmin));
  for (auto elt : poly) 
    bot.append(vec(r1 * elt.x, r1 * elt.y, zmax));
  return mesh_between_poly(bot, top);
}

Mesh taper_mesh(T len, T r0, T r1, Nested<TV2> contours) {
  if (contours.size() == 1) {
    if (true || contours.size() == 1)
      return taper_mesh(len, r0, r1, nested_elt(contours, 0));
    else
      return taper_mesh(len, r0, r1, nested_elt(contours, 1));
  } else {
    auto res = none_mesh();
    for (int i = 0; i < contours.size(); i++) {
      auto contour = nested_elt(contours, i);
      if (!is_clockwise(contour)) {
        // printf("OUTER %d\n", i);
        auto mesh = taper_mesh(len, r0, r1, contour);
        // pretty_mesh(mesh);
        res = union_add(res, mesh);
      }
    }
    for (int i = 0; i < contours.size(); i++) {
      auto contour = nested_elt(contours, i);
      if (is_clockwise(contour)) {
        // printf("INNER %d\n", i);  
        auto mesh = mul(scale_matrix(vec(1.0,1.0,2.0)), taper_mesh(len, r0, r1, contour));
        // pretty_print_mesh(mesh);
        res = union_add(res, mesh);
      }
    }
    return res;
  }
}

Mesh extrude(T len, Nested<TV2> poly) {
  return taper_mesh(len, 1.0, 1.0, poly);
}

Mesh revolve(int n, Array<TV2> poly) {
  Array<TV3> pts;

  for (int i = 0; i < n; i++) {
    T a = (2 * M_PI * i) / n;
    for (auto elt : poly) {
      auto rad = elt.x;
      pts.append(vec(rad * sin(a), elt.y, rad * cos(a)));
    }
  }

  int m = poly.size();

  Array<IV3> faces;

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
    return revolve(n, nested_elt(contours, 0));
  else {
    auto res = none_mesh();
    for (int i = 0; i < contours.size(); i++) {
      auto contour = nested_elt(contours, i);
      if (!is_clockwise(contour))
        res = union_add(res, revolve(n, contour));
    }
    for (int i = 0; i < contours.size(); i++) {
      auto contour = nested_elt(contours, i);
      if (is_clockwise(contour)) 
        res = union_add(res, revolve(n, contour));
    }
    return res;
    // return revolve(n, nested_elt(contours));
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
    return array_to_nested(points);
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
      res = union_add(res, edge);
    }
    j += 1;
  }
  // printf("THICKEN DONE\n");
  return res;
}

Mesh fat_triangle(T rad, TV p0, TV p1, TV p2) {
  auto n = rad * cross(p1 - p0, p2 - p0).normalized();
  Array<TV3> bot, top;
  // printf("NORMAL = %f,%f,%f\n", n.x, n.y, n.z);
  bot.append(p0 - n); bot.append(p1 - n); bot.append(p2 - n);
  top.append(p0 + n); top.append(p1 + n); top.append(p2 + n);
  return mesh_between_poly(top, bot);
}

Mesh fat_edge(int n, T rad, TV from, TV to) {
  auto v = to - from;
  auto res = extrude(magnitude(v), circle_poly(2 * rad, n));
  auto c = (to + from) * 0.5;
  return mul(translation_matrix(c), mul(rotation_matrix(vec(0.0, 0.0, 2*rad), v), res));
}

Nested<TV2> fat_edge(int n, T rad, TV2 from, TV2 to) {
  auto v = to - from;
  auto res = square_poly(vec(-rad, -magnitude(v)/2), vec(rad, magnitude(v)/2));
  auto c = (to + from) * 0.5;
  return mul(translation_matrix(c), mul(rotation_matrix(vec(0.0, magnitude(v)), v), res));
}

Nested<TV2> fat_dot(int n, T rad, TV2 pt) {
  auto circ = circle_poly(2 * rad, n);
  auto res = mul(translation_matrix(pt), circ);
  return res;
}

Nested<TV2> offset_poly(int n, T rad, Nested<TV2> poly) {
  Nested<TV2> res = poly;
  for (auto contour : poly) {
    for (int i = 0; i < contour.size(); i++) {
      res = union_add(res, union_add(fat_dot(n, rad, contour[i]), fat_edge(n, rad, contour[i], contour[(i+1)%contour.size()])));
    }
  }
  return res;
}

Nested<TV2> offset_polyline(int n, T rad, Nested<TV2> polyline) {
  Nested<TV2> res = none_poly();
  for (auto line : polyline) {
    for (auto pt : line) {
      res = union_add(res, fat_dot(n, rad, pt));
    }
    for (int i = 0; i < (line.size()-1); i++) {
      res = union_add(res, fat_edge(n, rad, line[i], line[(i+1)%line.size()]));
    }
  }
  return res;
}

Mesh thicken(int n, T rad, Mesh mesh) {
  Mesh res = const_mesh(sphere_mesh(n, mesh.points[0], rad));
  // printf("TH %f,%f,%f\n", line[0].x, line[0].y, line[0].z);

  for (int i = 1; i < mesh.points.size(); i++) {
    // printf("TH %f,%f,%f\n", line[i].x, line[i].y, line[i].z);
    res = union_add(res, const_mesh(sphere_mesh(n, mesh.points[i], rad)));
  }
  auto edges = mesh.soup->triangle_edges();
  for (auto edge : edges)
    printf("EDGE %d to %d\n", edge.x, edge.y);
  for (auto edge : edges) {
    res = union_add(res, fat_edge(16, rad, mesh.points[edge.x], mesh.points[edge.y]));
  }
  
  return res;
}

Mesh thicken(int n, T rad, Nested<TV3> polyline) {
  Mesh res = none_mesh();
  for (auto line : polyline) {
    // printf("TH %f,%f,%f\n", line[0].x, line[0].y, line[0].z);
    for (int i = 0; i < line.size(); i++) {
      // printf("TH %f,%f,%f\n", line[i].x, line[i].y, line[i].z);
      res = union_add(res, const_mesh(sphere_mesh(n, line[i], rad)));
    }
    for (int i = 1; i < line.size(); i++) {
      res = union_add(res, fat_edge(12, rad, line[i-1], line[i]));
    }
  }
  
  return res;
}

Mesh check_mesh(Mesh mesh) {
  auto topo = new_<TriangleTopology>(mesh.soup->elements);
  printf("VERTICES %d EDGES %d BOUNDARD-EDGES %d FACES %d IS-MANIFOLD %d HAS-BOUNDARY %d\n",
         topo->n_vertices(), topo->n_edges(), topo->n_boundary_edges(), topo->n_faces(),
         topo->is_manifold(), topo->has_boundary());
  topo->dump_internals();
  return mesh;
}

/*
Mesh offset_mesh(int n, T rad, Mesh mesh) {
  Array<TV3> pos(mesh.y);
  auto topo = new_<TriangleTopology>(mesh.x->elements);
  auto new_mesh = rough_offset_mesh(topo, RawField<const TV3,VertexId>(pos), rad);
  auto new_soup = new_mesh.x->face_soup().x;
  Array<TV3> new_points;
  for (auto point : new_mesh.y.flat)
    new_points.append(point);
  return tuple(const_mesh(new_soup), new_points);
  return mesh;
}

Mesh offset_mesh(int n, T rad, Mesh mesh) {
  // Mesh res = none_mesh();
  Mesh res = mesh;
  auto edges = mesh.soup->segment_soup();

  pretty_print_mesh(mesh);
  int k = 0;
  for (auto face : mesh.soup->elements) {
    printf("FACE %d/%d\n", k, mesh.soup->elements.size());
    // res = union_add(res, dither_mesh(const_mesh(sphere_mesh(n, pt, rad)), 1e-3));
    auto fatty = fat_triangle(rad, mesh.points[face.x], mesh.points[face.y], mesh.points[face.z]);
    // return fatty;
    res = union_add(res, dither_mesh(fatty, 1e-3));
    k += 1;
  }
  int j = 0;
  for (auto pt : mesh.points) {
    printf("POINT %d/%d\n", j, mesh.points.size());
    // res = union_add(res, dither_mesh(const_mesh(sphere_mesh(n, pt, rad)), 1e-3));
    res = union_add(res, const_mesh(sphere_mesh(n, pt, rad)));
    j += 1;
  }
  int i = 0;
  for (auto edge : edges->elements) {
    printf("EDGE %d/%d %d to %d\n", i, edges->elements.size(), edge.x, edge.y);
    // res = union_add(res, dither_mesh(fat_edge(8, rad, mesh.points[edge.x], mesh.points[edge.y]), 1e-3));
    res = union_add(res, fat_edge(8, rad, mesh.points[edge.x], mesh.points[edge.y]));
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
    res.extend(mul(tmat * smat, letter_outlines[(int)txt[i]]));
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

