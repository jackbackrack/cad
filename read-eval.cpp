#include "cad.h"
#include "geom.h"
#include "read-eval.h"

static bool is_init_tokenizer = false;

static int identChars[256];

void init_tokenizer( void ) {
  if (!is_init_tokenizer) {
    is_init_tokenizer = true;
    std::string chars("_");
    for (size_t i = 0; i < 256; i++) 
      identChars[i] = isalnum(i);
    for (size_t i = 0; i < chars.size(); i++)
      identChars[(size_t)chars[i]] = true;
  }
}

Token Tokenizer::do_get() {
  init_tokenizer();
  std::string str;
  bool isHex = false;
  int      c = ' ';
  while (isspace(c)) { // skip whitespace
    c = s->get();
  }
  if (c == '"') { // str?
    while (s->peek() > 0 && s->peek() != '"') {
      str += s->get();
    }
    s->get();
    return Token(str, true);
  }
  if (!(isdigit(c) || c == '.' || c == '-') && identChars[c]) { // sym?
    str += c;
    while (s->peek() > 0 && identChars[s->peek()]) {
      str += s->get();
    }
    return Token(str, false);
  }
  if (isdigit(c) || c == '.' || c == '-') {   // num?
    str += c;
    // TODO: REAL LEXER
    while (isdigit(s->peek()) || s->peek() == '.' || isalpha(s->peek())) { 
      char c = s->get();
      isHex = isHex || c == 'x';
      str += c;
    }
    double temp;
    if (isHex) {
      uint32_t hex = 0;
      for (uint32_t i = 2; i < str.size(); i++) {
        char c  = str[i];
        int dig = isdigit(c) ? (c - '0') : (c - 'a' + 10);
        hex = hex * 16 + dig;
      }
      temp = hex;
    } else {
      if (str == "-")
        return Token((TokenType)c);
      std::stringstream ssout(str);
      ssout >> temp;
    }
    return Token(temp);
  }
  if (c == ';') { // comment?
    do c = s->get();
    while (!s->eof() && c != '\n' && c != '\r');
    if (!s->eof())
      return get();
  }
  if (s->eof()) { // eof?
    return Token(tok_eof);
  }
  return Token((TokenType)c); // otherwise just return char as token
}

Token Tokenizer::get() {
  init_tokenizer();
  Token res;
  if (token.type != tok_null) {
    Token tok = token;
    token.clear();
    res = tok;
  } else {
    res = do_get();
  }
  // printf("GET = %s\n", res.to_str().c_str());
  return res;
}

std::string Tokenizer::get_sym() {
  auto t = get();
  if (t.type == tok_sym) {
    return t.sym;
  } else {
    error("GET-SYM: error");
    return "";
  }
}

double Tokenizer::get_double() {
  auto t = get();
  if (t.type == tok_num) {
    return t.num;
  } else {
    error("GET-DOUBLE: error");
    return -1;
  }
}

int Tokenizer::get_int() {
  return (int)get_double();
}

char Tokenizer::get_char() {
  auto tok = get();
  if (tok.type != tok_sym)
    error("EXPECT: FAILURE NOT CHAR");
  return tok.sym[0];
}

void Tokenizer::expect(char kind) {
  auto tok = get();
  if (tok.type != kind) {
    std::string kstr; kstr += kind;
    error("EXPECT: FAILURE NOT KIND ", kstr);
  }
}

void Tokenizer::expect(std::string name) {
  auto tok = get();
  if (tok.type != tok_sym || tok.sym != name)
    error("EXPECT: FAILURE NOT STRING ", name);
}

Token Tokenizer::peek() {
  init_tokenizer();
  if (token.type != tok_null) {
    return token;
  } else {
    return token = do_get();
  }
}


extern Geom* parse_expression(Tokenizer& s);

std::vector< Geom* > parse_args(Tokenizer& s, char end_char) {
  std::vector< Geom* > res;
  bool is_first = true;
  for (;;) {
    auto tok = s.peek();
    if (tok.type == end_char) {
      s.expect(end_char);
      return res;
    } else {
      if (!is_first)
        s.expect(',');
      res.push_back(parse_expression(s));
      is_first = false;
    }
  }
}

Geom* fold_bin_op (binary_geom_op_t op, uint32_t i, std::vector< Geom* > &args) {
  if (i == (args.size()-1))
    return args[i];
  else 
    return op(args[i], fold_bin_op(op, i+1 ,args));
}

Geom* parse_factor(Tokenizer& s) {
  auto tok1 = s.peek();
  // printf("PARSING FACTOR %s\n", tok1.to_str().c_str());
  if (tok1.type == '!') {
    s.get();
    return g_not(parse_factor(s));
    // } else if (tok1.type == '-') {
    // s.get();
    // return neg(parse_factor(s));
  } else if (tok1.type == '(') {
    s.expect('(');
    auto g = parse_expression(s);
    s.expect(')');
    return g;
  } else if (tok1.type == '[') {
    s.expect('[');
    auto args = parse_args(s, ']');
    if (args.size() == 0) {
      error("Illegal zero vec args\n"); return NULL;
    } else {
      if (all_args_kind(args, num_kind)) {
        if (args.size() == 3) {
          return new V3dGeom(vec(g_num_val(args[0]), g_num_val(args[1]), g_num_val(args[2])));
        } else if (args.size() == 2) {
          return new V2dGeom(vec(g_num_val(args[0]), g_num_val(args[1])));
        } else
          error("Wrong number of vec args\n");
      } else if (all_args_kind(args, v2d_kind)) {
        return g_array_v2d(args);
      } else if (all_args_kind(args, v3d_kind)) {
        return g_array_v3d(args);
      } else if (all_args_kind(args, v3i_kind)) {
        return g_array_v3i(args);
      } else if (all_args_kind(args, array_v2d_kind)) {
        return g_nested_v2d(args);
      } else if (all_args_kind(args, array_v3d_kind)) {
        return g_nested_v3d(args);
      } else {
        error("Illegal vec args\n"); return NULL;
      }
    }
  } else if (tok1.type == tok_str) {
    s.get();
    return new StringGeom(tok1.sym);
  } else if (tok1.type == tok_sym) {
    s.get();
    if (tok1.sym == "pi")
      return g_pi();
    else if (tok1.sym == "all") 
      return g_all();
    else if (tok1.sym == "none") 
      return g_none();
    else if (tok1.sym == "all2") 
      return g_all2();
    else if (tok1.sym == "none2") 
      return g_none2();
    else {
      auto tok2 = s.peek();
      Geom* g = NULL;
      if (tok2.type == '(') {
        s.get();
        {
          auto args = parse_args(s, ')');
          if (tok1.sym == "xrot") {
            g = g_xrot(args[0], args[1]);
            // g = g_reflect_yz(args[0]);
          } else if (tok1.sym == "mat") {
            if (args.size() == 9) {
              return new MatGeom(
                Matrix<T,4>(g_num_val(args[ 0]), g_num_val(args[ 1]), g_num_val(args[ 2]), 0,
                            g_num_val(args[ 3]), g_num_val(args[ 4]), g_num_val(args[ 5]), 0,
                            g_num_val(args[ 6]), g_num_val(args[ 7]), g_num_val(args[ 8]), 0,
                            g_num_val(args[ 9]), g_num_val(args[10]), g_num_val(args[11]), 1));
            } else if (args.size() == 16) {
              return new MatGeom(
                Matrix<T,4> (
                  g_num_val(args[ 0]), g_num_val(args[ 1]), g_num_val(args[ 2]), g_num_val(args[ 3]),
                  g_num_val(args[ 4]), g_num_val(args[ 5]), g_num_val(args[ 6]), g_num_val(args[ 7]),
                  g_num_val(args[ 8]), g_num_val(args[ 9]), g_num_val(args[10]), g_num_val(args[11]),
                  g_num_val(args[12]), g_num_val(args[13]), g_num_val(args[14]), g_num_val(args[15])));
            } else
              error("Wrong number of mat args\n");
          } else if (tok1.sym == "array_v3d") {
            return g_array_v3d(args);
          } else if (tok1.sym == "array_v3i") {
            return g_array_v3i(args);
          } else if (tok1.sym == "array_v2d") {
            return g_array_v2d(args);
          } else if (tok1.sym == "nested_v3d") {
            return g_nested_v3d(args);
          } else if (tok1.sym == "nested_v2d") {
            return g_nested_v2d(args);
          } else if (tok1.sym == "poly" || tok1.sym == "polygon") {
            return g_poly(args);
          } else if (tok1.sym == "mesh") {
            auto mesh = fab_mesh(g_array_v3i_val(args[1]), g_array_v3d_val(args[0]));
            return new MeshGeom(mesh);
          } else if (tok1.sym == "elt") {
            return g_elt(args[0], args[1]);
          } else if (tok1.sym == "letter") {
            return g_letter(args[0]);
          } else if (tok1.sym == "text") {
            return g_text(args[0]);
          } else if (tok1.sym == "thicken") {
            return g_thicken(args[0], args[1]);
          } else if (tok1.sym == "offset") {
            return g_offset(args[0], args[1]);
          } else if (tok1.sym == "hollow") {
            g = g_hollow(args[0], args[1]);
          } else if (tok1.sym == "rot") {
            if (args.size() == 2)
              g = g_rot(args[0], args[1]);
            else if (args.size() == 3)
              g = g_rot_from_to(args[0], args[1], args[2]);
            else {
              error("WRONG NUM OF ARGS FOR ROT");
            }
          } else if (tok1.sym == "yrot") {
            g = g_yrot(args[0], args[1]);
          } else if (tok1.sym == "zrot") {
            g = g_zrot(args[0], args[1]);
          } else if (tok1.sym == "mov") {
            g = g_mov(args[0], args[1]);
          } else if (tok1.sym == "xmov") {
            g = g_xmov(args[0], args[1]);
          } else if (tok1.sym == "ymov") {
            g = g_ymov(args[0], args[1]);
          } else if (tok1.sym == "zmov") {
            g = g_zmov(args[0], args[1]);
          } else if (tok1.sym == "mag") {
            g = g_mag(args[0], args[3]);
          } else if (tok1.sym == "mag1") {
            g = g_mag1(args[0], args[1]);
          } else if (tok1.sym == "xmag") {
            g = g_xmag(args[0], args[1]);
          } else if (tok1.sym == "ymag") {
            g = g_ymag(args[0], args[1]);
          } else if (tok1.sym == "zmag") {
            g = g_zmag(args[0], args[1]);
          } else if (tok1.sym == "slice") {
            g = g_slice(args[0], args[1]);
          } else if (tok1.sym == "bbox") {
            g = g_bbox(args[0]);
          } else if (tok1.sym == "dims") {
            g = g_dims(args[0]);
          } else if (tok1.sym == "center") {
            g = g_center(args[0]);
          } else if (tok1.sym == "centering") {
            g = g_centering(args[0]);
          } else if (tok1.sym == "simplify") {
            g = g_simplify(args[0]);
          } else if (tok1.sym == "cleanup") {
            g = g_cleanup(args[0]);
          } else if (tok1.sym == "extrude") {
            g = g_extrude(args[0], args[1]);
          } else if (tok1.sym == "cone") {
            g = g_cone(args[0], args[1]);
          } else if (tok1.sym == "sphere") {
            g = g_sphere(args[0]);
          } else if (tok1.sym == "circle") {
            g = g_circle(args[0]);
          } else if (tok1.sym == "print") {
            g = g_print(args[0]);
          } else if (tok1.sym == "check") {
            g = g_check(args[0]);
          } else if (tok1.sym == "pprint") {
            g = g_pretty_print(args[0]);
          } else if (tok1.sym == "square") {
            if (args.size() == 2)
              g = g_square_lo_hi(args[0], args[1]);
            else
              g = g_square(args[0]);
          } else if (tok1.sym == "cube") {
            if (args.size() == 2)
              g = g_cube_lo_hi(args[0], args[1]);
            else
              g = g_cube(args[0]);
          } else if (tok1.sym == "reflect_x") {
            g = g_reflect_x(args[0]);
          } else if (tok1.sym == "reflect_y") {
            g = g_reflect_y(args[0]);
          } else if (tok1.sym == "reflect_z") {
            g = g_reflect_z(args[0]);
          } else if (tok1.sym == "reflect_xy") {
            g = g_reflect_xy(args[0]);
          } else if (tok1.sym == "reflect_xz") {
            g = g_reflect_xz(args[0]);
          } else if (tok1.sym == "reflect_yz") {
            g = g_reflect_yz(args[0]);
          } else if (tok1.sym == "revolve") {
            g = g_revolve(args[0]);
          } else if (tok1.sym == "hull") {
            g = g_hull(args[0]);
          } else if (tok1.sym == "taper") {
            g = g_taper(args[0], args[1], args[2], args[3]);
          } else if (tok1.sym == "load") {
            g = g_load(args[0]);
          } else if (tok1.sym == "save") {
            g = g_save(args[0], args[1]);
            /*
          } else if (tok1.sym == "rect") {
            g = g_rect(args[0], args[1], args[2], args[3]);
          } else if (tok1.sym == "shear") {
            g = g_shear_x_z(args[0], args[1], args[2], args[3], args[4]);
          } else if (tok1.sym == "transform") {
            g = g_xform(args[0], args[1], args[2], args[3]);
          } else if (tok1.sym == "rect3") {
            g = g_rect3(args[0], args[1], args[2], args[3], args[4], args[5]);
          } else if (tok1.sym == "space") {
            g = g_space(args[0]);
          } else if (tok1.sym == "sin") {
            g = g_sin(args[0]);
          } else if (tok1.sym == "cos") {
            g = g_cos(args[0]);
          } else if (tok1.sym == "tan") {
            g = g_tan(args[0]);
          } else if (tok1.sym == "abs") {
            g = g_abs(args[0]);
          } else if (tok1.sym == "half") {
            g = g_half(args[0], args[1], args[2], args[3]);

          } else if (tok1.sym == "xbox") {
            g = fold_bin_op(g_xbox, 0, args);
          } else if (tok1.sym == "align_xmin") {
            g = fold_bin_op(g_align_xmin, 0, args);
          } else if (tok1.sym == "align_xmax") {
            g = fold_bin_op(g_align_xmax, 0, args);
          } else if (tok1.sym == "align_ymin") {
            g = fold_bin_op(g_align_ymin, 0, args);
          } else if (tok1.sym == "align_ymax") {
            g = fold_bin_op(g_align_ymax, 0, args);
          } else if (tok1.sym == "align_zmin") {
            g = fold_bin_op(g_align_zmin, 0, args);
          } else if (tok1.sym == "align_zmax") {
            g = fold_bin_op(g_align_zmax, 0, args);
          } else if (tok1.sym == "ybox") {
            g = fold_bin_op(g_ybox, 0, args);
          } else if (tok1.sym == "zbox") {
            g = fold_bin_op(g_zbox, 0, args);
            */
          } else
            error("UNKNOWN FUNCTION", tok1.sym);
        }
        // s.expect(')');
        return g;
      }
    }
  } else if (tok1.type == tok_num) {
    s.get();
    return g_num(tok1.num);
  } else if (tok1.type == '(') {
    s.expect('(');
    auto g = parse_expression(s);
    s.expect(')');
    return g;
  } else {
    error("factor: syntax error");
  }
  return NULL;
}

bool is_term_sym(Token tok) {
  return tok.type == '*' || tok.type == '/';
}

Geom* parse_term(Tokenizer& s) {
  Geom* g = parse_factor(s);
  auto tok = s.peek();
  // printf("PARSING TERM %s\n", tok.to_str().c_str());
  while (is_term_sym(tok)) {
    s.get();
    if (tok.type == '*')
      g = g_mul(g, parse_factor(s));
    else if (tok.type == '/')
      g = g_div(g, parse_factor(s));
    tok = s.peek();
  }
  return g;
}

bool is_expression_sym(Token tok) {
  return tok.type == '+' || tok.type == '-'  || tok.type == '\\' || tok.type == '&' || tok.type == '|';
}

Geom* parse_geom(Tokenizer& s) {
  return parse_expression(s);
}

Geom* parse_geom(std::istringstream &ss) {
  Tokenizer s(&ss, 0);
  return parse_geom(s);
}

Geom* parse_geom(std::string s) {
  std::istringstream ss(s);
  auto g  = parse_geom(ss);
  return g;
}

Geom* parse_expression(Tokenizer& s) {
  auto g = parse_term(s);
  auto tok = s.peek();
  // printf("PARSING EXPRESSION %s\n", tok.to_str().c_str());
  while (is_expression_sym(tok)) {
    s.get();
    if (tok.type == '+')
      g = g_add(g, parse_expression(s));
    else if (tok.type == '-') 
      g = g_sub(g, parse_expression(s));
    else if (tok.type == '\\')
      g = g_difference(g, parse_expression(s));
    else if (tok.type == '&')
      g = g_intersection(g, parse_expression(s));
    else if (tok.type == '|')
      g = g_union(g, parse_expression(s));
    else
      error("PARSE EXPR: error");
    tok = s.peek();
  }
  return g;
}

inline TV3 get_color(TV3 n) {
  return vec((n.x > 0.0 ? n.x : 0.0) + (n.y < 0.0 ? -0.5*n.y : 0.0) + (n.z < 0.0 ? -0.5*n.z : 0.0),
             (n.y > 0.0 ? n.y : 0.0) + (n.z < 0.0 ? -0.5*n.z : 0.0) + (n.x < 0.0 ? -0.5*n.x : 0.0),
             (n.z > 0.0 ? n.z : 0.0) + (n.x < 0.0 ? -0.5*n.x : 0.0) + (n.y < 0.0 ? -0.5*n.y : 0.0));
}

int display_mesh (Mesh mesh, bool is_show_lines, bool is_show_normals) {
  int dl = glGenLists(1);
  glNewList(dl, GL_COMPILE);
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glBegin(GL_TRIANGLES);
  for (auto t : mesh.soup->elements) {
    auto p0 = mesh.points(t[0]);
    auto p1 = mesh.points(t[1]);
    auto p2 = mesh.points(t[2]);
    auto n = cross(p1 - p0, p2 - p0);
    normalize(n);
    glNormal3d(n.x, n.y, n.z);
    auto d = get_color(n);
    glColor4f(d.x, d.y, d.z, 1.0);
    glVertex3d(p0.x, p0.y, p0.z);
    glVertex3d(p1.x, p1.y, p1.z);
    glVertex3d(p2.x, p2.y, p2.z);
  }
  glEnd();
  for (auto t : mesh.soup->elements) {
    auto p0 = mesh.points(t[0]);
    auto p1 = mesh.points(t[1]);
    auto p2 = mesh.points(t[2]);
    if (is_show_lines) {
      glColor4f(1.0, 0.0, 0.0, 1.0);
      glBegin(GL_LINE_LOOP);
      glVertex3d(p0.x, p0.y, p0.z);
      glVertex3d(p1.x, p1.y, p1.z);
      glVertex3d(p2.x, p2.y, p2.z);
      glEnd();
    }
    if (is_show_normals) {
      auto n = cross(p1 - p0, p2 - p0);
      normalize(n);
      auto c = (p0 + p1 + p2) * 0.3333;
      auto d = c + n;
      auto dc = get_color(n);
      // glColor4f(0.0, 1.0, 0.0, 1.0);
      glColor4f(dc.x, dc.y, dc.z, 1.0);
      glBegin(GL_LINES);
      glVertex3d(c.x, c.y, c.z);
      glVertex3d(d.x, d.y, d.z);
      glEnd();
    }
  }
  glPointSize(4);
  glColor3f(1, 1, 1);
  glBegin(GL_POINTS);
  for (auto p : mesh.points) {
    glVertex3d(p.x, p.y, p.z);
  }
  glEnd();
  glEndList();
  return dl;
}

// CLIP forces the number x into the range [min,max]
inline double CLIP(double x, double mn, double mx) { return max(mn, min(mx, x)); }

TV3 hsv_to_rgb (double h, double s, double v) {
  double rt, gt, bt;
  s = CLIP(s, (double)0, (double)1);
  if (s == 0.0) {
    rt = gt = bt = v;
  } else {
    double h_temp = (h == 360.0) ? 0.0 : h;
    double f, p, q, t; 
    h_temp /= 60.0;
    int i = (int)h_temp;
    f = h_temp - i;
    p = v*(1-s);
    q = v*(1-(s*f));
    t = v*(1-(s*(1-f)));
    switch (i) {
    case 0: rt = v; gt = t; bt = p; break;
    case 1: rt = q; gt = v; bt = p; break;
    case 2: rt = p; gt = v; bt = t; break;
    case 3: rt = p; gt = q; bt = v; break;
    case 4: rt = t; gt = p; bt = v; break;
    case 5: rt = v; gt = p; bt = q; break;
    }
  }
  return vec(rt, gt, bt);
}

int display_poly (Nested<TV2> poly) {
  int dl = glGenLists(1);
  glNewList(dl, GL_COMPILE);
  int tot = 0;
  for (auto contour : poly)
    tot += contour.size();
  int i = 0;
  for (auto contour : poly) {
    glBegin(GL_LINE_LOOP);
    for (auto p : contour) {
      auto c = hsv_to_rgb(((double)i / tot) * 360.0, 1.0, 1.0);
      glColor4f(c.x, c.y, c.z, 1.0);
      glVertex2d(p.x, p.y);
      i += 1;
    }
    glEnd();
  }
  glEndList();
  return dl;
}

int display_polyline2 (Nested<TV2> polyline) {
  int dl = glGenLists(1);
  glNewList(dl, GL_COMPILE);
  glColor4f(1.0, 1.0, 1.0, 1.0);
  for (auto line : polyline) {
    glBegin(GL_LINE_STRIP);
    for (auto p : line) {
      glVertex2d(p.x, p.y);
    }
    glEnd();
  }
  glEndList();
  return dl;
}

int display_polyline3 (Nested<TV3> polyline) {
  int dl = glGenLists(1);
  glNewList(dl, GL_COMPILE);
  glColor4f(1.0, 1.0, 1.0, 1.0);
  for (auto line : polyline) {
    glBegin(GL_LINE_STRIP);
    for (auto p : line) {
      glVertex3d(p.x, p.y, p.z);
    }
    glEnd();
  }
  glEndList();
  return dl;
}

int display_v3d (TV3 p) {
  int dl = glGenLists(1);
  glNewList(dl, GL_COMPILE);
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glBegin(GL_POINTS);
  glVertex3d(p.x, p.y, p.z);
  glEnd();
  glEndList();
  return dl;
}

int display_v2d (TV2 p) {
  int dl = glGenLists(1);
  glNewList(dl, GL_COMPILE);
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glBegin(GL_POINTS);
  glVertex2d(p.x, p.y);
  glEnd();
  glEndList();
  return dl;
}

int compile_geom (std::string expr, bool is_show_lines, bool is_show_normals) {
  // auto star  = star_poly(1.0, 2.0, 4);
  // auto shape = extrude(-8.0, 8.0, star);
  // auto shape = g_mesh(parse_geom("sphere(5.0)"));
  init_cad();
  // printf("EXPR = %s\n", expr.c_str());
  Geom* shape = parse_geom(expr);
  if (is_nested_v2d(shape))
    return display_polyline2(g_nested_v2d_val(shape));
  else if (is_nested_v3d(shape))
    return display_polyline3(g_nested_v3d_val(shape));
  else if (is_poly(shape))
    return display_poly(g_poly_val(shape));
  else if (shape->k == mesh_kind)
    return display_mesh(g_mesh_val(shape), is_show_lines, is_show_normals);
  else if (shape->k == v3d_kind)
    return display_v3d(g_v3d_val(shape));
  else if (shape->k == v2d_kind)
    return display_v2d(g_v2d_val(shape));
  else {
    error("UNDISPLAYABLE GEOM");
    return 0;
  }
}
