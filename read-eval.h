#ifndef __READ_EVAL__
#define __READ_EVAL__

#include <cstdio>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

typedef int TokenType;
const int tok_eof = -1;
const int tok_sym = -2;
const int tok_str = -3;
const int tok_num = -4;
const int tok_null = -5;

class Token {
 public:
  TokenType   type;
  std::string sym;
  double      num;
  void clear() { type = tok_null; sym.clear(); }
  Token () : type(tok_null) { }
  Token (TokenType typeArg) : type(typeArg) { }
  Token (std::string symArg, bool is_string) : type(is_string ? tok_str : tok_sym), sym(symArg) { }
  Token (double numArg) : type(tok_num), num(numArg) { }
  std::string to_str(void) { 
    if (type == tok_eof)
      return "EOF";
    else if (type == tok_sym)
      return sym;
    else if (type == tok_str)
      return "\"" + sym + "\"";
    else if (type == tok_num)
      return std::to_string(num);
    else if (type == tok_null)
      return "NULL";
    else {
      std::string s = { 'T', '=', (char)type };
      return s;
    }
  }
};

class Tokenizer {
 public:
  Token token;
  int   line_number;
  std::istream* s;
 Tokenizer( void ) : line_number(0), s(NULL) { }
 Tokenizer(std::istream* sArg, int line_number) : line_number(line_number), s(sArg) { }
  Token peek();
  Token get();
  Token unget(Token& tok) { return token = tok; };
  Token do_get();
  double get_double(void);
  int get_int(void);
  char get_char(void);
  std::string get_sym(void);
  void expect(char kind);
  void expect(std::string name);
};

typedef double flo_t;

inline flo_t sqr (flo_t x) { return x * x; }

extern int compile_geom (std::string expr, bool is_show_lines, bool is_show_normals);

extern int display_triangles_list (std::vector<Tri> &tris, bool is_lines = true, bool is_show_normals = false);

#endif
