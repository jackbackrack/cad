#ifndef __APP__
#define __APP__

#include "gui.h"
#include "gfx.h"
#include "viz.h"
#include "vec.h"
#include "lay.h"
#include "lisp.h"
#include "reader.h"

class sim_viz_t : public viz_t {
 public:
  int key_hit (int cmd, int modifiers);
  int exec (int is_pause);
  int render (int is_picking);
  int render_one(int is_picking);
  int render_frame_monitors ( void );
  int handle_drag(vec_t<3> pos);
  int process_picks (std::vector< int > picks);
  int open (int arg_offset, int argc, const char *argv[]);
  int close (void);
  bool is_single_step;

  void show_status (float x, float y, float w, float h);
};

class sim_t {
 public:
  int t;
  int list;
  void exec ( void );
  void open ( void );
  void render ( bool is_picking );
 sim_t(void) : t(0),list(-1) { }
};

#endif
