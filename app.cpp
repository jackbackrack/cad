#include "app.h"
#include "read-eval.h"

sim_t* sim = NULL;

extern std::string user_msg;

props_t* props;

defboolprop(is_show_status, false);
defboolprop(is_show_lines, true);
defboolprop(is_show_normals, false);
defboolprop(is_interpolating, false);
defboolprop(is_timing, false);
defnumpropmod(flo, radius,     8,   0.1,  0.0, 1000.0); 
defnumpropmod(flo, threshold,  0.1, 0.01, 0.0, 1.0); 
defnumpropmod(int, num_points, 10000, 100, 1, 1000000000); 
defstrprop(expr,  "cube(4.0) - xmov(2.5,cube(2.0))"); 

int sim_viz_t::exec (int is_pause) {
  sim->exec();
  return 1;
}

int sim_viz_t::handle_drag (vec_t<3> pos) {
  printf("DRAG %f,%f,%f\n", pos.x, pos.y, pos.z);
  return 1;
}

int sim_viz_t::process_picks (std::vector< int > picks) {
  return 1;
}
keys_t* top_keys;
cmds_t* top_cmds = new cmds_t();
cmds_t* cmds;

void install_keys (void) {
  now_keys = top_keys = new keys_t(key_not_found);
  top_keys->install("z", reset_view_cmd);
  top_keys->install("f", toggle_full_screen_cmd);
  top_keys->install("q", quit_cmd);
  top_keys->install("S-U", new_toggle_prop_cmd(is_show_status_var));
}

int sim_viz_t::key_hit (int cmd, int modifiers) {
  user_msg.clear();
  is_key_hit[cmd] = 1;
  key_modifiers = modifiers;
  now_keys->do_process_keys(cmd, 0, key_modifiers, sim);
  return 1;
}

void sim_viz_t::show_status (float x, float y, float w, float h) {
  char text[1000];
  glPushMatrix(); glPushAttrib(GL_CURRENT_BIT);
  glColor3f(1, 0, 1);
  char mode[10000];
  mode[0] = 0;
  sprintf(mode, "%d", sim->t);
  glTranslatef( x - w/4, y, 5); 
  // draw_text(100, 10, mode);
  draw_text(w, h, mode);
  // draw_text(w, h, text);
  glPopAttrib(); glPopMatrix();
}

int sim_viz_t::render_frame_monitors ( void ) {
  if (is_show_status) {
    show_status(viz->MAX_X-15, viz->MIN_Y+2, 30, 20);
  }
  return 1;
}

int sim_viz_t::render (int is_picking) {
  sim->render(is_picking);
  return 1;
}

viz_t* new_viz (void) {
  viz_t* viz = new sim_viz_t();
  return viz;
}

static int install_props (void) {
  props = new props_t();
  props->install(is_show_status_var);
  props->install(is_timing_var);
  props->install(is_interpolating_var);
  props->install(is_show_lines_var);
  props->install(is_show_normals_var);
  props->install(radius_var);
  props->install(threshold_var);
  props->install(num_points_var);
  props->install(expr_var);
  return 1;
}

static void sim_parse_args (int argc, const char *argv[]) {
  std::vector<const char*> args;
  for (int i = 1; i < (argc-1); i++) {
    args.push_back(argv[i]);
  }
  if (argc > 1) {
    args.push_back(":expr");
    args.push_back(argv[argc-1]);
  }
  std::vector<const char*>::iterator ap  = args.begin();
  std::vector<const char*>::iterator eap = args.end();
  props->parse_args(ap, eap, sim);
}

void sim_t::open ( void ) {
  printf("STR = %s\n", expr.c_str());
}

void sim_t::render ( bool is_picking ) {
  if (list != -1) {
    glPushMatrix();
    glScalef(8.0,8.0,8.0);
    glCallList(list);
    glPopMatrix();
  }
}

void sim_t::exec ( void ) {
  if (list == -1) {
    list = compile_geom(expr, is_show_lines, is_show_normals);
  }
  t += 1;
}

int sim_viz_t::open(int args_offset, int argc, const char** argv) {
  install_props();
  install_keys();
  sim_parse_args(argc, argv);
  sim = new sim_t();
  sim->open();
  return 1;
}

int sim_viz_t::close ( void ) {
  // sim->close();
  return 0;
}

