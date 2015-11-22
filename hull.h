#ifndef __IS_QUICK_HULL__
#define __IS_QUICK_HULL__

#include <vector>
#include <iostream>
#include <geode/vector/vector3d.h>
#include "cad.h"

using namespace geode;

extern Mesh quick_hull_mesh (Mesh mesh);
extern Nested<TV2> quick_hull_poly (Nested<TV2> poly);

#endif
