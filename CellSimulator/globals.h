#ifndef GLOBALS_H
#define	GLOBALS_H

#include <vector>
#include "cell.h"

using namespace std;

#define EPSILON 0.00001
#define PROTEINCANFUNCTIONRADIUSFRACTION 1/4
#define VELOCITYRETAINEDINCOLLISION 0.99

double square(double val);

extern vector<cell> allCells;

#endif	/* GLOBALS_H */

