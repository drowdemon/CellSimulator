#ifndef GLOBALS_H
#define	GLOBALS_H

#include <vector>
#include "cell.h"

using namespace std;

#define EPSILON 0.0001
#define PROTEINCANFUNCTIONRADIUSFRACTION 1/4

double square(double val);

extern vector<cell> allCells;

#endif	/* GLOBALS_H */

