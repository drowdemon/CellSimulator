#ifndef MULTIMOLECULE_H
#define	MULTIMOLECULE_H

#include <vector>
#include "dataStructs.h"

using namespace std;

class multiMolecule
{
public:
    vector<int> molecules;
    double mass;
    point velocity;
    point centerOfMass;
    multiMolecule(double m, point v, point cm);
    void collision(multiMolecule *m);
    point rotate(point v, point &rotateToAxes);
    point unRotate(point v, point &rotateToAxes);
};

#endif	/* MULTIMOLECULE_H */

