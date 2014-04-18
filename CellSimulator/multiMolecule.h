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
    int cellContainer;
    multiMolecule(double m, point v, point cm, int cc);
    void elasticCollision(multiMolecule *m);
    void inElasticCollision(point otherVelocity, double otherMass, double myMass);
    point rotate(point v, point &rotateToAxes);
    point unRotate(point v, point &rotateToAxes);
    void transport(point &displacement);
};

#endif	/* MULTIMOLECULE_H */

