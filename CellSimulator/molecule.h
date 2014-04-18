#ifndef MOLECULE_H
#define	MOLECULE_H

#include <vector>
#include <set>
#include "dataStructs.h"

using namespace std;

#define TYPE_PROTEIN 0
#define TYPE_DNA 1
#define TYPE_RNA 2

#define LISTVARMOLECULE \
    X(int, index) \
    X(point, position) \
    X(unsigned long long, id) \
    X(int, cellContainer) \
    X(int, moleculeType) \
    X(double, radius) \
    X(double, mass) \
    X(int, indexMultiMolecule)

class molecule
{
public:
    int index; //in allMolecules in the containing cell
    point position;
    unsigned long long id;
    int cellContainer; //what cell contains it
    int moleculeType; //0 = protein, 1 = DNA, 2 = RNA
    double radius; //all molecules are spheres. Yes, DNA is spherical.
    double mass;
    int indexMultiMolecule; //index of the containing multiMolecule
    
    vector<int> bondedWith; //indexes of bonds
    vector<set<unsigned long long> > bondsNaturallyWith; //a list of possible bonds. Each one represents one possible bond. The length of this is the total number of bonds there can be.
    vector<bool> usedBond; //matches size of bondsNaturallyWith, says whether that bond was filled.
#define X(type, name) \
    type p_ ## name,
    
    molecule(LISTVARMOLECULE bool extra);
#undef X
    virtual unsigned long long getID() = 0;
protected:
    vector<molecule*>* checkAround(double searchRadius);
public:
    int bond(molecule *m);
    void teleport(point p); //instantly teleport to the location p
    void getTouchingPointsHelper(vector<point> &p, pointArray center, int axis);
    void getTouchingPoints(vector<point> &p, point center);
};

#endif	/* MOLECULE_H */

