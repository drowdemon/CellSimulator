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
    X(point, velocity) \
    X(double, radius) \
    X(double, mass)

class molecule
{
public:
    int index;
    point position;
    unsigned long long id;
    int cellContainer;
    int moleculeType; //0 = protein, 1 = DNA, 2 = RNA
    point velocity;
    double radius; //all molecules are spheres. Yes, DNA is spherical.
    double mass;
    
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
    bool bond(int indexWithWhat);
    bool bond(molecule *m);
};

#endif	/* MOLECULE_H */

