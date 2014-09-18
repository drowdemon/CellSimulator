#ifndef DATASTRUCTS_H
#define	DATASTRUCTS_H

#include <set>

using namespace std;

class pointArray;
class molecule;

class point
{
public:
    double x;
    double y;
    double z;
    point(double px, double py, double pz);
    point();
    point(pointArray &p);
    double& operator[](int i);
};

class pointArray
{
public:
    double p[3];
    pointArray(double px, double py, double pz);
    pointArray(point from);
    pointArray();
    pointArray(double from[3]);
    pointArray(const pointArray &from);
    double& operator[](int i);
};

class collisionInfo
{
public:
    molecule* m; //collision with this molecule
    double t;
    molecule *whatCollided; //collision by this molecule - this one was the one in the loop, trying to move, and couldny
    collisionInfo(molecule *pm, double pt, molecule *wc);
    collisionInfo();
};

class bondSpot
{
public:
    set<unsigned long long> bondsNaturallyWith;
    bool inUse;
    point pointBond; //distance from center
    double radiusBond;
    bondSpot(set<unsigned long long> bnw, bool iu, point pb, double rb);
};

class bondData
{
public:
    int with;
    point contact; //vector from molecule center
    bondData(int w, point c);
};
#endif	/* DATASTRUCTS_H */