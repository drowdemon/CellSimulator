#ifndef DATASTRUCTS_H
#define	DATASTRUCTS_H

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
    molecule* m;
    double distanceAtMinPoint;
    double t;
    double sumDiffXVel;
    double sumVelSquared;
    double sumDiffSquared;
    collisionInfo(molecule *pm, double dmp, double pt, double sDxV, double sVS, double sDS);
    collisionInfo();
};
#endif	/* DATASTRUCTS_H */