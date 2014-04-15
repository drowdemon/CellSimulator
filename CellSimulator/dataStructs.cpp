#include "dataStructs.h"
#include <cstdlib>

point::point(double px, double py, double pz)
{
    x=px;
    y=py;
    z=pz;
}

point::point()
{
    x=y=z=0;
}

point::point(pointArray& p)
{
    x=p[0];
    y=p[1];
    z=p[2];
}
pointArray::pointArray(double px, double py, double pz)
{
    p[0]=px;
    p[1]=py;
    p[2]=pz;
}

pointArray::pointArray(point from)
{
    pointArray(from.x, from.y, from.z);
}

pointArray::pointArray()
{
    p[0]=p[1]=p[2]=0;
}

pointArray::pointArray(double from[])
{
    p[0]=from[0];
    p[1]=from[1];
    p[2]=from[2];
}

pointArray::pointArray(const pointArray& from)
{
    p[0]=from.p[0];
    p[1]=from.p[1];
    p[2]=from.p[2];
}

double& pointArray::operator [](int i)
{
    return p[i];
}

collisionInfo::collisionInfo(molecule* pm, double dmp, double pt, double sDxV, double sVS, double sDS)
{
    m = pm;
    distanceAtMinPoint = dmp;
    t = pt;
    sumDiffXVel = sDxV;
    sumVelSquared = sVS;
    sumDiffSquared = sDS;
}

collisionInfo::collisionInfo()
{
    m = NULL;
    distanceAtMinPoint = t = sumDiffSquared = sumDiffXVel = sumVelSquared = 0;
}
