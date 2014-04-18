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

double& point::operator [](int i)
{
    if(i==0)
        return x;
    else if(i==1)
        return y;
    else if(i==2)
        return z;
}

pointArray::pointArray(double px, double py, double pz)
{
    p[0]=px;
    p[1]=py;
    p[2]=pz;
}

pointArray::pointArray(point from)
{
    *this=pointArray(from.x, from.y, from.z);
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

collisionInfo::collisionInfo(molecule* pm, double pt, molecule *wc)
{
    m = pm;
    t = pt;
    whatCollided = wc;
}

collisionInfo::collisionInfo()
{
    m = whatCollided = NULL;
    t=9999999;
}
