#include <cmath>
#include "multiMolecule.h"

multiMolecule::multiMolecule(double m, point v, point cm)
{
    mass=m;
    velocity=v;
    centerOfMass = cm;
}

void multiMolecule::collision(multiMolecule *m)
{
    point axisOfImpulse(centerOfMass.x - m->centerOfMass.x, centerOfMass.y - m->centerOfMass.y, centerOfMass.z - m->centerOfMass.z);
    point VNewAxis = rotate(velocity,axisOfImpulse);
    point VOtherNewAxis = rotate(m->velocity,axisOfImpulse);
    point VNewAxisFinal = VNewAxis;
    
    VNewAxisFinal.x = (VNewAxis.x*(mass - m->mass) + 2*m->mass*VOtherNewAxis.x)/(mass + m->mass);
    VOtherNewAxis.x = (VOtherNewAxis.x*(m->mass - mass) + 2*mass*VNewAxis.x)/(mass + m->mass);
    velocity = unRotate(VNewAxisFinal, axisOfImpulse);
    m->velocity = unRotate(VOtherNewAxis, axisOfImpulse);
}

point multiMolecule::rotate(point v, point &rotateToAxes)
{
    point ret=v;
    //rotate in xy plane
    double angleA = -atan2(rotateToAxes.y, rotateToAxes.x);
    ret.x = v.x*cos(angleA) - v.y*sin(angleA);
    ret.y = v.x*sin(angleA) + v.y*cos(angleA);
    //rotate in xz plane
    angleA = -atan2(rotateToAxes.z, rotateToAxes.x);
    v.x = ret.x*cos(angleA) - ret.z*sin(angleA);
    v.z = ret.x*sin(angleA) + ret.z*cos(angleA);
    v.y=ret.y;
    
    return v;
}

point multiMolecule::unRotate(point v, point &rotateToAxes)
{
    point ret=v;
    //rotate in xz plane
    double angleA = atan2(rotateToAxes.z, rotateToAxes.x);
    ret.x = v.x*cos(angleA) - v.z*sin(angleA);
    ret.z = v.x*sin(angleA) + v.z*cos(angleA);
    
    //rotate in xy plane
    angleA = atan2(rotateToAxes.y, rotateToAxes.x);
    v.x = ret.x*cos(angleA) - ret.y*sin(angleA);
    v.y = ret.x*sin(angleA) + ret.y*cos(angleA);
    v.z=ret.z;
    return v;
}   