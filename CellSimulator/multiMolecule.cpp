#include <cmath>
#include "multiMolecule.h"
#include "globals.h"

multiMolecule::multiMolecule(double m, point v, point cm, int cc)
{
    mass=m;
    velocity=v;
    centerOfMass = cm;
    cellContainer=cc;
}

void multiMolecule::elasticCollision(multiMolecule *m)
{
    point axisOfImpulse(centerOfMass.x - m->centerOfMass.x, centerOfMass.y - m->centerOfMass.y, centerOfMass.z - m->centerOfMass.z);
    /*point VNewAxis = rotate(velocity,axisOfImpulse);
    point VOtherNewAxis = rotate(m->velocity,axisOfImpulse);
    point VNewAxisFinal = VNewAxis;
    
    VNewAxisFinal.x = VELOCITYRETAINEDINCOLLISION*(VNewAxis.x*(mass - m->mass) + 2*m->mass*VOtherNewAxis.x)/(mass + m->mass); //multiplied by constant so as to lose some energy upon collision
    VOtherNewAxis.x = VELOCITYRETAINEDINCOLLISION*(VOtherNewAxis.x*(m->mass - mass) + 2*mass*VNewAxis.x)/(mass + m->mass);
    velocity = unRotate(VNewAxisFinal, axisOfImpulse);
    m->velocity = unRotate(VOtherNewAxis, axisOfImpulse);*/
    double length = sqrt(square(axisOfImpulse.x) + square(axisOfImpulse.y) + square(axisOfImpulse.z));
    point normAxis(axisOfImpulse.x/length, axisOfImpulse.y/length, axisOfImpulse.z/length);
    double parallelComponent = normAxis.x*velocity.x+normAxis.y*velocity.y+normAxis.z*velocity.z;
    point velocityParallel(normAxis.x*parallelComponent, normAxis.y*parallelComponent, normAxis.z*parallelComponent);
    point velocityPerp(velocity.x-velocityParallel.x, velocity.y-velocityParallel.y, velocity.z-velocityParallel.z);
    
    double parallelComponentOther = normAxis.x*m->velocity.x+normAxis.y*m->velocity.y+normAxis.z*m->velocity.z;
    point velocityParallelOther(normAxis.x*parallelComponentOther, normAxis.y*parallelComponentOther, normAxis.z*parallelComponentOther);
    point velocityPerpOther(m->velocity.x-velocityParallelOther.x, m->velocity.y-velocityParallelOther.y, m->velocity.z-velocityParallelOther.z);
    
    velocity.x = VELOCITYRETAINEDINCOLLISION*((velocityParallel.x*(mass - m->mass) + 2*m->mass*velocityParallelOther.x)/(mass + m->mass)+velocityPerp.x); //multiplied by constant so as to lose some energy upon collision
    m->velocity.x = VELOCITYRETAINEDINCOLLISION*(velocityParallelOther.x*(m->mass - mass) + 2*mass*velocityParallel.x)/(mass + m->mass);
    velocity.y = VELOCITYRETAINEDINCOLLISION*((velocityParallel.y*(mass - m->mass) + 2*m->mass*velocityParallelOther.y)/(mass + m->mass)+velocityPerp.y);
    m->velocity.y = VELOCITYRETAINEDINCOLLISION*(velocityParallelOther.y*(m->mass - mass) + 2*mass*velocityParallel.y)/(mass + m->mass);
    velocity.z = VELOCITYRETAINEDINCOLLISION*((velocityParallel.z*(mass - m->mass) + 2*m->mass*velocityParallelOther.z)/(mass + m->mass)+velocityPerp.z);
    m->velocity.z = VELOCITYRETAINEDINCOLLISION*(velocityParallelOther.z*(m->mass - mass) + 2*mass*velocityParallel.z)/(mass + m->mass);
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
void multiMolecule::inElasticCollision(point otherVelocity, double otherMass, double myMass) //so much easier
{    
    velocity.x = (velocity.x*myMass + otherMass*otherVelocity.x)/(myMass + otherMass);
    velocity.y = (velocity.y*myMass + otherMass*otherVelocity.y)/(myMass + otherMass);
    velocity.z = (velocity.z*myMass + otherMass*otherVelocity.z)/(myMass + otherMass);
}

void multiMolecule::transport(point &displacement)
{
    centerOfMass.x += displacement.x;
    centerOfMass.y += displacement.y;
    centerOfMass.z += displacement.z;
    for(unsigned int j=0; j<molecules.size(); j++)
    {
        point p(allCells[cellContainer].allMolecules[molecules[j]]->position);
        p.x += displacement.x;
        p.y += displacement.y;
        p.z += displacement.z;
        allCells[cellContainer].allMolecules[molecules[j]]->teleport(p);
    }
}