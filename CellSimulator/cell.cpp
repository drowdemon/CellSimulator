#include <cmath>
#include <cstdlib>
#include <iostream>
#include "cell.h"
#include "globals.h"
#include "multiMolecule.h"

cell::cell(int i)
{
    index=i;
}

void cell::timeStep() //TODO finish
{
    for(unsigned int i=0; i<allMultiMolecules.size(); i++)
    {
        if(!allMultiMolecules[i])
            continue;
        if(abs(allMultiMolecules[i]->velocity.x)<EPSILON && abs(allMultiMolecules[i]->velocity.y)<EPSILON && abs(allMultiMolecules[i]->velocity.z)<EPSILON) //velocity is 0
            continue; //not moving! thank god.
        checkCollisionAndMove(i, 1.0);
    }
}

void cell::solve2dLine(double intercept1, double intercept2, double start, double velocity1, double velocity2, double velocity3, vector<pointArray>& out)
{//example: y=m(x-x0)+y0  y=velocity2/velocity1*(x-start)+intercept; gets y values at integer x coordinates
    if((int)start==(int)(start+velocity1)) //did not leave grid cube in this dimension
        return;
    for(int i=(int)(start/*+(velocity1/abs(velocity1)*EPSILON)*/); ((i-(int)(start+velocity1))==0) || (((i-(int)(start+velocity1))>0) == (-velocity1>0)); i+=(int)(velocity1/abs(velocity1))) //starts where it does to prevent i=.9999999 is i=0 errors. I think. //The code the previous comment refers to has been commented out.
    {
        out.push_back(pointArray(i,velocity2/velocity1*(i-start)+intercept1,velocity3/velocity1*(i-start)+intercept2));
    }
}

bool cell::addGridCubes(unsigned int counter, int flagWhere, pointArray vel, pointArray curLoc, vector<point> &gridCubes)
{
    while(flagWhere<3 || (counter&(1<<flagWhere))==0)
    {
        flagWhere++;
    }
    if(flagWhere>=3)
        return false;
    if(!addGridCubes(counter, flagWhere, vel, curLoc, gridCubes)) //need to actually add to gridCubes, base case
    {
        gridCubes.push_back(curLoc);
        curLoc[flagWhere]+=vel[flagWhere]/abs(vel[flagWhere]*EPSILON);
        gridCubes.push_back(curLoc);
    }
    else //recursion will/has handle(d) it, just recurse into the next level
    {
        curLoc[flagWhere]+=vel[flagWhere]/abs(vel[flagWhere]*EPSILON);
        addGridCubes(counter, flagWhere, vel, curLoc, gridCubes);
    }
    return true;
}

void cell::getGridCubesHelper(point* vel, point* pos, vector<point> *gridCubes)
{
    vector<pointArray> intersections[3];
    solve2dLine(pos->y, pos->z, pos->x, vel->x, vel->y, vel->z, intersections[0]); //integer x: interesect yz plane
    solve2dLine(pos->x, pos->z, pos->y, vel->y, vel->x, vel->z, intersections[1]); //integer y: xz plane
    for(unsigned int i=0; i<intersections[1].size(); i++)
    {
        double temp = intersections[1][i][0];
        intersections[1][i][0] = intersections[1][i][1];
        intersections[1][i][1] = temp;
    }
    solve2dLine(pos->y, pos->x, pos->z, vel->z, vel->y, vel->x, intersections[2]); //integer z: xy plane
    for(unsigned int i=0; i<intersections[2].size(); i++)
    {
        double temp = intersections[2][i][0];
        intersections[2][i][0] = intersections[2][i][2];
        intersections[2][i][2] = temp;
    }
    gridCubes->push_back(*pos);
    for(int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<intersections[i].size(); j++)
        {
            unsigned int counter=0;
            for(int k=0; k<3; k++)
            {
                if(abs(intersections[i][j][k]-(int)intersections[i][j][k])<EPSILON)
                    counter |= 1<<k;
            }
            if(counter!=0)
                addGridCubes(counter, 0, *vel, *pos, *gridCubes);
        }
    }
}

vector<point>* cell::getGridCubes(point *vel, point *pos, molecule *m)
{
    vector<point> *gridCubes = new vector<point>;
    vector<point> temp;
    m->getTouchingPoints(temp, *pos);
    for(unsigned int i=0; i<temp.size(); i++)
    {
        getGridCubesHelper(vel, &temp[i], gridCubes);
    }
    for(unsigned int i=0; i<gridCubes->size(); i++) //remove duplicates
    {
        for(unsigned int j=i+1; j<gridCubes->size(); j++)
        {
            if((int)(*gridCubes)[i].x==(int)(*gridCubes)[j].x && (int)(*gridCubes)[i].y==(int)(*gridCubes)[j].y && (int)(*gridCubes)[i].z==(int)(*gridCubes)[j].z) //same grid cube
            {
                gridCubes->erase(gridCubes->begin()+j);
                j--;
            }
        }
    }
    return gridCubes;
}

void cell::checkCollisionAndMove(unsigned int i, double timeRemaining, int prevCollidedWith)
{
    collisionInfo closestCollision;
    point vel = allMultiMolecules[i]->velocity;
    vel.x*=timeRemaining;
    vel.y*=timeRemaining;
    vel.z*=timeRemaining;
    pointArray velA(vel);
    for(unsigned int a=0; a<allMultiMolecules[i]->molecules.size(); a++)
    {
        if(allMolecules[allMultiMolecules[i]->molecules[a]]->position.z+vel.z>=0 && allMolecules[allMultiMolecules[i]->molecules[a]]->position.z+vel.z<volume.size() && allMolecules[allMultiMolecules[i]->molecules[a]]->position.y+vel.y>=0 && allMolecules[allMultiMolecules[i]->molecules[a]]->position.y+vel.y<volume[(int)(allMolecules[allMultiMolecules[i]->molecules[a]]->position.z+vel.z)].size() && allMolecules[allMultiMolecules[i]->molecules[a]]->position.x+vel.x>=0 && allMolecules[allMultiMolecules[i]->molecules[a]]->position.x+vel.x<volume[(int)(allMolecules[allMultiMolecules[i]->molecules[a]]->position.z+vel.z)][(int)(allMolecules[allMultiMolecules[i]->molecules[a]]->position.y+vel.y)].size()) //within cell
        {           
            point *pos = &allMolecules[allMultiMolecules[i]->molecules[a]]->position;
            pointArray posA(*pos);
            vector<point> *gridPoints = getGridCubes(&vel, pos, allMolecules[allMultiMolecules[i]->molecules[a]]); //get all of the grid cubes the molecule touches in this movement
            //vector<collisionInfo> collisions;
            set<int> moleculesTested;
            for(unsigned int j=0; j<gridPoints->size(); j++) //search for collisions
            {
                for(unsigned int k=0; k<volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x].size(); k++)
                {
                    if(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x].size()>1)
                        cout << "h" << endl;
                    int testingIndex = volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x][k];
                    bool dead = false;
                    for(unsigned int h=0; h<allMultiMolecules[i]->molecules.size(); h++)
                    {
                        if(testingIndex==allMolecules[allMultiMolecules[i]->molecules[h]]->index)
                        {
                             dead=true; //same index as someone in my multimolecule... same multimolecule... do not test for collisions with self. Multimolecule is allowed to collide with itself within a timestep, since all molecules in it have the same velocity.
                             break;
                        }   
                    }
                    if(dead)
                        continue;
                    if(moleculesTested.find(testingIndex)!=moleculesTested.end()) //already checked this in another gridsquare
                        continue;
                    moleculesTested.insert(testingIndex);
                    //Assumption: all molecules currently in legal positions, no collisions occuring at the start
                    bool good = true;
                    for(int h=0; h<3; h++)
                    {
                        if(abs(posA[h]+velA[h]*EPSILON-((pointArray)allMolecules[testingIndex]->position)[h]) > abs(posA[h]-((pointArray)allMolecules[testingIndex]->position)[h])) //moving away from this molecule
                        {
                            good = false;
                            break;
                        }
                    }
                    if(!good) //moving away: doesn't intersect
                        continue;
                    //http://stackoverflow.com/questions/2062286/testing-whether-a-line-segment-intersects-a-sphere contains the math behind this. Not very complicated math, but makes this look less like magic
                    pointArray diff(allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]]->position.x-pos->x,allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]]->position.y-pos->y,allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]]->position.z-pos->z);
                    double minDistT = 0;
                    double sumDiffXVel = 0;
                    double sumVelSquared = 0;
                    double sumDiffSquared =0;
                    for(int h=0; h<3; h++)
                    {
                        sumDiffXVel += diff[h]*velA[h];
                        sumVelSquared += square(velA[h]);
                        sumDiffSquared += square(diff[h]);
                    }
                    minDistT = sumDiffXVel / sumVelSquared;
                    if(minDistT>timeRemaining) // closest point of approach is past the movement, use the endpoint instead
                        minDistT=timeRemaining;
                    if(minDistT<0)
                    {
                        minDistT=0; //shouldn't really ever reach this line.
                        cout << "Huh? MinDistT < 0. -12" << endl;
                    }
                    double distanceBetweenSquared = sumVelSquared*square(minDistT) - 2*minDistT*sumDiffXVel + sumDiffSquared;
                    double neededDistSquared = square(allMolecules[allMultiMolecules[i]->molecules[a]]->radius+allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]]->radius);
                    if(distanceBetweenSquared <= neededDistSquared) //intersects with sphere of other molecule
                    {
                        //solve the quadratic for when the molecules just barely touch
                        double quadraticSqrt = sqrt(4*(square(sumDiffXVel) - sumVelSquared*(sumDiffSquared - neededDistSquared)));
                        double answerT = 2*sumDiffXVel;
                        if(answerT-quadraticSqrt > 0)
                            answerT-=quadraticSqrt;
                        else
                            answerT+=quadraticSqrt;
                        answerT/=(2*sumVelSquared); //at this t, the molecules are tangent
                        if(answerT<0 || answerT>timeRemaining)
                        {
                            cout << "What the hell?? answerT out of bounds" << endl;
                            exit(-11);
                        }
                        if(answerT<closestCollision.t)
                            closestCollision = collisionInfo(allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]], answerT, allMolecules[allMultiMolecules[i]->molecules[a]]);
                    }
                }
            }
            delete gridPoints;
        }
        else //collision with membrane
        {
            //TODO code here
        }
    }
    //now working with first collision, if it exists
    if(closestCollision.m!=NULL) //if closest collision exists
    {
        //Finally, handling the actual collision
        point otherV(allMultiMolecules[closestCollision.m->indexMultiMolecule]->velocity); // for inelastic
        double otherM = allMultiMolecules[closestCollision.m->indexMultiMolecule]->mass; // for inelastic
        double myMass = allMultiMolecules[i]->mass;
        point disp = vel;
        disp.x *= closestCollision.t;
        disp.y *= closestCollision.t;
        disp.z *= closestCollision.t;
        allMultiMolecules[i]->transport(disp);
        int retBond = closestCollision.whatCollided->bond(closestCollision.m);
        if(retBond!=-1) //bonded successfully, inelastic collision
        {
            allMultiMolecules[retBond]->inElasticCollision(otherV, otherM, myMass);
            for(int j=0; j<3; j++)
            {
                if(((pointArray)(allMultiMolecules[retBond]->velocity))[j]<EPSILON)
                    ((pointArray)(allMultiMolecules[retBond]->velocity))[j] = 0;
            }
            checkCollisionAndMove(retBond, timeRemaining-closestCollision.t);
        }
        else
        {
            if(closestCollision.t<EPSILON && prevCollidedWith==closestCollision.m->indexMultiMolecule) //no actual movement could occur
            { //do this only if I already collided with this molecule in this frame
                cout << "collided twice" << endl;
                return; //just stop moving, allow the other molecule to move first
            }
            allMultiMolecules[i]->elasticCollision(allMultiMolecules[closestCollision.m->indexMultiMolecule]);
            for(int j=0; j<3; j++)
            {
                if(abs(allMultiMolecules[i]->velocity[j])<EPSILON)
                    allMultiMolecules[i]->velocity[j] = 0;
                if(abs(allMultiMolecules[closestCollision.m->indexMultiMolecule]->velocity[j])<EPSILON)
                    allMultiMolecules[closestCollision.m->indexMultiMolecule]->velocity[j] = 0;
            }
            checkCollisionAndMove(i, timeRemaining-closestCollision.t, closestCollision.m->indexMultiMolecule);
        }
    }
    else //just move
    {
        allMultiMolecules[i]->transport(vel);
    }
}