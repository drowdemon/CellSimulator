#include <cmath>
#include <cstdlib>
#include "cell.h"
#include "globals.h"

cell::cell(int i)
{
    index=i;
}

void cell::timeStep() //TODO finish
{
    for(unsigned int i=0; i<allMolecules.size(); i++)
    {
        if(allMolecules[i]->position.z+allMolecules[i]->velocity.z>=0 && allMolecules[i]->position.z+allMolecules[i]->velocity.z<volume.size() && allMolecules[i]->position.y+allMolecules[i]->velocity.y>=0 && allMolecules[i]->position.y+allMolecules[i]->velocity.y<volume[(int)(allMolecules[i]->position.z+allMolecules[i]->velocity.z)].size() && allMolecules[i]->position.x+allMolecules[i]->velocity.x>=0 && allMolecules[i]->position.x+allMolecules[i]->velocity.x<volume[(int)(allMolecules[i]->position.z+allMolecules[i]->velocity.z)][(int)(allMolecules[i]->position.y+allMolecules[i]->velocity.y)].size()) //within cell
        {           
            point *vel = &allMolecules[i]->velocity;
            point *pos = &allMolecules[i]->position;
            pointArray velA(*vel);
            pointArray posA(*pos);
            vector<point> *gridPoints = getGridCubes(vel, pos); //get all of the grid cubes the molecule touches in this movement
            //vector<collisionInfo> collisions;
            collisionInfo closestCollision;
            for(unsigned int j=0; j<gridPoints->size(); j++) //search for collisions
            {
                for(unsigned int k=0; k<volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x].size(); k++)
                {
                    if(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x][k]-1==(int)i)
                        continue; //same index... same molecule... do not test for collisions with self. bad idea.
                    //Assumption: all molecules currently in legal positions, no collisions occuring at the start
                    bool good = true;
                    for(int h=0; h<3; h++)
                    {
                        if(abs(posA[h]+velA[h]*EPSILON-((pointArray)(*gridPoints)[j])[h]) > abs(posA[h]-((pointArray)(*gridPoints)[j])[h])) //moving away from this molecule
                        {
                            good = false;
                            break;
                        }
                    }
                    if(!good) //moving away: doesn't intersect
                        continue;
                    //http://stackoverflow.com/questions/2062286/testing-whether-a-line-segment-intersects-a-sphere contains the math behind this. Not very complicated math, but makes this look less like magic
                    pointArray diff(pos->x-allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]-1]->position.x,pos->y-allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]-1]->position.y,pos->z-allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]-1]->position.z);
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
                    if(minDistT>1) // closest point of approach is past the movement, use the endpoint instead
                        minDistT=1;
                    if(minDistT<0)
                        minDistT=0; //shouldn't really ever reach this line. Maybe in bonded (tangent) molecules.
                    double distanceBetweenSquared = sumVelSquared*square(minDistT) - 2*minDistT*sumDiffXVel + sumDiffSquared;
                    if(distanceBetweenSquared <= square(allMolecules[i]->radius+allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]-1]->radius)) //intersects with sphere of other molecule
                    {
                        /*int pos = 0;
                        int max = collisions.size()+1;
                        int min = -1;
                        while(min+1<max)
                        {
                            pos = (max-min)/2 + min;
                            if(minDistT<collisions[pos].t)
                            {
                                max = pos;
                            }
                            else
                            {
                                min = pos;
                            }
                        }
                        collisions.insert(collisions.begin()+pos, collisionInfo(allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]-1], distanceBetweenSquared, minDistT, sumDiffXVel, sumVelSquared, sumDiffSquared)); //instead of just pushing back, sort by minDistT*/
                        if(minDistT<closestCollision.t)
                            closestCollision = collisionInfo(allMolecules[(volume[(int)(*gridPoints)[j].z][(int)(*gridPoints)[j].y][(int)(*gridPoints)[j].x])[k]-1], distanceBetweenSquared, minDistT, sumDiffXVel, sumVelSquared, sumDiffSquared);
                    }
                }
            }
            delete gridPoints;
            //now working with first collision, if it exists
            if(closestCollision.m!=NULL) //if closest collision exists
            {
                double quadraticSqrt = sqrt(4*(square(closestCollision.sumDiffXVel) - closestCollision.sumVelSquared*(closestCollision.sumDiffSquared - closestCollision.distanceAtMinPoint)));
                double answerT = 2*closestCollision.sumDiffXVel;
                if(answerT-quadraticSqrt > 0)
                    answerT-=quadraticSqrt;
                else
                    answerT+=quadraticSqrt;
                answerT/=(2*closestCollision.sumDiffXVel); //at this t, the molecules are tangent
                //Finally, handling the actual collision
                allMolecules[i]->position.x = vel->x*answerT;
                allMolecules[i]->position.y = vel->y*answerT;
                allMolecules[i]->position.z = vel->z*answerT;
                if(!allMolecules[i]->bond(closestCollision.m)) //failed to bond
                    ;//allMolecules[i]->collide(closestCollision.m);
            }
            /*while(abs(toMove.x)>=0.0001 || abs(toMove.x)>=0.0001 || abs(toMove.x)>=0.0001)
            {
                point result(allMolecules[i]->position.x+toMove.x/abs(toMove.x),allMolecules[i]->position.y+toMove.y/abs(toMove.y),allMolecules[i]->position.z+toMove.z/abs(toMove.z));
                if(volume[result.x][result.y][result.z]>0) //collision
                {
                    allMolecules[i]->bond(volume[result.x][result.y][result.z]-1); //may or may not succeed, function checks if bond can be made, this function ensured the position/velocity is correct
                    break;
                }
                else //empty space
                {
                    allMolecules[i]->position=result;
                    toMove.x-=toMove.x/abs(toMove.x);
                    toMove.y-=toMove.y/abs(toMove.y);
                    toMove.z-=toMove.z/abs(toMove.z);
                }
            }*/
        }
        else //collision with membrane
        {
            
        }
    }
}

void cell::solve2dLine(double intercept1, double intercept2, double start, double velocity1, double velocity2, double velocity3, vector<pointArray>& out)
{//example: y=m(x-x0)+y0  y=velocity2/velocity1*(x-start)+intercept; gets y values at integer x coordinates
    if((int)start==(int)(start+velocity1)) //did not leave grid cube in this dimension
        return;
    for(int i=(int)start+(int)(velocity1/abs(velocity1)*EPSILON); i<=(int)start+velocity1; i+=(int)(velocity1/abs(velocity1)*EPSILON))
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

vector<point>* cell::getGridCubes(point* vel, point* pos)
{
    vector<pointArray> intersections[3];
    vector<point> *gridCubes = new vector<point>;
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