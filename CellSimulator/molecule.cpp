#include <cstdlib>
#include <cmath>

#include "molecule.h"
#include "globals.h"
#include "multiMolecule.h"

molecule::molecule(int p_index, point p_position, unsigned long long p_id, int p_cellContainer, int p_moleculeType, double p_radius, double p_mass, int p_indexMultiMolecule, bool extra)
{
#define X(type, name) \
    name = p_ ## name;
    LISTVARMOLECULE
#undef X
}

vector<molecule*>* molecule::checkAround(double searchRadius) //return molecules whose volume intersects with searchRadius
{
    vector<molecule*>* ret = new vector<molecule*>;
    for(int i=(((int)floor(position.z-searchRadius)!=(int)floor(position.z)) ? -(int)floor(searchRadius) : 0); i<=(((int)floor(position.z+searchRadius)!=(int)floor(position.z)) ? (int)ceil(searchRadius) : 0); i++)
    {
        for(int j=(((int)floor(position.y-searchRadius)!=(int)floor(position.y)) ? -(int)floor(searchRadius) : 0); j<=(((int)floor(position.y+searchRadius)!=(int)floor(position.y)) ? (int)ceil(searchRadius) : 0); j++)
        {
            for(int k=(((int)floor(position.x-searchRadius)!=(int)floor(position.x)) ? -(int)floor(searchRadius) : 0); k<=(((int)floor(position.x+searchRadius)!=(int)floor(position.x)) ? (int)ceil(searchRadius) : 0); k++)
            {
                if(position.z+i>=0 && position.z+i<allCells[cellContainer].volume.size() && position.y+j>=0 && position.y+j<allCells[cellContainer].volume[(int)(position.z+i)].size() && position.x+k>=0 && position.x+k<allCells[cellContainer].volume[(int)(position.z+i)][(int)(position.y+j)].size()) //within cell
                {
                    if(allCells[cellContainer].volume[(int)position.z+i][(int)position.y+j][(int)position.x+k].size()>0)
                    {
                        vector<int>* tempCells = &allCells[cellContainer].volume[(int)position.z+i][(int)position.y+j][(int)position.x+k];
                        for(unsigned int h=0; h<allCells[cellContainer].volume[(int)position.z+i][(int)position.y+j][(int)position.x+k].size(); h++)
                        {
                            if(square(searchRadius+allCells[cellContainer].allMolecules[(*tempCells)[h]]->radius) < square(allCells[cellContainer].allMolecules[(*tempCells)[h]]->position.x-position.x) + square(allCells[cellContainer].allMolecules[(*tempCells)[h]]->position.y-position.y) + square(allCells[cellContainer].allMolecules[(*tempCells)[h]]->position.z-position.z)) //intersect
                            {
                                ret->push_back(allCells[cellContainer].allMolecules[(*tempCells)[h]]);
                            }
                        }
                    }
                }
            }
        }
    }
    return ret;
}
int molecule::bond(molecule *m, point contact) //make sure this isn't called when it shouldn't be. Use with great care. It should only be called when molecules actually collide, or when an enzymatic protein bonds things together. It makes no checks on the legality of the bond based on position/velocity, only on type of bond
{
    for(unsigned int i=0; i<possibleBonds.size(); i++)
    {
        if(!possibleBonds[i].inUse && possibleBonds[i].bondsNaturallyWith.find(m->id)!=possibleBonds[i].bondsNaturallyWith.end()) //can bond with this type of molecule, and didn't already
        {
            if(square(contact.x-(position.x+possibleBonds[i].pointBond.x))+square(contact.y-(position.y+possibleBonds[i].pointBond.y))+square(contact.z-(position.z+possibleBonds[i].pointBond.z)) < square(possibleBonds[i].radiusBond)) //within area that can be bonded to
            {
                for(unsigned int j=0; j<m->possibleBonds.size(); j++)
                {
                    if(!m->possibleBonds[j].inUse && m->possibleBonds[j].bondsNaturallyWith.find(id)!=m->possibleBonds[j].bondsNaturallyWith.end()) //the molecule I'm bonding with has an empty slot
                    {
                        if(square(contact.x-(m->position.x+m->possibleBonds[j].pointBond.x))+square(contact.y-(m->position.y+m->possibleBonds[j].pointBond.y))+square(contact.z-(m->position.z+m->possibleBonds[j].pointBond.z)) < square(m->possibleBonds[j].radiusBond)) //within area that can be bonded to
                        {
                            //start bonding
                            bondedWith.push_back(bondData(m->index, point(contact.x-position.x, contact.y-position.y, contact.z-position.z)));
                            m->bondedWith.push_back(bondData(index, point(contact.x-m->position.x, contact.y-m->position.y, contact.z-m->position.z)));
                            possibleBonds[i].inUse=true;
                            m->possibleBonds[i].inUse=true;
                            molecule *smaller = m;
                            molecule *larger = this;
                            if(allCells[cellContainer].allMultiMolecules[indexMultiMolecule]->molecules.size() < allCells[cellContainer].allMultiMolecules[m->indexMultiMolecule]->molecules.size())
                            {//my multimolecule is smaller than the other cells
                                smaller = this;
                                larger = m;
                            }
                            multiMolecule *smallMulti = allCells[cellContainer].allMultiMolecules[smaller->indexMultiMolecule];
                            multiMolecule *largeMulti = allCells[cellContainer].allMultiMolecules[larger->indexMultiMolecule];
                            largeMulti->centerOfMass.x*=mass;
                            largeMulti->centerOfMass.y*=mass;
                            largeMulti->centerOfMass.z*=mass;
                            for(unsigned int k=0; k<smallMulti->molecules.size(); k++)
                            {
                                largeMulti->molecules.push_back(smallMulti->molecules[k]);
                                largeMulti->mass+=allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->mass;
                                largeMulti->centerOfMass.x+=allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->position.x*allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->mass;
                                largeMulti->centerOfMass.y+=allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->position.y*allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->mass;
                                largeMulti->centerOfMass.z+=allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->position.z*allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->mass;
                            }
                            largeMulti->centerOfMass.x/=largeMulti->mass;
                            largeMulti->centerOfMass.y/=largeMulti->mass;
                            largeMulti->centerOfMass.z/=largeMulti->mass;
                            
                            int smallMultiIndex = smaller->indexMultiMolecule;
                            for(unsigned int k=0; k<smallMulti->molecules.size(); k++)
                            {
                                allCells[cellContainer].allMolecules[smallMulti->molecules[k]]->indexMultiMolecule=larger->indexMultiMolecule;
                            }
                            
                            //allCells[cellContainer].allMultiMolecules.erase(allCells[cellContainer].allMultiMolecules.begin()+smaller->indexMultiMolecule);
                            delete smallMulti;//allCells[cellContainer].allMultiMolecules[smaller->indexMultiMolecule];
                            allCells[cellContainer].allMultiMolecules[smallMultiIndex] = NULL;
                            allCells[cellContainer].deletedMultiMolecules.push(smallMultiIndex);
                            return larger->indexMultiMolecule;
                        }
                    }
                }
            }
        }
    }
    return -1;
}
void molecule::teleport(point p)
{
    //pointArray temp(0,0,0);
    vector<point> currPts;
    vector<point> newPts;
    /*for(int i=0; i<3; i++)
    {
        for(int j=-1; j<=1; j+=2)
        {
            temp[i]=j;
            currPts.push_back(point(position.x+radius*temp[0], position.y+radius*temp[1], position.z+radius*temp[2]));
            newPts.push_back(point(p.x+radius*temp[0], p.y+radius*temp[1], p.z+radius*temp[2]));
        }
        temp[i]=0;
    }
    
    for(unsigned int i=0; i<currPts.size(); i++) //remove duplicates
    {
        for(unsigned int j=i+1; j<currPts.size(); j++)
        {
            if((int)currPts[i].x==(int)currPts[j].x && (int)currPts[i].y==(int)currPts[j].y && (int)currPts[i].z==(int)currPts[j].z) //same grid cube
            {
                currPts.erase(currPts.begin()+j);
                j--;
            }
        }
    }
    for(unsigned int i=0; i<newPts.size(); i++) //remove duplicates
    {
        for(unsigned int j=i+1; j<newPts.size(); j++)
        {
            if((int)newPts[i].x==(int)newPts[j].x && (int)newPts[i].y==(int)newPts[j].y && (int)newPts[i].z==(int)newPts[j].z) //same grid cube
            {
                newPts.erase(newPts.begin()+j);
                j--;
            }
        }
    }*/
    getTouchingPoints(currPts, position);
    getTouchingPoints(newPts, p);
    
    for(unsigned int i=0; i<currPts.size(); i++) //remove self from parts of map that I left
    {
        bool contained=false;
        for(unsigned int j=0; j<newPts.size(); j++)
        {
            if((int)newPts[j].x==(int)currPts[i].x && (int)newPts[j].y==(int)currPts[i].y && (int)newPts[j].z==(int)currPts[i].z) //same pt
            {
                contained = true;
                break;
            }
        }
        if(!contained) // the molecule has left the gridpoint currPts[i]
        {
            unsigned int a=0;
            for(; a<allCells[cellContainer].volume[(int)currPts[i].z][(int)currPts[i].y][(int)currPts[i].x].size(); a++)
            {
                if(allCells[cellContainer].volume[(int)currPts[i].z][(int)currPts[i].y][(int)currPts[i].x][a]==index)
                {
                    break;
                }
            }
            allCells[cellContainer].volume[(int)currPts[i].z][(int)currPts[i].y][(int)currPts[i].x].erase(allCells[cellContainer].volume[(int)currPts[i].z][(int)currPts[i].y][(int)currPts[i].x].begin()+a);
        }
    }
    for(unsigned int i=0; i<newPts.size(); i++) //insert self into parts of map that I entered
    {
        bool contained=false;
        for(unsigned int j=0; j<currPts.size(); j++)
        {
            if((int)newPts[i].x==(int)currPts[j].x && (int)newPts[i].y==(int)currPts[j].y && (int)newPts[i].z==(int)currPts[j].z) //same pt
            {
                contained = true;
                break;
            }
        }
        if(!contained) // the molecule has entered the gridpoint newPts[i]
        {
            allCells[cellContainer].volume[(int)newPts[i].z][(int)newPts[i].y][(int)newPts[i].x].push_back(index);
        }
    }
    position.x = p.x;
    position.y = p.y;
    position.z = p.z;
    /*if((int)position.x!=(int)(p.x) || (int)position.y!=(int)(p.y) || (int)position.z!=(int)(p.z))
    {
        unsigned int i=0;
        for(; i<allCells[cellContainer].volume[(int)position.z][(int)position.y][(int)position.x].size(); i++)
        {
            if(allCells[cellContainer].volume[(int)position.z][(int)position.y][(int)position.x][i]==index)
            {
                break;
            }
        }
        allCells[cellContainer].volume[(int)position.z][(int)position.y][(int)position.x].erase(allCells[cellContainer].volume[(int)position.z][(int)position.y][(int)position.x].begin()+i);
        position.x = p.x;
        position.y = p.y;
        position.z = p.z;
        allCells[cellContainer].volume[(int)position.z][(int)position.y][(int)position.x].push_back(index);
    }
    else
    {
        position.x = p.x;
        position.y = p.y;
        position.z = p.z;
    }*/
}

void molecule::getTouchingPoints(vector<point> &p, point center)
{
    pointArray temp(0,0,0);
    for(int i=0; i<3; i++)
    {
        for(int j=-1; j<=1; j+=2)
        {
            temp[i]=j;
            getTouchingPointsHelper(p, pointArray(center.x+radius*temp[0], center.y+radius*temp[1], center.z+radius*temp[2]), 0);          
        }
        temp[i]=0;
    }
    
    for(unsigned int i=0; i<p.size(); i++) //remove duplicates
    {
        for(unsigned int j=i+1; j<p.size(); j++)
        {
            if((int)p[i].x==(int)p[j].x && (int)p[i].y==(int)p[j].y && (int)p[i].z==(int)p[j].z) //same grid cube
            {
                p.erase(p.begin()+j);
                j--;
            }
        }
    }
}

void molecule::getTouchingPointsHelper(vector<point> &p, pointArray center, int axis)
{
    if(abs(((pointArray)center)[axis]-(int)(((pointArray)center)[axis]))<EPSILON) //on an edge
    {
        if(axis<2) //not last axis
        {
            center[axis]-=EPSILON;
            getTouchingPointsHelper(p, center, axis+1);
            center[axis]+=2*EPSILON;
            getTouchingPointsHelper(p, center, axis+1);
        }
        else
        {
            center[axis]-=EPSILON;
            p.push_back((point)center);
            center[axis]+=2*EPSILON;
            p.push_back((point)center);
        }
    }
    else
    {
        if(axis<2)
            getTouchingPointsHelper(p,center,axis+1);
        else
            p.push_back((point)center);
    }
}