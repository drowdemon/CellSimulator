#include <cstdlib>
#include <cmath>

#include "molecule.h"
#include "globals.h"

molecule::molecule(int p_index, point p_position, unsigned long long p_id, int p_cellContainer, int p_moleculeType, point p_velocity, double p_radius, double p_mass, bool extra)
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
                            if(square(searchRadius+allCells[cellContainer].allMolecules[(*tempCells)[h]-1]->radius) < square(allCells[cellContainer].allMolecules[(*tempCells)[h]-1]->position.x-position.x) + square(allCells[cellContainer].allMolecules[(*tempCells)[h]-1]->position.y-position.y) + square(allCells[cellContainer].allMolecules[(*tempCells)[h]-1]->position.z-position.z)) //intersect
                            {
                                ret->push_back(allCells[cellContainer].allMolecules[(*tempCells)[h]-1]);
                            }
                        }
                    }
                }
            }
        }
    }
    return ret;
}
bool molecule::bond(int indexWithWhat) //make sure this isn't called when it shouldn't be. Use with great care. It should only be called when molecules actually collide, or when an enzymatic protein bonds things together. It makes no checks on the legality of the bond based on position/velocity, only on type of bond
{
    for(unsigned int i=0; i<bondsNaturallyWith.size(); i++)
    {
        if(!usedBond[i] && bondsNaturallyWith[i].find(allCells[cellContainer].allMolecules[indexWithWhat]->id)!=bondsNaturallyWith[i].end()) //can bond with this type of molecule, and didn't already
        {
            for(unsigned int j=0; j<allCells[cellContainer].allMolecules[indexWithWhat]->bondsNaturallyWith.size(); j++)
            {
                if(!allCells[cellContainer].allMolecules[indexWithWhat]->usedBond[j] && allCells[cellContainer].allMolecules[indexWithWhat]->bondsNaturallyWith[j].find(id)!=allCells[cellContainer].allMolecules[indexWithWhat]->bondsNaturallyWith[j].end()) //the molecule I'm bonding with has an empty slot
                {
                    //start bonding
                    bondedWith.push_back(indexWithWhat);
                    allCells[cellContainer].allMolecules[indexWithWhat]->bondedWith.push_back(index);
                    usedBond[i]=true;
                    allCells[cellContainer].allMolecules[indexWithWhat]->usedBond[j]=true;
                    velocity.x=velocity.y=velocity.z=0;
                    allCells[cellContainer].allMolecules[indexWithWhat]->velocity=velocity;
                    return true;
                }
            }
        }
    }
    return false;
}
bool molecule::bond(molecule *m) //make sure this isn't called when it shouldn't be. Use with great care. It should only be called when molecules actually collide, or when an enzymatic protein bonds things together. It makes no checks on the legality of the bond based on position/velocity, only on type of bond
{
    for(unsigned int i=0; i<bondsNaturallyWith.size(); i++)
    {
        if(!usedBond[i] && bondsNaturallyWith[i].find(m->id)!=bondsNaturallyWith[i].end()) //can bond with this type of molecule, and didn't already
        {
            for(unsigned int j=0; j<m->bondsNaturallyWith.size(); j++)
            {
                if(!m->usedBond[j] && m->bondsNaturallyWith[j].find(id)!=m->bondsNaturallyWith[j].end()) //the molecule I'm bonding with has an empty slot
                {
                    //start bonding
                    bondedWith.push_back(m->index);
                    m->bondedWith.push_back(index);
                    usedBond[i]=true;
                    m->usedBond[j]=true;
                    velocity.x=velocity.y=velocity.z=0;
                    m->velocity=velocity;
                    return true;
                }
            }
        }
    }
    return false;
}