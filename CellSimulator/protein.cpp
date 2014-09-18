#include <cstdlib>
#include "molecule.h"
#include "protein.h"
#include "cell.h"
#include "globals.h"
#include "multiMolecule.h"

using namespace std;

#define X(type, name) \
    p_ ## name, 
protein::protein(int p_type, bool p_defaultActive, point p_specialFunction, int p_index, point p_position, unsigned long long p_id, int p_cellContainer, int p_moleculeType, double p_radius, double p_mass, int p_indexMultiMolecule, bool extra) : molecule(LISTVARMOLECULE extra)
{
#undef X
    
#define X(type, name) \
    name = p_ ## name;
    
    LISTVARPROTEIN
#undef X
    
    active = defaultActive;
    act[0] = &protein::structuralAct;
    act[1] = &protein::activatorAct;
    act[2] = &protein::repressorAct;
    act[3] = &protein::transporterAct;
    act[4] = &protein::bondCreatorAct;
    act[5] = &protein::bondDestroyerAct;
    act[6] = &protein::channelAct;
    act[7] = &protein::gatedChannelAct;
    act[8] = &protein::digestorAct;
    
    act[9] = &protein::DNAPolymeraseAct;
    act[10] = &protein::RNAPolymeraseAct;
    act[11] = &protein::ribosomeAct;
}

unsigned long long protein::getID()
{
    return 0;
}

void protein::structuralAct()
{
}
void protein::activatorAct()
{
    if(!active) //uh uh
        return;
    //vector<molecule*> *p = checkAround(radius*(1+PROTEINCANFUNCTIONRADIUSFRACTION));
    /*molecule* ret = pickAdjacentProtein(p, actOnID);
    if(p->size()>0)
        ((protein*)ret)->active=true; //activate!
    delete p;*/
    molecule* ret = pickProtein(bondedWith, actOnID);
    if(ret)
        ((protein*)ret)->active=true; //activate!
}
void protein::repressorAct()
{
    if(!active) //uh uh
        return;
    /*vector<molecule*> *p = checkAround(radius*(1+PROTEINCANFUNCTIONRADIUSFRACTION));
    molecule* ret = pickProtein(p, actOnID);
    if(p->size()>0)
        ((protein*)ret)->active=false; //deactivate!
    delete p;*/
    molecule* ret = pickProtein(bondedWith, actOnID);
    if(ret)
        ((protein*)ret)->active=false; //deactivate!
}
void protein::transporterAct()
{
}
void protein::bondCreatorAct()
{
    if(!active)
        return;
    molecule *m1 = pickProtein(bondedWith, actOnID);
    if(m1==NULL)
        return;
    vector<bondData> copy = bondedWith;
    for(unsigned int i=0; i<copy.size(); i++) //make sure same molecule isn't chosen twice
    {
        if(copy[i].with==m1->index)
        {
            copy.erase(copy.begin()+i);
            break;
        }
    }
    molecule *m2 = pickProtein(copy, actOnIDB);
    if(m2==NULL)
        return;
}
void protein::bondDestroyerAct()
{
    if(!active)
        return;
    molecule *m1 = pickProtein(bondedWith, actOnID);
    if(m1==NULL)
        return;
    molecule *m2 = pickProtein(m1->bondedWith, actOnIDB);
    if(m2==NULL)
        return;
    bool didOne=false;
    point vel = allCells[cellContainer].allMultiMolecules[m1->indexMultiMolecule]->velocity;
    for(unsigned int i=0; i<m1->bondedWith.size(); i++)
    {
        if(m1->bondedWith[i].with==index || m1->bondedWith[i].with==m2->index)
        {
            m1->bondedWith.erase(m1->bondedWith.begin()+i);
            if(didOne)
                break;
            else
                didOne=true;
            i--;
        }
    }
    didOne=false;
    for(unsigned int i=0; i<m2->bondedWith.size(); i++)
    {
        if(m2->bondedWith[i].with==index || m2->bondedWith[i].with==m1->index)
        {
            m2->bondedWith.erase(m2->bondedWith.begin()+i);
            if(didOne)
                break;
            else
                didOne=true;
            i--;
        }
    }
    didOne=false;
    for(unsigned int i=0; i<bondedWith.size(); i++)
    {
        if(bondedWith[i].with==m2->index || bondedWith[i].with==m1->index)
        {
            bondedWith.erase(bondedWith.begin()+i);
            if(didOne)
                break;
            else
                didOne=true;
            i--;
        }
    }
    
    set<int> contents1;
    multiMolecule *mm1 = outerBuildMultimolecule(contents1, vel, m1);
    int thisOneIsSeparate = 0;
    if(contents1.find(m2->index)!=contents1.end()) //if the second molecule is in the first one's multimolecule chain anyway, i.e. they're still in a multimolecule, just broke up a cycle a bit.
    {
        if(contents1.find(index)!=contents1.end()) // if this molecule is still in the multimolecule...
        {
            return; //I'm done here, it's /still/ all one big happy multimolecule
        }
        else
            thisOneIsSeparate++;
    }
    else //the second molecule has a second, separate multimolecule
    {
        delete allCells[cellContainer].allMultiMolecules[indexMultiMolecule];
        allCells[cellContainer].allMultiMolecules[indexMultiMolecule] = NULL;
        allCells[cellContainer].deletedMultiMolecules.push(indexMultiMolecule);
        addMultimolecule(mm1);
        set<int> contents2;
        multiMolecule *mm2 = outerBuildMultimolecule(contents2, vel, m2);
        if(contents2.find(index)!=contents2.end()) //this one is in the multimolecule
            thisOneIsSeparate++; //this molecule is part of something, no need to create it on its pwn
        addMultimolecule(mm2);
        
        //randomize velocity
        double rnd = (double)(rand()%100)/1000.0 + 0.9;
        if(rand()%2)
        {
            multiMolecule *temp = mm2;
            mm2=mm1;
            mm1=temp;
        }
        mm1->velocity.x*=rnd;
        mm1->velocity.y*=rnd;
        mm1->velocity.z*=rnd;
        mm2->velocity.x*=(2-rnd-(EPSILON*50));
        mm2->velocity.y*=(2-rnd-(EPSILON*50));
        mm2->velocity.z*=(2-rnd-(EPSILON*50));
    }
    if(thisOneIsSeparate == 0)
    {
        set<int> contentsThis;
        multiMolecule *mmThis = outerBuildMultimolecule(contentsThis, vel, this);
        addMultimolecule(mmThis);
        //randomize velocity
        double rnd = (double)(rand()%100)/1000.0 + 0.9;
        mmThis->velocity.x*=rnd;
        mmThis->velocity.y*=rnd;
        mmThis->velocity.z*=rnd;
    }
}
void protein::channelAct()
{
}
void protein::gatedChannelAct()
{
}
void protein::digestorAct()
{
}
void protein::DNAPolymeraseAct()
{
}
void protein::RNAPolymeraseAct()
{
}
void protein::ribosomeAct()
{
}

molecule* protein::pickProtein(vector<molecule*> *p, set<unsigned long long> &matchID)
{
    for(unsigned int i=0; i<p->size(); i++)
    {
        if((*p)[i]->moleculeType!=TYPE_PROTEIN || matchID.find((*p)[i]->id)==matchID.end())
        {
            p->erase(p->begin()+i);
            i--;
        }
    }
    return (*p)[rand()%p->size()];
}

molecule* protein::pickProtein(vector<bondData> p, set<unsigned long long> &matchID)
{
    for(unsigned int i=0; i<p.size(); i++)
    {
        if(allCells[0].allMolecules[p[i].with]->moleculeType!=TYPE_PROTEIN || matchID.find(allCells[0].allMolecules[p[i].with]->id)==matchID.end())
        {
            p.erase(p.begin()+i);
            i--;
        }
    }
    if(p.size()>0)
        return allCells[0].allMolecules[p[rand()%p.size()].with];
    else
        return NULL;
}

int protein::buildMultimolecule(set<int> &checkIndexes, molecule* act, int callerIndex, multiMolecule* build)
{
    build->molecules.push_back(index);
    build->centerOfMass.x+=position.x;
    build->centerOfMass.y+=position.y;
    build->centerOfMass.z+=position.z;
    build->mass+=mass;
    int sum = 0;
    for(unsigned int i=0; i<act->bondedWith.size(); i++)
    {
        if(callerIndex==act->bondedWith[i].with) //this is the caller, did it already
            continue;
        if(checkIndexes.insert(act->bondedWith[i].with).second) //this molecule was new
        {
            sum += buildMultimolecule(checkIndexes, allCells[cellContainer].allMolecules[act->bondedWith[i].with], index, build);
        }
    }
    return sum+1;
}
multiMolecule* protein::outerBuildMultimolecule(set<int> &contents, point &vel, molecule *m)
{
    multiMolecule *mm = new multiMolecule(0, point(), point(), cellContainer);
    int size = buildMultimolecule(contents, m, -1, mm);
    mm->mass/=size;
    mm->centerOfMass.x/=size;
    mm->centerOfMass.y/=size;
    mm->centerOfMass.z/=size;
    mm->velocity = vel;
    return mm;
}
void protein::addMultimolecule(multiMolecule* mm)
{
    if(!allCells[cellContainer].deletedMultiMolecules.empty())
    {
        int placementIndex = allCells[cellContainer].deletedMultiMolecules.top();
        allCells[cellContainer].deletedMultiMolecules.pop();
        for(unsigned int i=0; i<mm->molecules.size(); i++)
        {
            allCells[cellContainer].allMolecules[mm->molecules[i]]->indexMultiMolecule=placementIndex;
        }
        allCells[cellContainer].allMultiMolecules[placementIndex]=mm;
    }
    else
    {
        int placementIndex = allCells[cellContainer].allMultiMolecules.size();
        allCells[cellContainer].allMultiMolecules.push_back(mm);
        for(unsigned int i=0; i<mm->molecules.size(); i++)
        {
            allCells[cellContainer].allMolecules[mm->molecules[i]]->indexMultiMolecule=placementIndex;
        }
    }
}