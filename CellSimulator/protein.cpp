#include <cstdlib>
#include "molecule.h"
#include "protein.h"
#include "cell.h"
#include "cytosol.h"
#include "globals.h"

using namespace std;

#define X(type, name) \
    p_ ## name, 
protein::protein(int p_type, bool p_active, int p_index, point p_position, unsigned long long p_id, int p_cellContainer, int p_moleculeType, point p_velocity, double p_radius, double p_mass, bool extra) : molecule(LISTVARMOLECULE extra)
{
#undef X
    
#define X(type, name) \
    name = p_ ## name;
    
    LISTVARPROTEIN
#undef X
    
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
    vector<molecule*> *p = checkAround(radius*(1+PROTEINCANFUNCTIONRADIUSFRACTION));
    molecule* ret = pickAdjacentProtein(p, actOnID);
    if(p->size()>0)
        ((protein*)ret)->active=true; //activate!
    delete p;
}
void protein::repressorAct()
{
    if(!active) //uh uh
        return;
    vector<molecule*> *p = checkAround(radius*(1+PROTEINCANFUNCTIONRADIUSFRACTION));
    molecule* ret = pickAdjacentProtein(p, actOnID);
    if(p->size()>0)
        ((protein*)ret)->active=false; //deactivate!
    delete p;
}
void protein::transporterAct()
{
}
void protein::bondCreatorAct()
{
}
void protein::bondDestroyerAct()
{
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

molecule* protein::pickAdjacentProtein(vector<molecule*> *p, set<unsigned long long> &matchID)
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