#ifndef PROTEIN_H
#define	PROTEIN_H

#include <vector>
#include <set>
#include "dataStructs.h"
#include "molecule.h"

using namespace std;

#define LISTVARPROTEIN \
    X(int, type) \
    X(bool, active) 

class protein : molecule
{
public:
    int type; //0=structural, 1=activator, 2=repressor, 3=transporter, 4=create bond, 5=destroy bond, 6=channel, 7=gated channel, 8=digestor, SPECIAL: 9=DNA polymerase, 10=RNA polymerase, 11=ribosome
    //int actOnID;
    //int actOnIDB;
    bool active;
    
    vector<int> aminoAcids; //stores full information, the rest is derived from this, but useful to have in a nicer format
    void (protein::*act[12])();
    set<unsigned long long> actOnID;
    set<unsigned long long> actOnIDB;
    //constructor
#define X(type, name) \
    type p_ ## name,
    
    protein(LISTVARPROTEIN LISTVARMOLECULE bool extra);
#undef X
    //end constructor
    unsigned long long getID();
    //9 act functions
    void structuralAct();   //0
    void activatorAct();    //1
    void repressorAct();    //2
    void transporterAct();  //3
    void bondCreatorAct();  //4
    void bondDestroyerAct();//5
    void channelAct();      //6
    void gatedChannelAct(); //7
    void digestorAct();     //8
    //3 special act functions
    void DNAPolymeraseAct();//9
    void RNAPolymeraseAct();//10
    void ribosomeAct();     //11
    //end act functions
private:
    molecule* pickAdjacentProtein(vector<molecule*> *p, set<unsigned long long> &matchID);
};

#endif	/* PROTEIN_H */