#ifndef WORLD_H
#define	WORLD_H

#include <vector>
#include "protein.h"
//#include "cytosol.h"

using namespace std;

class multiMolecule;

class cell
{
public:
    vector<molecule*> allMolecules;
    vector<multiMolecule*> allMultiMolecules;
    vector<vector<vector<vector<int> > > > volume; //3 spatial dimensions and then all of the molecules in that 1x1x1 grid cube
    //vector<vector<vector<cytosol> > > volume;
    int index;
    cell(int i);
    void timeStep(); //TODO finish
private:
    void solve2dLine(double intercept1, double intercept2, double start, double velocity1, double velocity2, double velocity3, vector<pointArray> &out);
    bool addGridCubes(unsigned int counter, int flagWhere, pointArray vel, pointArray curLoc, vector<point> &gridCubes);
    void getGridCubesHelper(point *vel, point *pos, vector<point> *gridCubes); //must delete return value when this is called!!!
    vector<point>* getGridCubes(point *vel, point *pos, molecule *m);
    void checkCollisionAndMove(unsigned int i, double timeRemaining, int prevCollidedWith=-1);
};

#endif	/* WORLD_H */

