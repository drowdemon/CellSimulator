#include <cstdlib>
#include <ctime>
#include <fstream>
#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include "cell.h"
#include "globals.h"
#include "multiMolecule.h"

using namespace std;

//TODO  finish cell::timeStep
//TODO finish molecule::collide
    //NOTE: when heavy molecule A collides with light molecule B, molecule B will have a velocity like molecule A, A's velocity will be mostly unchanged, but B will not have moved, so A will collide with it again... and again... and again... ugh
    //Perhaps a pseudo-bond: inelastic collision style solution, if they collide twice in the same timestep they are pseudo-bonded and move together as a unit for that timestep
//TODO molecule mvmt does not account for molecules that are bonded together. Shit.
    //TODO create class multiMolecule. It'll have a list of bonded molecules. These are allowed to intersect DURING A TIMESTEP, i.e. during movement, but at no other time. No way are they moving in formation neatly without intersecting...
        //multiMolecules get to move before molecules, and all the molecules in a multiMolecule move one after the other without any other molecules moving in between. If one hits an obstacle, the rest may have to go back and recalculate their movement.
//TODO make bondsNaturallyWith contain information on where on the molecule the bond can occur, which would consist of a vector of points (from the center of the molecule) and a vector of radiuses around that point within which a bond can be made. 
    //If molecule A can bond with 2 molecules of B, just have B in the entry twice, like you would now, with 1 pt each, not with both pts on each. But if A bonds with 1 C molecule at either of 2 locations, but both in the vector. Anywhere on the cell: pt=(0,0,0), radius = molecule_radius + epsilon

void createProtein(double mass, point velocity, point position, double radius, unsigned long long id);

int WIDTH;
int HEIGHT;
int mainWindow;

class RGB
{
public:
    unsigned char r;
    unsigned char g;
    unsigned char b;
    RGB(unsigned char pr, unsigned char pg, unsigned char pb)
    {
        r=pr;
        g=pg;
        b=pb;
    }
};

double square(double val)
{
    return val*val;
}

void resetTransformations()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-10, 10, -10, 10, -1000, 1000);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    /*glTranslated(20,20,20);
    gluLookAt(30.0f, 30.0f, 0.0f,
        20.0f, 20.0f, 20.0f,
        0.0f, 1.0f, 0.0f);*/
    gluLookAt(25.774f, 25.774f, 25.774f, 
        20.0f, 20.0f, 20.0f,
        0.0f, 1.0f, 0.0f);
        
}

void renderScene()
{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Reset transformations
    resetTransformations();

    // Just to see some triangles
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glColor3f(0.0, 0.0, 0.0);  //color = black
    /*glPushMatrix();
        glTranslated(0, 0, 0);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
    glPushMatrix();
        glTranslated(2, 0, 0);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
     glPushMatrix();
        glTranslated(2, 2, 0);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
    glPushMatrix();
        glTranslated(0, 2, 0);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
    glPushMatrix();
        glTranslated(0, 0, 2);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
    glPushMatrix();
        glTranslated(posX, 0, 2);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
     glPushMatrix();
        glTranslated(2, 2, 2);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();
    glPushMatrix();
        glTranslated(0, 2, 2);
        glutSolidSphere(.5, 50, 50);
    glPopMatrix();*/
    for(unsigned int i=0; i<allCells[0].volume.size(); i++)
    {
        for(unsigned int j=0; j<allCells[0].volume[i].size(); j++)
        {
            for(unsigned int k=0; k<allCells[0].volume[i][j].size(); k++)
            {
                for(unsigned int h=0; h<allCells[0].volume[i][j][k].size(); h++)
                {
                    glPushMatrix();
                        glTranslated(allCells[0].allMolecules[allCells[0].volume[i][j][k][h]]->position.x, allCells[0].allMolecules[allCells[0].volume[i][j][k][h]]->position.y, allCells[0].allMolecules[allCells[0].volume[i][j][k][h]]->position.z);
                        //glTranslated(0,0,0);
                        glutSolidSphere(allCells[0].allMolecules[allCells[0].volume[i][j][k][h]]->radius,50,50);
                        //glutSolidSphere(2, 50, 50);
                    glPopMatrix();
                }
            }
        }
    }

    glFlush();
    glutSwapBuffers();
}

void timerProc(int arg)
{
    glutTimerFunc(50,timerProc,0);
    
    for(unsigned int i=0; i<allCells.size(); i++)
    {
        allCells[0].timeStep();
    }
    
    glutPostRedisplay();
}

void init()
{
    allCells.push_back(cell(0));
    allCells[0].volume.resize(100);
    for(unsigned int i=0; i<allCells[0].volume.size(); i++)
    {
        allCells[0].volume[i].resize(100);
        for(unsigned int j=0; j<allCells[0].volume[i].size(); j++)
        {
            allCells[0].volume[i][j].resize(100);
        }
    }
    
    createProtein(1,point(-.1,0,0),point(22,20,20),.3,0);
    allCells[0].allMolecules[0]->possibleBonds.push_back(bondSpot(set<unsigned long long>(), false, point(-allCells[0].allMolecules[0]->radius*.8,0,0), allCells[0].allMolecules[0]->radius*.4));
    allCells[0].allMolecules[0]->possibleBonds[0].bondsNaturallyWith.insert(0);
    
    createProtein(1,point(.1,0,0),point(20,20,20),.3,0); //top left
    allCells[0].allMolecules[1]->possibleBonds.push_back(bondSpot(set<unsigned long long>(), false, point(allCells[0].allMolecules[0]->radius*.8,0,0), allCells[0].allMolecules[1]->radius*.4));
    allCells[0].allMolecules[1]->possibleBonds[0].bondsNaturallyWith.insert(0);
    allCells[0].allMolecules[1]->possibleBonds.push_back(bondSpot(set<unsigned long long>(), false, point(0,0,0), allCells[0].allMolecules[1]->radius*(1+EPSILON)));
    allCells[0].allMolecules[1]->possibleBonds[1].bondsNaturallyWith.insert(1);
    
    createProtein(1,point(0,0,-.15),point(20.6,20,25),.3,1);
    ((protein*)(allCells[0].allMolecules[2]))->active=true;
    ((protein*)(allCells[0].allMolecules[2]))->type=5;
    ((protein*)(allCells[0].allMolecules[2]))->actOnID.insert(0);
    ((protein*)(allCells[0].allMolecules[2]))->actOnIDB.insert(0);
    allCells[0].allMolecules[2]->possibleBonds.push_back(bondSpot(set<unsigned long long>(), false, point(0,0,0), allCells[0].allMolecules[2]->radius*(1+EPSILON)));
    allCells[0].allMolecules[2]->possibleBonds[0].bondsNaturallyWith.insert(0);
    
    for(unsigned int iii=0; iii<allCells[0].allMolecules.size(); iii++)
    {
        vector<point> currPts;
        allCells[0].allMolecules[iii]->getTouchingPoints(currPts, allCells[0].allMolecules[iii]->position);
        
        for(unsigned int k=0; k<currPts.size(); k++)
        {
            allCells[0].volume[(int)currPts[k].z][(int)currPts[k].y][(int)currPts[k].x].push_back(iii);
        }
    }
}

void createProtein(double mass, point velocity, point position, double radius, unsigned long long id)
{
    allCells[0].allMolecules.push_back(new protein(0, true, point(), allCells[0].allMolecules.size(), position, id, 0, TYPE_PROTEIN, radius, mass, allCells[0].allMultiMolecules.size()));
    allCells[0].allMultiMolecules.push_back(new multiMolecule(mass, velocity, position, 0));
    allCells[0].allMultiMolecules[allCells[0].allMultiMolecules.size()-1]->molecules.push_back(allCells[0].allMolecules.size()-1);
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    init();
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(250,30);
    glutInitWindowSize(glutGet(GLUT_SCREEN_HEIGHT)-60,glutGet(GLUT_SCREEN_HEIGHT)-60);
    //glutInitWindowSize(700,700);

    mainWindow=glutCreateWindow("Game Engine");

    //glutFullScreen();
    WIDTH=glutGet(GLUT_WINDOW_WIDTH);
    HEIGHT=glutGet(GLUT_WINDOW_HEIGHT);

    glutDisplayFunc(renderScene);
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    glutTimerFunc(50,timerProc,0);
    glutMainLoop();
    return 0;
}