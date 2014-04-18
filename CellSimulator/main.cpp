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
    
    allCells[0].allMolecules.push_back(new protein(0, true, 0, point(22, 20, 24), 0, 0, TYPE_PROTEIN, .3, .5, 0));
    //allCells[0].allMolecules[0]->bondsNaturallyWith.push_back(set<unsigned long long>());
    //allCells[0].allMolecules[0]->bondsNaturallyWith[0].insert(0);
    //allCells[0].allMolecules[0]->usedBond.push_back(false);
    allCells[0].allMultiMolecules.push_back(new multiMolecule(.5, point(-.1,0,-.2), point(22,20,24), 0));
    allCells[0].allMultiMolecules[0]->molecules.push_back(0);
    
    allCells[0].allMolecules.push_back(new protein(0, true, 1, point(20, 20, 20), 0, 0, TYPE_PROTEIN, .3, 1, 1));
    allCells[0].allMolecules[1]->bondsNaturallyWith.push_back(set<unsigned long long>());
    allCells[0].allMolecules[1]->bondsNaturallyWith[0].insert(0);
    allCells[0].allMolecules[1]->usedBond.push_back(false);
    allCells[0].allMultiMolecules.push_back(new multiMolecule(1, point(0,0,0), point(20,20,20), 0));
    allCells[0].allMultiMolecules[1]->molecules.push_back(1);
    
    for(unsigned int iii=0; iii<allCells[0].allMolecules.size(); iii++)
    {
        //pointArray temp(0,0,0);
        vector<point> currPts;
        /*for(int i=0; i<3; i++)
        {
            for(int j=-1; j<=1; j+=2)
            {
                temp[i]=j;
                currPts.push_back(point(allCells[0].allMolecules[iii]->position.x+allCells[0].allMolecules[iii]->radius*temp[0], allCells[0].allMolecules[iii]->position.y+allCells[0].allMolecules[iii]->radius*temp[1], allCells[0].allMolecules[iii]->position.z+allCells[0].allMolecules[iii]->radius*temp[2]));
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
        }*/
        allCells[0].allMolecules[iii]->getTouchingPoints(currPts, allCells[0].allMolecules[iii]->position);
        
        for(unsigned int k=0; k<currPts.size(); k++)
        {
            allCells[0].volume[(int)currPts[k].z][(int)currPts[k].y][(int)currPts[k].x].push_back(iii);
        }
    }
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

/*int main()
{
    srand(time(NULL));
    return 0;
}*/