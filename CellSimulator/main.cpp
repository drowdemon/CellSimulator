#include <cstdlib>
#include <ctime>
#include <fstream>
#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>

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

double square(double val)
{
    return val*val;
}

/*void renderScene()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
            // Set the camera
    gluLookAt(500.0f, 500.0f, 500.0f,
              0.0f, 0.0f, 0.0f,
              0.0f, 1.0f, 0.0f);

    glutSwapBuffers();
}

void timerProc(int arg)
{
    glutTimerFunc(50,timerProc,0);
    
    // Reset transformations
    glLoadIdentity();
    // Set the camera
    gluLookAt(500.0f, 500.0f, 500.0f,
              0.0f, 0.0f, 0.0f,
              0.0f, 1.0f, 0.0f);
    
    RGB black(0,0,0);
    
    for(int i=0; i<map.size(); i++)
    {
        for(int j=0; j<map[i].size(); j++)
        {
            for(int k=0; k<map[i][j].size(); k++)
            {
                if(map[i][j][k]>0)
                {
                    glPushMatrix();
                        glTranslated(allMolecules[map[i][j][k]-1]->position.x, allMolecules[map[i][j][k]-1]->position.y, allMolecules[map[i][j][k]-1]->position.z);
                        glutSolidSphere(allMolecules[map[i][j][k]-1]->radius, 16,16);
                    glPopMatrix();
                }
            }
        }
    }
    
    glutSwapBuffers();
}

int main(int argc, char **argv)
{
    map.resize(200);
    for(int i=0; i<200; i++)
    {
        map[i].resize(200);
    }
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(50,30);
    glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)-80,glutGet(GLUT_SCREEN_HEIGHT)-60);
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
}*/

int main()
{
    srand(time(NULL));
    return 0;
}