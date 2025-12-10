#include <string.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <GL/glui.h>
#include <GL/glut.h>

#include "geometry.h"


extern double point_rad;// = 5;
extern int sx; // = 800;
extern int sy; // = 600;

float angles[2]{0.0, 0.0};

std::vector<Point> points = {{300, 400}, {600, 400}};
bool pressed[2] = {false, false};

Point& A = points[0];
Point& B = points[1];
Digon D(A, B, M_PI / 4);
Point mouse;

int main_window;


/**************************************** myGlutKeyboard() **********/
void myGlutKeyboard(unsigned char Key, int x, int y) {
    switch(Key) {
        case 27:
        case 'q':
            exit(0);
        break;
    };
    glutPostRedisplay();
}


/***************************************** myGlutMenu() ***********/
void myGlutMenu(int value) {
    myGlutKeyboard(value, 0, 0);
}

/***************************************** myGlutMouse() **********/
void myGlutMouse(int button, int button_state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN) {
        //Point Q = {mouse.x, sy - mouse.y};    
        for (int i = 0; i < points.size(); ++i) {
            if (dist(mouse, points[i]) <= point_rad) {
                pressed[i] = true;
            }       
        }
    }
}


/***************************************** myGlutMotion() **********/
void myGlutMotion(int x, int y) {
    for (int i = 0; i < points.size(); ++i) {
        if (pressed[i]) {
            points[i].x += x - mouse.x;
            points[i].y += y - mouse.y;
        }
    }
    
    //std::cout << x << "  " << y << std::endl;

    mouse.x = x;
    mouse.y = y;

    glutPostRedisplay();
}

void myGlutPassiveMotion(int x, int y) {
    mouse.x = x;
    mouse.y = y;
    
    for (int i = 0; i < points.size(); ++i) {
        pressed[i] = false;
    }

    glutPostRedisplay();
}

/**************************************** myGlutReshape() *************/
void myGlutReshape( int x, int y ) {
    glViewport(0, 0, x, y);
    
    mouse.y += y - sy;
    sx = x;
    sy = y; 
    
    glutPostRedisplay();
}

/***************************************** myGlutDisplay() *****************/

void myGlutDisplay() {
    /*glClearColor( .9f, .9f, .9f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum( -xy_aspect*.08, xy_aspect*.08, -.08, .08, .1, 15.0 );
    

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glTranslatef( 0.0f, 0.0f, -1.6f );
    glRotatef( rotationY, 0.0, 1.0, 0.0 );
    glRotatef( rotationX, 1.0, 0.0, 0.0 );

    /*** Now we render object, using the variables 'obj', 'segments', and
    'wireframe'.  These are _live_ variables, which are transparently
    updated by GLUI ***/
    
    glClearColor( .9f, .9f, .9f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glFrustum( -xy_aspect*.08, xy_aspect*.08, -.08, .08, .1, 15.0 );
    
 
    // Given the coordinates
    //gluOrtho2D(0.0, 800.0, 0.0, 600.0);
    //gluOrtho2D(0.0, sx, 0.0, sy);
    gluOrtho2D(0, sx, 0, sy);    
    glViewport(0, 0, sx, sy);    
    
    
    D.draw();
    for (const Point& P : points) {
        //Point Q = {mouse.x, sy - mouse.y};    
		if (dist(P, mouse) <= point_rad) {    
		    glColor3d(1.0, 0.3, 0.3);
		} else {
		    glColor3d(0.3, 0.3, 1.0);
		}
		P.draw();   
    }
    glColor3d(0.3, 0.3, 1.0);
    
    glFlush();
    glutSwapBuffers();
}

void bisector_cb(int control) {
    D.bisector = !D.bisector;
    glutPostRedisplay();
}

void alpha_trans(int control) {
    static int prev = 0;
    static int angle = 0;
    static int digon = 45;
    
    angle = static_cast<int>(angles[0]);
    digon += angle - prev;
    prev = angle;
    if (digon <= 0) {
        digon = 1;
    }    
    if (digon >= 360) {
        digon = 359;
    }
    D.alpha = M_PI / 180 * digon;
    
    glutPostRedisplay();
}

void orientation_trans(int control) {
    static int prev = 0;
    static int angle = 0;
    static int orient = 0;
    
    angle = static_cast<int>(angles[1]);
    orient += angle - prev;
    prev = angle;    
    if (orient >= 180) {
        orient = 179;
    }    
    if (orient <= -180) {
        orient = -179;
    }
    D.orientation = M_PI / 180 * orient;    
    
    glutPostRedisplay();
}



/**************************************** main() ********************/

int main(int argc, char* argv[]) {
    /****************************************/
    /*   Initialize GLUT and create window  */
    /****************************************/

    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE /*| GLUT_DEPTH */);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(sx, sy);

    main_window = glutCreateWindow("Moebius Digon");
    glutDisplayFunc( myGlutDisplay );
    glutReshapeFunc( myGlutReshape );
    glutKeyboardFunc( myGlutKeyboard );
    glutMotionFunc( myGlutMotion );
    glutPassiveMotionFunc( myGlutPassiveMotion );    
    glutMouseFunc( myGlutMouse );

    /****************************************/
    /*       Set up OpenGL lights           */
    /****************************************/

    //GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
    //GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
    //GLfloat light0_position[] = {1.0f, 1.0f, 1.0f, 0.0f};

    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    //glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    //glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    //glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

    /****************************************/
    /*          Enable z-buferring          */
    /****************************************/

    //glEnable(GL_DEPTH_TEST);

    /****************************************/
    /*         Here's the GLUI code         */
    /****************************************/

    GLUI *glui = GLUI_Master.create_glui( "GLUI", 0, 800, 50 ); /* name, flags, x, and y */
    
    new GLUI_Checkbox(glui, "Bisector", 0, 0, bisector_cb);
    new GLUI_Translation(glui, "Digon angle", GLUI_TRANSLATION_Y, angles, 0, alpha_trans);    
    new GLUI_Translation(glui, "Digon position", GLUI_TRANSLATION_X, &angles[1], 0, orientation_trans); 
    
    new GLUI_Button(glui, "Quit", 0, (GLUI_Update_CB)exit);

    glui->set_main_gfx_window(main_window);

    /* We register the idle callback with GLUI, *not* with GLUT */
    //GLUI_Master.set_glutIdleFunc( myGlutIdle );
    GLUI_Master.set_glutIdleFunc( NULL );

    glutMainLoop();

    return EXIT_SUCCESS;
}
