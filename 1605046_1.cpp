// #include<stdio.h>
// #include<stdlib.h>
// #include<math.h>
#include <windows.h>
#include <GL/glut.h>
#include <bits/stdc++.h>

using namespace std;

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

struct point{
	double x,y,z;
};

int stepUnit = 5;
float stepAngle = 2;

class Vector{
public:
	double x, y, z;
	Vector(){
	}
	Vector(double px, double py, double pz){
		x = px, y = py, z = pz;
	}
	void copyIt(Vector a){
		x = a.x, y = a.y, z = a.z;
	}

	Vector cross(Vector a){
		return Vector(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	}
	Vector add(Vector a){
		return Vector(x + a.x, y + a.y, z + a.z);
	}

	Vector multiply(double val){
		return Vector(val*x, val*y, val*z);
	}

	void shiftInDirection(Vector to, float steps){
		x += (to.x)*steps;
		y += (to.y)*steps;
		z += (to.z)*steps;
	}

	void rotate(Vector per, float angle){
		Vector t = cross(per);
		Vector m = *this;
		m = m.multiply(cos(angle*pi/180.0));
		t = t.multiply(sin(angle*pi/180.0));
		m = t.add(m);
		copyIt(m);
	}

	void print(){
		cout << x << " " << y << " " << z << "\n";
	}
};

Vector pos, u, r, l;

void drawAxes()
{
	int range = 500;
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( range,0,0);
			glVertex3f(-range,0,0);

			glVertex3f(0,-range,0);
			glVertex3f(0, range,0);

			glVertex3f(0,0, range);
			glVertex3f(0,0,-range);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawHalfCircle(float radius, int hSlices, int vSlices){
	vector< vector< point > > points(hSlices+1, vector<point>(vSlices+1) );
	for(int i = 0; i <= hSlices; i++){
		double angle = (pi/2) * i/(double)hSlices;
		double r = radius * cos(angle);
		double z = radius * sin(angle);

		for(int j = 0; j <= vSlices; j++){
			angle = (2*pi) * j/(double) vSlices;
			points[i][j].x = r * cos(angle);
			points[i][j].y = r * sin(angle);
			points[i][j].z = z;
		}
	}

	for(int i = 0; i < (int)points.size() - 1; i++){
		for(int j = 0; j < (int)points[i].size() - 1; j++){
			int white = (j%2);
		    glColor3f( white , white, white);
			glBegin(GL_QUADS);{
				glVertex3f(points[i][j].x , points[i][j].y , points[i][j].z);				
				glVertex3f(points[i][j+1].x , points[i][j+1].y , points[i][j+1].z);				
				glVertex3f(points[i+1][j+1].x , points[i+1][j+1].y , points[i+1][j+1].z);				
				glVertex3f(points[i+1][j].x , points[i+1][j].y , points[i+1][j].z);				
			}glEnd();			
		}
	}	
}

void drawCylinder(float radius, int height, int hSlices, int vSlices){
	vector< vector< point > > points(hSlices+2, vector<point>(vSlices+1) );

	for(int j = 0; j <= vSlices; j++){
		angle = (2*pi) * j/(double) vSlices;
		points[0][j].x = radius * cos(angle);
		points[0][j].y = radius * sin(angle);
		points[0][j].z = 0;
	}

	for(int i = 0; i <= hSlices; i++){
		double angle = (pi/2) * i/(double)hSlices;
		double r = 2*radius - radius * cos(angle);
		double z = height +  radius * sin(angle);

		for(int j = 0; j <= vSlices; j++){
			angle = (2*pi) * j/(double) vSlices;
			points[i+1][j].x = r * cos(angle);
			points[i+1][j].y = r * sin(angle);
			points[i+1][j].z = z;
		}
	}

	for(int i = 0; i < (int)points.size() - 1; i++){
		for(int j = 0; j < (int)points[i].size() - 1; j++){
			int white = j%2;
		    glColor3f( white , white, white);
			glBegin(GL_QUADS);{
				glVertex3f(points[i][j].x , points[i][j].y , points[i][j].z);				
				glVertex3f(points[i][j+1].x , points[i][j+1].y , points[i][j+1].z);				
				glVertex3f(points[i+1][j+1].x , points[i+1][j+1].y , points[i+1][j+1].z);				
				glVertex3f(points[i+1][j].x , points[i+1][j].y , points[i+1][j].z);				
			}glEnd();			
		}
	}

}


int figAngle = 0;
int halfFigAngle = 0;
int cylAngle = 0;
int cylRotAngle = 0;

int stripes = 50;
int circleRad = 30;
int cylRad = 15;
int cylHeight = 70;
int slices = 20;
int squareSide = 70;
int squareDist = 200;
int bulletSize = 2;
bool showAxes= true;

vector<point> shots;

void drawSS()
{
    glColor3f(1,0,0);

	glRotatef(90, 1, 0, 0);

    glPushMatrix();

	glRotatef(figAngle, 0, 1, 0);

	drawHalfCircle(circleRad , slices, stripes);
	glRotatef(180, 0, 1, 0);

	glRotatef(-halfFigAngle, 1, 0, 0);

	drawHalfCircle(circleRad, slices, stripes);

	glTranslatef(0, 0, circleRad + cylRad);
	glRotatef(180, 0, 1, 0);

	glRotatef(cylAngle, 1, 0, 0);
	glRotatef(cylRotAngle, 0, 0, 1);

	drawHalfCircle(cylRad, slices, stripes);
	glRotatef(180, 0, 1, 0);
	drawCylinder(cylRad, cylHeight, slices, stripes);

    glPopMatrix();

    glColor3f(.2,.2,.2);

	glBegin(GL_QUADS);{
		glVertex3f(-squareSide , -squareSide , -cylHeight - squareDist);
		glVertex3f(squareSide , -squareSide , -cylHeight - squareDist);
		glVertex3f(squareSide , squareSide , -cylHeight - squareDist);
		glVertex3f(-squareSide , squareSide, -cylHeight - squareDist);
	}glEnd();

    glColor3f(1, 0, 0);
	for(int i = 0; i < (int)shots.size(); i++){
		glBegin(GL_QUADS);{
			glVertex3f(shots[i].x + bulletSize, shots[i].y + bulletSize, shots[i].z);
			glVertex3f(shots[i].x + bulletSize, shots[i].y -bulletSize, shots[i].z);
			glVertex3f(shots[i].x -bulletSize, shots[i].y -bulletSize, shots[i].z);
			glVertex3f(shots[i].x -bulletSize, shots[i].y + bulletSize, shots[i].z);
		}glEnd();
	}
    glColor3f(.2,.2,.2);
}

bool inLimit(int x){
	return x <= 45 && x >= -45;
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			// drawgrid=1-drawgrid;
			l.rotate(u, -stepAngle);			
			r.rotate(u, -stepAngle);			
			break;
		case '2':
			// drawgrid=1-drawgrid;
			l.rotate(u, stepAngle);			
			r.rotate(u, stepAngle);			
			break;
		case '3':
			// drawgrid=1-drawgrid;
			l.rotate(r, -stepAngle);			
			u.rotate(r, -stepAngle);			
			break;
		case '4':
			// drawgrid=1-drawgrid;
			l.rotate(r, stepAngle);			
			u.rotate(r, stepAngle);			
			break;
		case '5':
			// drawgrid=1-drawgrid;
			u.rotate(l, stepAngle);
			r.rotate(l, stepAngle);
			break;
		case '6':
			// drawgrid=1-drawgrid;
			u.rotate(l, -stepAngle);
			r.rotate(l, -stepAngle);
			break;
		case 'q':
			if( figAngle < 45)
				figAngle += stepAngle;
			break;
		case 'w':
			if( figAngle > -45)
				figAngle -= stepAngle;		
			break;
		case 'e':
			if( halfFigAngle < 45)		
				halfFigAngle += stepAngle/2;
			break;
		case 'r':
			if( halfFigAngle > -45)		
				halfFigAngle -= stepAngle/2;
			break;
		case 'a':
			if( cylAngle < 45)
				cylAngle += stepAngle/2;
			break;
		case 's':
			if( cylAngle > -45)
				cylAngle -= stepAngle/2;
			break;
		case 'd':
			cylRotAngle += stepAngle/2;
			break;
		case 'f':
			cylRotAngle -= stepAngle/2;
			break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraHeight -= 3.0;
			pos.shiftInDirection(l, -stepUnit);
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;
			pos.shiftInDirection(l, stepUnit);
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			// pos.x -= stepUnit;
			pos.shiftInDirection(r, stepUnit);
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			pos.shiftInDirection(r, -stepUnit);
			break;

		case GLUT_KEY_PAGE_UP:
			pos.shiftInDirection(u, stepUnit);
			break;
		case GLUT_KEY_PAGE_DOWN:
			pos.shiftInDirection(u, -stepUnit);
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			pos.print();
			l.print();
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){
				Vector lineDir = Vector(tan( figAngle*pi/180.0 ), tan( (halfFigAngle + cylAngle)*pi/180.0 ), 1);
				double zDist = -(cylHeight + squareDist);
				double xDist = lineDir.x * zDist;
				double yDist = lineDir.y * zDist;
				if(xDist  <= squareSide && xDist >= -squareSide && yDist <= squareSide && yDist >= -squareSide){
					point p;
					p.x = xDist, p.y = -yDist, p.z = zDist+1;
					shots.push_back(p);
				}
			}
			break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	// gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	// gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(pos.x, pos.y, pos.z,
			pos.x + l.x, pos.y + l.y, pos.z + l.z,	
			u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	//clear the screen
	glClearColor(0,0,0,0);


	u = Vector(0, 0, 1);
	r = Vector(-1/sqrt(2), 1/sqrt(2), 0);
	l = Vector(-1/sqrt(2), -1/sqrt(2), 0);
	pos = Vector(100, 100, 0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
