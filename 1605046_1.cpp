#include "bitmap_image.hpp"
#include "1605046_Classes.h"

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
int stepUnit = 5;
float stepAngle = 2;

int imageWidth = 768;
int imageHeight = 768;
int level;
// int imageWidth = 200;
// int imageHeight = 200;

Vector3D pos, u, r, l;

void initialize(){

	objects.push_back( new Floor(1000, 20) );

	int n;
	ifstream inp;
	inp.open("scene.txt");
	inp >> level;

	inp >> imageWidth;
	imageHeight = imageWidth;
	inp >> n;
	for(int i = 0; i < n; i++){
		string figure;
		inp >> figure;
		if(figure == "sphere"){
			Vector3D point;
			double rad, shine;
			vector<double> col(3), light(4);
			inp >> point.x >> point.y >> point.z;
			inp >> rad;
			for(int j = 0; j < 3; j++) inp >> col[j];
			for(int j = 0; j < 4; j++) inp >> light[j];
			inp >> shine;			
			objects.push_back( new Sphere(point, rad, col, light, shine) );
		}else if(figure == "triangle"){
			vector<Vector3D> points(3);
			double shine;
			vector<double> col(3), light(4);
			for(int j = 0; j < 3; j++){
				inp >> points[j].x >> points[j].y >> points[j].z ;
			}
			for(int j = 0; j < 3; j++) inp >> col[j];
			for(int j = 0; j < 4; j++) inp >> light[j];
			inp >> shine;
			objects.push_back( new Triangle(points, col, light, shine) );
		}else if(figure == "general"){
			vector<double> val(10);
			Vector3D point;
			double shine, l, w, h;
			vector<double> col(3), light(4);
			for(int j = 0; j < 10; j++) inp >> val[j];

			inp >> point.x >> point.y >> point.z ;
			inp >> l >> w >> h;
			for(int j = 0; j < 3; j++) inp >> col[j];
			for(int j = 0; j < 4; j++) inp >> light[j];
			inp >> shine;
			objects.push_back( new Quadric(val[0], val[1], val[2], val[3], val[4], 
					val[5], val[6], val[7], val[8], val[9], point, l, w, h, col, light, shine) );
		}

	}

	inp >> n;
	for(int i = 0; i < n; i++){
		Vector3D point;
		vector<double> col(3);
		inp >> point.x >> point.y >> point.z;
		for(int j = 0; j < 3; j++) inp >> col[j];
		lights.push_back(Light(point, col));
	}

	inp.close();


	// objects.push_back( new Sphere(Vector3D(40.0, 0.0, 10.0), 10.0, {0.0, 1.0, 0.0}, {0.4, 0.2, 0.2, 0.2}, 10 ) ); 
	// objects.push_back( new Sphere(Vector3D(-30.0, 60.0, 20.0), 20.0, {0.0, 0.0, 1.0}, {0.2, 0.2, 0.4, 0.2}, 15 ) );
	// objects.push_back( new Sphere(Vector3D(-15.0, 15.0, 45.0), 15.0, {1.0, 1.0, 0.0}, {0.4, 0.3, 0.1, 0.2}, 5 ) );

	// objects.push_back( new Triangle({Vector3D(50, 30, 0), 
	// 								Vector3D(70, 60, 0), Vector3D(50, 45, 50)},
	// 								{1.0, 0.0, 0.0}, {0.4, 0.2, 0.1, 0.3}, 5) );

	// objects.push_back( new Triangle({Vector3D(70, 60, 0), 
	// 								Vector3D(30, 60, 0), Vector3D(50, 45, 50)},
	// 								{0.0, 1.0, 0.0}, {0.4, 0.2, 0.1, 0.3}, 5) );

	// objects.push_back( new Triangle({Vector3D(30, 60, 0), 
	// 								Vector3D(50, 30, 0), Vector3D(50, 45, 50)},
	// 								{0.0, 0.0, 1.0}, {0.4, 0.2, 0.1, 0.3}, 5) );

	// objects.push_back( new Quadric(0.0625, 0.04, 0.04, 0, 0, 0, 0, 0, 0, -36, Vector3D(0, 0, 0), 
	// 	0 ,0 ,15, {1.0, 0.0, 0.0}, {0.4, 0.2, 0.1, 0.3}, 15 ) );

	// lights.push_back( Light(Vector3D(70.0, 70.0, 70.0), {1.0, 0.0, 0.0}) );
	// lights.push_back( Light(Vector3D(-70, 70, 70), {0.0, 0.0, 1.0}) );
	// lights.push_back( Light(Vector3D(70, -70, 70), {1, 0, 0.0}));
	// lights.push_back( Light(Vector3D(-70, -70, 70), {0, 1.0, 0}));


}

void capture(){
	bitmap_image image(imageWidth, imageHeight);
	vector<vector< vector<int> >> img(imageWidth, vector<vector<int>>(imageHeight, vector<int>(3, 0)) );

	double planeDistance = 5;
	double windowWidth = 10;
	double windowHeight = 10;

	pos.print();
	l.print();
	u.print();
	r.print();

	Vector3D eye = pos;
	Vector3D right = r.multiply(windowWidth/2).multiply(-1);
	Vector3D up = u.multiply(windowHeight/2);
	Vector3D look = l.multiply(planeDistance);

	Vector3D topLeft = eye.add(look).add(right).add(up);

	double du = windowWidth/imageWidth;
	double dv = windowHeight/imageHeight;

	topLeft = topLeft.add(r.multiply(0.5*du)).add( u.multiply( - 0.5*dv) );

	int counter = 0;

	cout << "Started\n";
	auto start = chrono::high_resolution_clock::now();

	double col[] = {1, 1, 1};
	for(int x = 0; x < imageWidth; x++){
		for(int y = 0; y < imageHeight; y++){
			Vector3D cur = topLeft;
			cur = cur.add(r.multiply(x*du)).add( u.multiply( - y*dv) );
			Ray ray = Ray( eye, cur.add(eye.multiply(-1)).normalize() );
			double t_min = 1111111;
			for(int i = 0; i < (int)objects.size(); i++){
				double t = objects[i]->intersect(ray, col, level);
				if(t > 0 && t < t_min){
					t_min = t;
					for(int i = 0; i < 3; i++){
						col[i] = min(col[i], 1.0);
					}
					img[x][y] = {(int)(col[0]*255), (int)(col[1]*255), (int)(col[2]*255)};
				}
			}
			if(t_min < 10000) counter++;
		}
	}
	cout << "finished\n";
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
	cout << duration.count() << endl;	

	for(int i = 0; i < imageWidth; i++){
		for(int j = 0; j < imageHeight; j++){
			image.set_pixel(i,j, img[i][j][0] ,img[i][j][1], img[i][j][2]);
		}
	}

	image.save_image("output.bmp");
	cout << counter << "\n";
}

void drawSS()
{
	for(int i = 0; i < (int)objects.size(); i++){
		objects[i]->draw();
	}
	for(int i = 0; i < (int)lights.size(); i++){
		lights[i].draw();
	}

	// drawHalfCircle(30 , 20, 50);
}


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

bool showAxes= true;


bool inLimit(int x){
	return x <= 45 && x >= -45;
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			l.rotate(u, -stepAngle);			
			r.rotate(u, -stepAngle);			
			break;
		case '2':
			l.rotate(u, stepAngle);			
			r.rotate(u, stepAngle);			
			break;
		case '3':
			l.rotate(r, -stepAngle);			
			u.rotate(r, -stepAngle);			
			break;
		case '4':
			l.rotate(r, stepAngle);			
			u.rotate(r, stepAngle);			
			break;
		case '5':
			u.rotate(l, stepAngle);
			r.rotate(l, stepAngle);
			break;
		case '6':
			u.rotate(l, -stepAngle);
			r.rotate(l, -stepAngle);
			break;
		case '0':
			capture();
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
			break;

		case GLUT_RIGHT_BUTTON:
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


	u = Vector3D(0, 0, 1);
	r = Vector3D(0.422618, 0.906308, 0);
	l = Vector3D(-0.906308, 0.422618, 0);
	pos = Vector3D(143.633, -51.5565, 45);
	// u = Vector3D(0, 0, 1);
	// r = Vector3D(-1/sqrt(2), 1/sqrt(2), 0);
	// l = Vector3D(-1/sqrt(2), -1/sqrt(2), 0);
	// pos = Vector3D(100, 100, 0);

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
	initialize();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
