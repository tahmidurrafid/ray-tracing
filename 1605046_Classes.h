#include <windows.h>
#include <GL/glut.h>
#include <bits/stdc++.h>

#define pi (2*acos(0.0))
using namespace std;


struct point{
	double x,y,z;
};

class Vector3D{
public:
	double x, y, z;
	Vector3D(){
	}
	Vector3D(double px, double py, double pz){
		x = px, y = py, z = pz;
	}
	void copyIt(Vector3D a){
		x = a.x, y = a.y, z = a.z;
	}

	Vector3D cross(Vector3D a){
		return Vector3D(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	}
    double dotProduct(Vector3D a){
        return x*a.x + y*a.y + z*a.z; 
    }

	Vector3D add(Vector3D a){
		return Vector3D(x + a.x, y + a.y, z + a.z);
	}

	Vector3D multiply(double val){
		return Vector3D(val*x, val*y, val*z);
	}

	void shiftInDirection(Vector3D to, float steps){
		x += (to.x)*steps;
		y += (to.y)*steps;
		z += (to.z)*steps;
	}

	void rotate(Vector3D per, float angle){
		Vector3D t = cross(per);
		Vector3D m = *this;
		m = m.multiply(cos(angle*pi/180.0));
		t = t.multiply(sin(angle*pi/180.0));
		m = t.add(m);
		copyIt(m);
	}

    Vector3D normalize(){
        double mul = sqrt(x*x + y*y + z*z );
        x = x/mul;
        y = y/mul;
        z = z/mul;
        return *this;
    }

	void print(){
		cout << x << " " << y << " " << z << "\n";
	}
};

class Ray{
public:
    Vector3D start;
    Vector3D dir;
    Ray(Vector3D st, Vector3D d){
        start = st;
        dir = d;
        dir.normalize();
    }

};


class Object{
public:
    Vector3D reference_point;
    double height, width, length;
    vector<double> color;
    vector<double> coEfficients;
    int shine;

    Object(Vector3D rp, double h, double w, double l, vector<double>  col, vector<double>  coEff, double sh){
        reference_point = rp;
        height = h, width = w, length = l;
        color = col;
        coEfficients = coEff;
        shine = sh;
    }

    virtual void draw(){
        cout << "Draw method not overriddent\n";
    }

    virtual double intersect(Ray ray, vector<double> &cols, int level){

    }

    void setColor(){
    }

    void setShine(){
    }

    void setCoefficient(){
    }
};

class Light{
public:
    Vector3D light_pos;
    vector<double> color;

    Light(Vector3D pos, vector<double> col){
        light_pos = pos;
        color = col;
    }

    void draw(){
        int a = 2;
        glColor3f(color[0], color[1], color[2]);        
        glBegin(GL_QUADS);{
            glVertex3f( light_pos.x + a, light_pos.y + a, light_pos.z+2);
            glVertex3f( light_pos.x + a,light_pos.y-a,light_pos.z+2);
            glVertex3f( light_pos.x -a,light_pos.y-a,light_pos.z+2);
            glVertex3f( light_pos.x -a,light_pos.y+a,light_pos.z+2);
    	}glEnd();
    }
};

class Sphere : public Object{
public:
    Sphere(Vector3D rp, double radius, vector<double>  col, vector<double>  coEff, double sh) 
        : Object(rp, radius, radius, radius, col, coEff, sh)
    {}

    double getRadius(){
        return height;
    }

    void draw(){
        int hSlices = 20;
        int vSlices = 20;
        vector< vector< point > > points(2*hSlices+1, vector<point>(vSlices+1) );
        for(int i = -hSlices; i <= hSlices; i++){
            double angle = (pi/2) * i/(double)hSlices;
            double r = getRadius() * cos(angle);
            double z = getRadius() * sin(angle);

            for(int j = 0; j <= vSlices; j++){
                angle = (2*pi) * j/(double) vSlices;
                points[i + hSlices][j].x = r * cos(angle) + reference_point.x;
                points[i + hSlices][j].y = r * sin(angle) + reference_point.y;
                points[i + hSlices][j].z = z + reference_point.z;
            }
        }

        glColor3f(color[0], color[1], color[2]);
        for(int i = 0; i < (int)points.size() - 1; i++){
            for(int j = 0; j < (int)points[i].size() - 1; j++){
                // int white = (j%2);
                // glColor3f( white , white, white);
                glBegin(GL_QUADS);{
                    glVertex3f(points[i][j].x , points[i][j].y , points[i][j].z);				
                    glVertex3f(points[i][j+1].x , points[i][j+1].y , points[i][j+1].z);				
                    glVertex3f(points[i+1][j+1].x , points[i+1][j+1].y , points[i+1][j+1].z);				
                    glVertex3f(points[i+1][j].x , points[i+1][j].y , points[i+1][j].z);				
                }glEnd();			
            }
        }	
    }

    double intersect(Ray ray, vector<double> &cols, int level){
        Vector3D r0 = ray.start.add( reference_point.multiply(-1) );
        Vector3D rd = ray.dir;
        double a = 1;
        double b = 2*rd.dotProduct(r0);
        double c = r0.dotProduct(r0) - getRadius()*getRadius();
        double d = b*b - 4*a*c;

        for(int i = 0; i < 3; i++) cols[i] = color[i];
        if(d < 0){
            return -1;
        }
        d = sqrt(d);
        return (-b - d)/(2.0*a);
    }

};

class Floor : public Object{
public:
    Floor(double floorWidth, double tileWidth) 
        : Object(Vector3D(-floorWidth/2, -floorWidth/2, 0), 
            floorWidth, floorWidth, tileWidth, {0,0,0}, {0,0,0,0}, 1)
    {}

    void draw(){
        glColor3f(1, 1, 1);
        int x = width/length;
        int y = width/length;
        // glBegin(GL_QUADS);{
        //     glVertex3f(-width/2 , -width/2 , 0);				
        //     glVertex3f(-width/2 , width/2 , 0);				
        //     glVertex3f(width/2 , width/2 , 0);				
        //     glVertex3f(width/2 , -width/2 , 0);				
        // }glEnd();
        for(int i = 0; i < x; i++){
            for(int j = 0; j < y; j++){
                int white = ( (i + j) %2);
                double minx = reference_point.x + length*i;
                double miny = reference_point.y + length*j;
                glColor3f( white , white, white);
                glBegin(GL_QUADS);{
                    glVertex3f( minx, miny , reference_point.z);				
                    glVertex3f( minx + length, miny , reference_point.z);
                    glVertex3f( minx + length, miny + length , reference_point.z);
                    glVertex3f( minx, miny + length , reference_point.z);
                }glEnd();
            }
        }
    }

    double intersect(Ray ray, vector<double> &cols, int level){
        // Vector3D r0 = ray.start.add( reference_point.multiply(-1) );
        // Vector3D rd = ray.dir;
        // double a = 1;
        // double b = 2*rd.dotProduct(r0);
        // double c = r0.dotProduct(r0) - getRadius()*getRadius();
        // double d = b*b - 4*a*c;

        // for(int i = 0; i < 3; i++) cols[i] = color[i];
        // if(d < 0){
        //     return -1;
        // }
        // d = sqrt(d);
        // return (-b - d)/(2.0*a);
    }

};
