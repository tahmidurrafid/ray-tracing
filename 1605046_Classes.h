#include <windows.h>
#include <GL/glut.h>
#include <bits/stdc++.h>

#define pi (2*acos(0.0))
using namespace std;

enum COEFF{AMB = 0, DIFF = 1, SPEC = 2, RECUR = 3};

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
        return -1;
    }

    void calculateColor(Vector3D point, Vector3D V, Vector3D N, vector<double> &cols, vector<double> color,int level);

    void setColor(){
    }

    void setShine(){
    }

    void setCoefficient(){
    }
};


vector<Light> lights;
vector<Object*> objects;

void Object::calculateColor(Vector3D point, Vector3D V, Vector3D N, vector<double> &cols,vector<double> color, int level){
    for(int i = 0; i < 3; i++) cols[i] = color[i]*coEfficients[AMB];

    for(int i = 0; i < (int)lights.size(); i++){
        Light light = lights[i];
        Vector3D L = point.add(light.light_pos.multiply(-1));
        L.normalize();
        Vector3D R = N.multiply( - 2 * L.dotProduct(N)).add( L );
        R.normalize();
        for(int j = 0; j < 3; j++){
            cols[j] += light.color[j]*coEfficients[DIFF]*max(0.0, -L.dotProduct(N) )*color[j];
            cols[j] += light.color[j]*coEfficients[SPEC]*max(pow( -R.dotProduct(V), shine ), 0.0)*color[j];
        }
    }

    Vector3D recurDir = N.multiply( -2* N.dotProduct(V) ).add(V).normalize();

    Ray recur = Ray(point.add(recurDir.multiply(.001)), recurDir );
    // Ray recur = Ray(point , recurDir );

    vector<double> reColor(3, 0);
    vector<double> colorReflected;
    double t_min = 111111;

    for(int i = 0; i < (int)objects.size(); i++){
        double t = objects[i]->intersect(recur, reColor, level - 1);
        if(t > 0 && t < t_min){
            t_min = t;
            colorReflected = reColor;
        }
    }
    if(t_min < 10000){
        for(int i = 0; i < 3; i++){
            cols[i] += colorReflected[i] * coEfficients[RECUR];
        }
    }
}


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
        for(int i = 0; i < 3; i++) cols[i] = 0;
        if(d < 0){
            return -1;
        }
        d = sqrt(d);

        double t = (-b - d)/(2.0*a);
        if(level == 0){
            return t;
        }
        Vector3D point = ray.start.add( ray.dir.multiply(t) );
        Vector3D N = point.add(reference_point.multiply(-1));
        N.normalize();
        Vector3D V = rd;

        calculateColor(point, V, N, cols, color,level);
        
        return t;
    }

};

class Floor : public Object{
public:
    Floor(double floorWidth, double tileWidth) 
        : Object(Vector3D(-floorWidth/2, -floorWidth/2, 0), 
            floorWidth, floorWidth, tileWidth, {0,0,0}, {.5,.2, .2, .2}, 10)
    {}

    void draw(){
        glColor3f(1, 1, 1);
        int x = width/length;
        int y = width/length;

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
        Vector3D r0 = ray.start;
        Vector3D rd = ray.dir;
        if(rd.z == reference_point.z)
            return -1;
        double t = - r0.z/rd.z;
        Vector3D point = r0.add(rd.multiply(t));
        if(point.x >= reference_point.x && point.x <= reference_point.x + width && 
            point.y >= reference_point.y && point.y <= reference_point.y + width )
        {
            int u = (point.x - reference_point.x)/length;
            int v = (point.y - reference_point.y)/length;
            int parity = (u + v)%2;
            vector<double> color = {(double)parity, (double)parity, (double)parity};
            if(level == 0)
                return t;

            Vector3D N = Vector3D(0, 0, 1);
            N.normalize();
            Vector3D V = rd;
            calculateColor(point, V, N, cols, color, level );
            return t;
        }
        return -1;
    }
};

class Triangle : public Object{
public:

    vector<Vector3D> points;

    Triangle(vector<Vector3D> p, vector<double> col, vector<double> coeff, double sh) 
        : Object(Vector3D(0, 0, 0), 
            0, 0, 0, col, coeff, sh)
    {
        points = p;
    }

    void draw(){
        glColor3f( color[0] ,  color[1],  color[2]);
        glBegin(GL_TRIANGLES);
        {
            for(int i = 0; i < 3; i++){
                glVertex3f(points[i].x, points[i].y , points[i].z);                
            }
        }glEnd();
    }

    double intersect(Ray ray, vector<double> &cols, int level){
        return -1;
        // Vector3D r0 = ray.start;
        // Vector3D rd = ray.dir;
        // if(rd.z == reference_point.z)
        //     return -1;
        // double t = - r0.z/rd.z;
        // Vector3D point = r0.add(rd.multiply(t));
        // if(point.x >= reference_point.x && point.x <= reference_point.x + width && 
        //     point.y >= reference_point.y && point.y <= reference_point.y + width )
        // {
        //     int u = (point.x - reference_point.x)/length;
        //     int v = (point.y - reference_point.y)/length;
        //     int parity = (u + v)%2;
        //     vector<double> color = {(double)parity, (double)parity, (double)parity};
        //     if(level == 0)
        //         return t;

        //     Vector3D N = Vector3D(0, 0, 1);
        //     N.normalize();
        //     Vector3D V = rd;
        //     calculateColor(point, V, N, cols, color, level );
        //     return t;
        // }
        // return -1;
    }
};
