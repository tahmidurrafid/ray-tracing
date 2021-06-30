#include <windows.h>
#include <GL/glut.h>
#include <bits/stdc++.h>

#define pi (2*acos(0.0))
using namespace std;

int state = 0;

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

	Vector3D cross(const Vector3D &a){
		return Vector3D(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	}
    double dotProduct(const Vector3D &a){
        return x*a.x + y*a.y + z*a.z; 
    }

	Vector3D add(const Vector3D &a){
		return Vector3D(x + a.x, y + a.y, z + a.z);
	}

	Vector3D multiply(double val){
		return Vector3D(val*x, val*y, val*z);
	}

	void shiftInDirection(const Vector3D &to, float steps){
		x += (to.x)*steps;
		y += (to.y)*steps;
		z += (to.z)*steps;
	}

	void rotate(const Vector3D &per, float angle){
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
    double color[3];
    vector<double> coEfficients;
    int shine;

    Object(Vector3D rp, double h, double w, double l, vector<double>  col, vector<double>  coEff, double sh){
        reference_point = rp;
        height = h, width = w, length = l;
        color[0] = col[0];
        color[1] = col[1];
        color[2] = col[2];
        coEfficients = coEff;
        shine = sh;
    }

    virtual void draw(){
        cout << "Draw method not overriddent\n";
    }

    virtual double intersect(Ray &ray, double cols[], int level){
        return -1;
    }

    void calculateColor(Vector3D &point, Vector3D &V, Vector3D &N, double cols[], double color[],int level);

    void setColor(){
    }

    void setShine(){
    }

    void setCoefficient(){
    }
};


vector<Light> lights;
vector<Object*> objects;

void Object::calculateColor(Vector3D &point, Vector3D &V, Vector3D &N, double cols[], double color[], int level){
    for(int i = 0; i < 3; i++) cols[i] = color[i]*coEfficients[AMB];

    for(int i = 0; i < (int)lights.size(); i++){
        Light &light = lights[i];
        // Vector3D L = point.add(light.light_pos.multiply(-1));
        Vector3D L = Vector3D(point.x - light.light_pos.x, point.y - light.light_pos.y, 
                    point.z - light.light_pos.z);
        L.normalize();
        // Vector3D R = N.multiply( - 2 * L.dotProduct(N)).add( L );
        double dot = - 2 * L.dotProduct(N);
        Vector3D R = Vector3D(
            N.x*dot + L.x, N.y*dot + L.y, N.z*dot + L.z
        );

        R.normalize();
        for(int j = 0; j < 3; j++){
            cols[j] += light.color[j]*coEfficients[DIFF]*max(0.0, -L.dotProduct(N) )*color[j];
            cols[j] += light.color[j]*coEfficients[SPEC]*max(pow( -R.dotProduct(V), shine ), 0.0)*color[j];
        }
    }

    double dot = -2* N.dotProduct(V);
    // Vector3D recurDir = N.multiply( dot ).add(V).normalize();
    Vector3D recurDir = Vector3D(N.x*dot + V.x, N.y*dot + V.y, N.z*dot + V.z);
    recurDir.normalize();

    Ray recur = Ray(point.add(recurDir.multiply(.001)), recurDir );

    double reColor[] = {0, 0 ,0};

    double colorReflected[3];
    double t_min = 111111;

    for(int i = 0; i < (int)objects.size(); i++){
        double t = objects[i]->intersect(recur, reColor, level - 1);
        if(t > 0 && t < t_min){
            t_min = t;
            colorReflected[0] = reColor[0];
            colorReflected[1] = reColor[1];
            colorReflected[2] = reColor[2];
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

    double intersect(Ray &ray, double cols[], int level){
        Vector3D r0 = ray.start.add( reference_point.multiply(-1) );
        Vector3D &rd = ray.dir;
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
        Vector3D &V = rd;

        calculateColor(point, V, N, cols, color,level);
        
        return t;
    }

};

class Floor : public Object{
public:
    Floor(double floorWidth, double tileWidth) 
        : Object(Vector3D(-floorWidth/2, -floorWidth/2, 0), 
            floorWidth, floorWidth, tileWidth, {0,0,0}, {.2,.5, .5, .5}, 10)
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

    double intersect(Ray &ray, double cols[], int level){
        Vector3D &r0 = ray.start;
        Vector3D &rd = ray.dir;
        if(rd.z == reference_point.z)
            return -1;
        double t = - r0.z/rd.z;
        // Vector3D point = r0.add(rd.multiply(t));
        Vector3D point = Vector3D(r0.x + rd.x*t, r0.y + rd.y*t, r0.z + rd.z*t);

        if(point.x >= reference_point.x && point.x <= reference_point.x + width && 
            point.y >= reference_point.y && point.y <= reference_point.y + width )
        {
            int u = (point.x - reference_point.x)/length;
            int v = (point.y - reference_point.y)/length;
            int parity = (u + v)%2;
            double color[] = {(double)parity, (double)parity, (double)parity};
            if(level == 0)
                return t;

            Vector3D N = Vector3D(0, 0, 1);
            Vector3D &V = rd;
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

    double determinant(double mat[][3]){
        return mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2]) -
                mat[0][1]*(mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0]) +
                mat[0][2]*(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
    }

    double intersect(Ray &ray, double cols[], int level){
        Vector3D p0 = ray.start.add(points[0].multiply(-1));
        Vector3D v1 = points[1].add(points[0].multiply(-1));
        Vector3D v2 = points[2].add(points[0].multiply(-1));
        Vector3D v3 = ray.dir.multiply(-1);
        // vector<vector<double>> dt = vector<vector<double>>(3, vector<double>(3, 0)), a1, a2, a3;
        double dt [3][3], a1[3][3], a2[3][3], a3[3][3];


        dt[0][0] = v1.x, dt[0][1] = v2.x, dt[0][2] = v3.x;
        dt[1][0] = v1.y, dt[1][1] = v2.y, dt[1][2] = v3.y;
        dt[2][0] = v1.z, dt[2][1] = v2.z, dt[2][2] = v3.z;

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                a1[i][j] = a2[i][j] = a3[i][j] = dt[i][j];  
            }
        }

        a1[0][0] = p0.x, a1[1][0] = p0.y, a1[2][0] = p0.z;
        a2[0][1] = p0.x, a2[1][1] = p0.y, a2[2][1] = p0.z;
        a3[0][2] = p0.x, a3[1][2] = p0.y, a3[2][2] = p0.z;
        double Ddt = determinant(dt), 
                Da1 = determinant(a1), Da2 = determinant(a2), Da3 = determinant(a3);
        if(fabs(Ddt) < .0000001){
            return -1;
        }
        double alpha = Da1/Ddt, beta = Da2/Ddt, t = Da3/Ddt;
        Vector3D point = ray.start.add(ray.dir.multiply(t));
        if(!(alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && alpha + beta <= 1)){
            return -1;
        }

        if(level == 0){
            return t;
        }
        Vector3D N = v1.cross(v2);
        if(ray.dir.dotProduct(N) > 0){
            N = N.multiply(-1);
        }
        N.normalize();
        Vector3D V = ray.dir;

        calculateColor(point, V, N, cols, color,level);
        
        return t;

    }
};


class Quadric : public Object{
public:
    vector<Vector3D> points;
    double A, B, C, D, E, F, G, H, I, J;

    Quadric(double a, double b,double c,double d,double e, double f, double g, double h, double i, double j,
    Vector3D point, double len, double wid, double hei,
    vector<double> col, vector<double> coeff, double sh) 
        : Object(Vector3D(0, 0, 0), 
            hei, wid, len, col, coeff, sh)
    {
        A = a;
        B = b;
        C = c;
        D = d;
        E = e;
        F = f;
        G = g;
        H = h;
        I = i;
        J = j;        
    }

    void draw(){

    }

    double determinant(double mat[][3]){
        return mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2]) -
                mat[0][1]*(mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0]) +
                mat[0][2]*(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
    }

    bool isInside(Vector3D &point){
        if(fabs(length) >= .0000001 && !(point.x >= reference_point.x && point.x <= reference_point.x + length) ){
            return false;            
        }
        if(fabs(width) >= .0000001 && !(point.y >= reference_point.y && point.y <= reference_point.y + width) ){
            return false;
        }
        if(fabs(height) >= .0000001 && !(point.z >= reference_point.z && point.z <= reference_point.z + height) ){
            return false;
        }
        return true;
    }

    double intersect(Ray &ray, double cols[], int level){
        Vector3D &R = ray.start;
        Vector3D &d = ray.dir;
        double a = A*d.x*d.x + B*d.y*d.y + C*d.z*d.z + D*d.x*d.y + E*d.x*d.z + F*d.y*d.z;
        double b = 2*( A*R.x*d.x + B*R.y*d.y + C*R.z*d.z ) + D*R.x*d.y + D*R.y*d.x + E*R.x*d.z + E*R.z*d.x +
                F*R.y*d.z + F*R.z*d.y + G*d.x + H*d.y + I*d.z;
        double c = A*R.x*R.x + B*R.y*R.y + C*R.z*R.z + D*R.x*R.y + E*R.x*R.z + F*R.y*R.z + 
                G*R.x + H*R.y + I*R.z + J;
        double discriminant = b*b - 4*a*c;
        if(discriminant < 0){
            return -1;
        }
        double t1 = (-b - sqrt(discriminant))/(2.0*a);
        double t2 = (-b + sqrt(discriminant))/(2.0*a);
        Vector3D point1 = R.add(d.multiply(t1));
        Vector3D point2 = R.add(d.multiply(t2));
        double t = -1;
        Vector3D point;
        if(t1 >= 0 && isInside(point1)){
            t = t1;
            point = point1;
        }else if(t2 >= 0 && isInside(point2)){
            t = t2;
            point = point2;
        }
        if(t < 0) 
            return -1;
        if(level == 0){
            return t;
        }
        cols[0] = cols[1] = cols[2] = .3;
        Vector3D N(2*A*point.x + D*point.y + E*point.z + G,
                2*B*point.y + D*point.x + F*point.z + H,
                2*C*point.z + E*point.x + F*point.y + I);
        
        N.normalize();
        if(N.dotProduct(ray.dir) > 0){
            N = N.multiply(-1);
        }
        Vector3D &V = ray.dir;
        state = 1;
        calculateColor(point, V, N, cols, color,level);        
        state = 0;
        return t;
        // if(level == 0){
        //     return t;
        // }
        // Vector3D N = v1.cross(v2);
        // if(ray.dir.dotProduct(N) > 0){
        //     N = N.multiply(-1);
        // }
        // N.normalize();
        // Vector3D V = ray.dir;

        // calculateColor(point, V, N, cols, color,level);
        
        // return t;

    }
};
