#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265358979323846

#include <iostream>
#include <random>

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);

class Vector{
public:
	explicit Vector(double x = 0., double y = 0., double z = 0.){
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
	};
	double operator[](int i) const {return coords[i];};
	double &operator[](int i){return coords[i];};
	double sqrNorm() const {
		return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
	}
	Vector get_normalized(){
		double n = sqrt(sqrNorm());
		return Vector(coords[0]/n, coords[1]/n, coords[2]/n);
	}
	Vector gamma_correct(){
		double gamma = 2.2;
		return Vector(
			std::pow(coords[0]/255, 1./gamma)*255, 
			std::pow(coords[1]/255, 1./gamma)*255, 
			std::pow(coords[2]/255, 1./gamma)*255
		);
	}
private:
	double coords[3];
};

Vector operator+(const Vector &a, const Vector &b){
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector &a, const Vector &b){
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(double a, const Vector &b){
	return Vector(a*b[0], a*b[1], a*b[2]);
}

Vector operator*(const Vector &a, double b){
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}

Vector operator/(const Vector &a, double b){
	return Vector(a[0]/b, a[1]/b, a[2]/b);
}

double dot(const Vector &a, const Vector &b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

class Ray{
public:
	Ray(const Vector &u, const Vector &C): u(u), C(C) {
	}
	Vector u, C;
};

class Sphere{
public:
	Sphere(const Vector &O, double R, const Vector &sphere_color): O(O), R(R), sphere_color(sphere_color) {
	}
	bool intersect(const Ray &r, Vector &P, Vector &N, double &t) {
		// solves a*t^2 + b*t + c = 0
		double a = 1;
		double b = 2*dot(r.u, r.C - O);
		double c = (r.C - O).sqrNorm() - R*R;
		double delta = b*b - 4*a*c;

		if (delta < 0)
        {
            return false;
        }
        double sqDelta = sqrt(delta);
        double t2 = (-b + sqDelta) / (2*a);

        if (t2 < 0)
            return false;

        // double t;
        double t1 = (-b - sqDelta) / (2*a);

        if (t1 > 0)
        {
            t = t1;
        }
        else
        {
            t = t2;
        }

        P = r.C + t*r.u;
        N = (P - O).get_normalized();

        return true;
	}
	Vector O;
	double R;
	Vector sphere_color;
};

class Scene{
public:
	Scene(){};
	std::vector<Sphere>objects;
	bool intersect(const Ray &r, Vector &P, Vector &N, Vector &sphere_color){
		double t = 1E10;
		bool has_inter = false;
		for (int i = 0; i < objects.size(); i++){
			Vector localP, localN;
			double localt;
			if (objects[i].intersect(r, localP, localN, localt) && localt < t){
				t = localt;
				has_inter = true;
				sphere_color = objects[i].sphere_color;
				P = localP;
				N = localN;
			}
		}
		return has_inter;
	};
};

int main() {
	int W = 512;
	int H = 512;

	//Vector O(0, 0, 0);
	//double R = 10;
	Vector C(0, 0, 55);
	Scene scene;
	
	Sphere S1(Vector(-10, 0, 0), 10, Vector(1., 0., 0.));
    Sphere S2(Vector(10, 0, 0), 10, Vector(0., 0., 1.));
    Sphere Ssol(Vector(0, -1000, 0), 990, Vector(0., 1., 1.));
	Sphere SmurBehind(Vector(0, 0, 1000), 940, Vector(1., 0., 1.));
	Sphere SmurFront(Vector(0, 0, -1000), 940, Vector(0., 1., 0.));
	Sphere Splafond(Vector(0, 1000, 0), 940, Vector(1., 0., 0.));

    scene.objects.push_back(S1);
    scene.objects.push_back(S2);
    scene.objects.push_back(Ssol);
    scene.objects.push_back(SmurBehind);
    scene.objects.push_back(SmurFront);
    scene.objects.push_back(Splafond);

	double fov = 60* M_PI / 180;
	double I = 1E7;
	//Vector rho(1,0,0);
	Vector Light(-10,20,40);
	
	std::vector<unsigned char> image(W*H*3, 0);

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector u(j - W/2, i - H/2, -W/(2.*tan(fov/2)));
			u = u.get_normalized();
			Ray r(u, C);
			Vector P, N, sphere_color;
			bool inter = scene.intersect(r, P, N, sphere_color);
			Vector color(0, 0, 0);
			if (inter) {
				Vector PLight = Light - P;
				double d = sqrt(PLight.sqrNorm());
				// color = Vector(255, 255, 255);
				color = I/(4*M_PI*d*d) * (sphere_color/M_PI) * std::max(0., dot(N, PLight/d));
				color.gamma_correct();
			}

			image[((H - i - 1)*W + j)*3 + 0] = color[0];
			image[((H - i - 1)*W + j)*3 + 1] = color[1];
			image[((H - i - 1)*W + j)*3 + 2] = color[2];
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}