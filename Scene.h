#ifndef SCENE_H
#define SCENE_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;

struct Color {
  double red;
  double green;
  double blue;
};

struct Shading {
  Shading();
  Shading(double Kd_input, double Ks_input, double shine_input, double transmittance_input, double index_of_refraction_input);
  double Kd;
  double Ks;
  double shine;
  double transmittance;
  double index_of_refraction;
};

struct HitRecord {
  float t;
  Color c;
  Eigen::Vector3d normal;
  Shading shading;
};

struct Light {
  Light(Eigen::Vector3d inputPosition, double red=1, double green=1, double blue =1);
  Eigen::Vector3d position;
  Color color;
};

struct Fragment {
  int i;
  int j;
  double z;
  Color color;
};

class Camera {
  public:
  Camera();
  Camera(Eigen::Vector3d fromInput, Eigen::Vector3d atInput, Eigen::Vector3d upInput, double angleInput, double hitherInput, int resxInput, int resyInput);
  Eigen::Vector3d from;
  Eigen::Vector3d at;
  Eigen::Vector3d up;
  double angle;
  double hither;
  int resx;
  int resy;
  double left;
  double right;
  double top;
  double bottom;
  double d;
  double aspectRatio;
  Eigen::Vector3d u;
  Eigen::Vector3d v;
  Eigen::Vector3d w;
};

struct Ray {
  Ray();
  Ray(double i, double j, Camera camera);
  Eigen::Vector3d origin; //eye
  Eigen::Vector3d direction; //direction
};

class Shape {
  public:
  Shape();
};

class Triangle {
  public:
  Triangle(); //compute vertex colors in constructor
  Triangle(Eigen::Vector3d firstVertex, Eigen::Vector3d secondVertex, Eigen::Vector3d thirdVertex, bool hasVertexNormals);
  void transformVertices(Eigen::Matrix<double, 4, 4> M);
  std::vector<Fragment> rasterize(vector<vector<array<double, 3>>> &pixels, std::vector<Light> lights, Camera camera, Shading shading, Color color); //pass pixels by reference so it can be updated
  friend ostream& operator<<(ostream& os, const Triangle& triangle);
  bool intersect(const Ray &r, double t0, double t1, HitRecord *hr);
  void shadeVertices(std::vector<Light> lights, Camera camera, Color color, Shading shading);
  Color shadeVertex(Eigen::Vector3d vector, Eigen::Vector3d normal, std::vector<Light> lights, Camera camera, Shading shading, Color color);
  double getCameraDistance(Camera camera);
  Eigen::Vector3d getLightDirection(Light lights, Eigen::Vector3d vertex);
  Eigen::Vector3d getViewDirection(Eigen::Vector3d from, Eigen::Vector3d vertex);
  Eigen::Vector3d vertex1;
  Eigen::Vector3d vertex2;
  Eigen::Vector3d vertex3;
  Eigen::Vector3d vertex1Normal;
  Eigen::Vector3d vertex2Normal;
  Eigen::Vector3d vertex3Normal;
  Eigen::Vector3d surface_normal;
  Eigen::Vector3d view_direction;
  Eigen::Vector3d transV1; //transformed vertices 
  Eigen::Vector3d transV2;
  Eigen::Vector3d transV3;
  Color vertex1Color;
  Color vertex2Color;
  Color vertex3Color;
  double cameraDistance(Camera camera);
};

class Polygon: public Shape {
  public:
  Polygon(std::vector<Eigen::Vector3d> inputVertices);
  Polygon(std::vector<Eigen::Vector3d> inputVertices, std::vector<Eigen::Vector3d> inputNormals);
  void transformVertices(Eigen::Matrix<double, 4, 4> M);
  std::vector<Fragment> rasterize(std::vector<std::vector<array<double, 3>>> &pixels, std::vector<Light> lights, Camera camera);
  void shadeVertices(std::vector<Light> lights, Camera camera);
  void getCameraDistance(Camera camera);
  bool intersect(const Ray &r, double t0, double t1, HitRecord *hr);
  //write a destructor
  double cameraDistance;
  std::vector<Eigen::Vector3d> vertices;
  std::vector<Eigen::Vector3d> normals;
  std::vector<Triangle> triangles;
  friend ostream& operator<<(ostream& os, const Polygon& polygon);
  Color color;
  Shading shading;
};

class Scene {
  public:
  Scene(string filename);
  void calcIntersections();
  void processVertices();
  void rasterize();
  void WriteToPPM();
  friend ostream& operator<<(ostream& os, const Scene& scene);
  Camera camera;
  std::vector<Light> lights;
  std::vector<Polygon> shapes;
  Color backgroundColor;
  Eigen::Vector3d xDir;
  Eigen::Vector3d yDir;
  Eigen::Vector3d zDir;
  vector<Fragment> fragments;
  vector<vector<array<double, 3>>> pixels;
  //vector<int[3]> pixels;
};



#endif
