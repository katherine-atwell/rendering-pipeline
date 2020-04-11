#include "Scene.h"
#include <cmath>
#include "math.h"

using namespace std;

Light::Light(Eigen::Vector3d inputPosition, double red, double green, double blue) {
  position = inputPosition;
  color.red = red;
  color.green = green;
  color.blue = blue;
}

Ray::Ray() {}

//sets origin and direction based on pixel's location in image
Ray::Ray(double i, double j, Camera camera) {
  double d = 1;
  double h = (camera.top - camera.bottom);
  double x = camera.left + (j * h + h/2)/camera.resx*d; 
  double y = camera.top - (i * h + h/2)/camera.resy*d;
  double z = -1;
  Eigen::Vector3d s = x*camera.u + y*camera.v + z*camera.w;
  origin = camera.from;
  direction = (s).normalized();
}

Camera::Camera() {}

//sets camera parameters and calculates vectors u, v, and w
Camera::Camera(Eigen::Vector3d fromInput, Eigen::Vector3d atInput, Eigen::Vector3d upInput, double angleInput, double hitherInput, int resxInput, int resyInput) {
  from = fromInput;
  at = atInput;
  up = upInput;
  angle = angleInput;
  hither = hitherInput;
  resx = resxInput;
  resy = resyInput;
  d = 1;
  aspectRatio = resx/resy;
  left = -tan(angle/2)/d;
  right = tan(angle/2)/d;
  bottom = -tan(angle/2)/d/aspectRatio;
  top = tan(angle/2)/d/aspectRatio;
  w = (from - at).normalized();
  u = (up.cross(w)).normalized();
  v = (w.cross(u)).normalized();
}

Shading::Shading() {}

Shading::Shading(double Kd_input, double Ks_input, double shine_input, double transmittance_input, double index_of_refraction_input) {
  Kd = Kd_input;
  Ks = Ks_input;
  shine = shine_input;
  transmittance = transmittance_input;
  index_of_refraction = index_of_refraction_input;
}

//parses file based on passed in filename
Scene::Scene(string filename) {
  //parse input from file whose name is given as argument
  ifstream inputFile;
  inputFile.open(filename);
  string backgroundColorChar;
  getline(inputFile, backgroundColorChar, ' ');
  inputFile >> backgroundColor.red;
  inputFile >> backgroundColor.green;
  inputFile >> backgroundColor.blue;
  inputFile.ignore(200, '\n');
  
  string viewingChar;
  getline(inputFile, viewingChar, '\n');
    
  string fromString;
  getline(inputFile, fromString, ' ');
  double x, y, z;
  inputFile >> x;
  inputFile >> y;
  inputFile >> z;
  Eigen::Vector3d from;
  from = Eigen::Vector3d(x, y, z);

  string atString;
  
  getline(inputFile, atString, ' ');
  
  inputFile >> x;
  inputFile >> y;
  inputFile >> z;

  Eigen::Vector3d at;
  at = Eigen::Vector3d(x, y, z);
  string upString;
  getline(inputFile, upString, ' ');
  
  inputFile >> x;
  inputFile >> y;
  inputFile >> z;
  Eigen::Vector3d up;
  up = Eigen::Vector3d(x, y, z);
    
  string angleString;
  getline(inputFile, angleString, ' ');
  
  double angleInput;
  double angle;
  inputFile >> angleInput;
  angle = angleInput * M_PI / 180;

  string hitherString;
  getline(inputFile, hitherString, ' ');
  
  double hither;
  inputFile >> hither;

  string resString;
  getline(inputFile, resString, ' ');
  
  int resx, resy;
  inputFile >> resx;
  inputFile >> resy;
  inputFile.ignore(200, '\n');
  string shapeString;
  getline(inputFile, shapeString, ' ');
  double lightX, lightY, lightZ;
  
  while (shapeString == "l") {
    inputFile >> lightX;
    inputFile >> lightY;
    inputFile >> lightZ;
    // TODO: support colors other than white
    inputFile.ignore(200, '\n');
    Eigen::Vector3d lightPosition(lightX, lightY, lightZ);
    Light newLight(lightPosition);
    lights.push_back(newLight);
    getline(inputFile, shapeString, ' ');
  } 
  
  /*for (int i=0; i < (int) lights.size(); i++) {
    cout << lights[i].position << endl;
  }*/
  
  Color currColor;
  inputFile >> currColor.red;
  inputFile >> currColor.green;
  inputFile >> currColor.blue;
  Shading currShading;
  inputFile >> currShading.Kd;
  inputFile >> currShading.Ks;
  inputFile >> currShading.shine;
  inputFile >> currShading.transmittance;
  inputFile >> currShading.index_of_refraction;
  inputFile.ignore(200, '\n');
  getline(inputFile, shapeString, ' ');
  Camera newCamera(from, at, up, angle, hither, resx, resy);
  camera = newCamera;
  
  while ((shapeString == "p" || shapeString == "s" || shapeString == "f" || shapeString == "pp") && !inputFile.eof()) {
    if (shapeString == "p") {
      int numVertices;
      inputFile >> numVertices;
      std::vector<Eigen::Vector3d> vertices;
      for (int i=0; i < numVertices; i++) {
        double x, y, z;
        inputFile >> x;
        inputFile >> y;
        inputFile >> z;
        Eigen::Vector3d newVertex(x, y, z);
        vertices.push_back(newVertex);
      }
      Polygon newPolygon(vertices);
      newPolygon.color = currColor;
      newPolygon.shading = currShading;
      shapes.push_back(newPolygon);
      inputFile.ignore(200, '\n');
      getline(inputFile, shapeString, ' ');
    }
    else if (shapeString == "f") {
      inputFile >> currColor.red;
      inputFile >> currColor.green;
      inputFile >> currColor.blue;
      double Kd,Ks, shine, transmittance, index_of_refraction;
      inputFile >> Kd;
      inputFile >> Ks;
      inputFile >> shine;
      inputFile >> transmittance;
      inputFile >> index_of_refraction;
      
      inputFile.ignore(200, '\n');
      getline(inputFile, shapeString, ' ');
    }
    else if (shapeString == "pp") {
      int numVertices;
      inputFile >> numVertices;
      std::vector<Eigen::Vector3d> vertices;
      std::vector<Eigen::Vector3d> normals;
      for (int i=0; i < numVertices; i++) {
        double normX, normY, normZ;
        double x, y, z;
        inputFile >> x;
        inputFile >> y;
        inputFile >> z;
        inputFile >> normX;
        inputFile >> normY;
        inputFile >> normZ;
        Eigen::Vector3d newVertex(x, y, z);
        vertices.push_back(newVertex);
        Eigen::Vector3d newNormal(normX, normY, normZ);
        normals.push_back(newNormal);
      }
      
      Polygon newPolygon(vertices, normals);
      newPolygon.color = currColor;
      newPolygon.shading = currShading;
      shapes.push_back(newPolygon);
      inputFile.ignore(200, '\n');
      getline(inputFile, shapeString, ' ');
    }
  }
  
  //populates pixels array with background color
  for (int i=0; i<camera.resy; i++) {
    vector<array<double,3>> row;
    for (int j=0; j<camera.resx; j++)  {
      //remember to bounds check later; look at project 2
      row.push_back({backgroundColor.red * 255, backgroundColor.green * 255, backgroundColor.blue* 255});
    }
    pixels.push_back(row);
  }
}

void Scene::processVertices() {
  Eigen::Matrix<double, 4, 4> Mvp;
  Eigen::Matrix<double, 4, 4> Mper;
  Eigen::Matrix<double, 4, 4> cam1;
  Eigen::Matrix<double, 4, 4> cam2;
  Eigen::Matrix<double, 4, 4> Mcam;
  Eigen::Matrix<double, 4, 4> M;
  
  double n = (double) camera.hither; //n is z coordinate of from I think? Will have to confirm
  double f = (double) camera.hither * 100;   //f is z coordinate of at??
  double r = camera.right;
  double l = camera.left;
  double t = camera.top;
  double b = camera.bottom;
  Eigen::Vector3d u = camera.u;
  Eigen::Vector3d v = camera.v;
  Eigen::Vector3d w = camera.w;
  Eigen::Vector3d e = camera.from;
  
  Mvp << camera.resx/2, 0, 0, (camera.resx-1)/2,
         0, camera.resy/2, 0, (camera.resy-1)/2, 
         0, 0, 1, 0, 
         0, 0, 0, 1;
  
  Mper << -2*n/(r - l), 0, -(r+l)/(r-l), 0,
          0, -2*n/(t-b), (t+b)/(t-b), 0,  
          0, 0, (f+n)/(n-f), -2*f*n/(f-n),
          0, 0, 1, 0;
  
  cam1 << u.x(), u.y(), u.z(), 0,
          v.x(), v.y(), v.z(), 0, 
          w.x(), w.y(), w.z(), 0, 
          0, 0, 0, 1;
  cam2 << 1, 0, 0, -1*e.x(),
          0, 1, 0, -1*e.y(),
          0, 0, 1, -1*e.z(),
          0, 0, 0, 1;
 
  Mcam = cam1 * cam2;
  M = Mvp * Mper * Mcam;
  
  for (int i=0; i < (int) shapes.size(); i++) {
    shapes[i].transformVertices(M);
  } 

  for (int i=0; i < (int) shapes.size(); i++) {
    shapes[i].shadeVertices(lights, camera);
  }
}

//calls rasterize() on all polygons contained within image
void Scene::rasterize() {
  processVertices();

  //loops through from end (greatest distance from camera) to beginning to update color of pixels
  for (int i=0; i < (int) shapes.size(); i++) {
    vector<Fragment> polygonFragments;
    polygonFragments = shapes[i].rasterize(pixels, lights, camera);
    for (int j=0; j < (int) polygonFragments.size(); j++) {
      fragments.push_back(polygonFragments[j]);
    }
  }
  
  //sort fragments
  sort(fragments.begin(), fragments.end());
  
  for (int k= 0; k < (int) fragments.size(); k++) {
    pixels[fragments[k].i][fragments[k].j] = {fragments[k].color.red, fragments[k].color.green, fragments[k].color.blue};
  }
}

//determines whether a triangle intersects with a ray
bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord *hr) {
  Eigen::Matrix<double, 3, 3> M;
  M << (r.direction).x(), (vertex1-vertex2).x(), (vertex1-vertex3).x(),
    (r.direction).y(), (vertex1-vertex2).y(), (vertex1-vertex3).y(),
    (r.direction).z(), (vertex1-vertex2).z(), (vertex1-vertex3).z();

  Eigen::Vector3d multipliedMatrix((vertex1-r.origin).x(), (vertex1-r.origin).y(), (vertex1-r.origin).z());
  Eigen::Matrix<double, 3, 3> mInverse;
  mInverse = M.inverse();
  Eigen::Vector3d x;
  x = mInverse * multipliedMatrix;
  double t = x.x();
  double beta = x.y();
  double gamma = x.z();

  if (t < t0 || t > t1) {
    return false;
  }
  if (gamma < 0 || gamma > 1) {
    return false;
  }
  if (beta < 0 || beta + gamma > 1) {
    return false;
  }
  hr->t = t;
  hr->normal = surface_normal;
  return true;
}

//determines whether a convex polygon intersects with a triangle. This is based on whether or not                                                                            
bool Polygon::intersect(const Ray &r, double t0, double t1, HitRecord *hr) {
  bool hit = false;
  for (int i=0; i < (int) triangles.size(); i++) {
    if (triangles[i].intersect(r, t0, t1, hr)) {                                                                                              
      t1 = hr->t;
      hr->c = color;
      hr->shading = shading;
      hit = true;
    }
  }
  return hit;
}


//writes calculated vector of pixels to a PPM file
void Scene::WriteToPPM() {
  unsigned char pixelArray[3];
  FILE *f = fopen("hide.ppm","wb");
  fprintf(f, "P6\n%d %d\n%d\n", camera.resx, camera.resy, 255);
  
  for(int i = 0; i < camera.resy; ++i) {
     for(int j = 0; j < camera.resy; ++j) {
       pixelArray[0] = (unsigned char) pixels[i][j][0];
       pixelArray[1] = (unsigned char) pixels[i][j][1];
       pixelArray[2] = (unsigned char) pixels[i][j][2];
       fwrite(pixelArray, 1, 3, f);
     }
  }
 fclose(f);
}

Shape::Shape() {}

void Triangle::shadeVertices(vector<Light> lights, Camera camera, Color color, Shading shading) {
  vertex1Color.red = shadeVertex(vertex1, vertex1Normal, lights, camera, shading, color).red * 255;
  vertex1Color.green = shadeVertex(vertex1, vertex1Normal, lights, camera, shading, color).green * 255;
  vertex1Color.blue = shadeVertex(vertex1, vertex1Normal, lights, camera, shading, color).blue * 255;
  vertex2Color.red = shadeVertex(vertex2, vertex2Normal, lights, camera, shading, color).red * 255;
  vertex2Color.green = shadeVertex(vertex2, vertex2Normal, lights, camera, shading, color).green * 255;
  vertex2Color.blue = shadeVertex(vertex2, vertex2Normal, lights, camera, shading, color).blue * 255;
  vertex3Color.red = shadeVertex(vertex3, vertex3Normal, lights, camera, shading, color).red * 255;
  vertex3Color.green = shadeVertex(vertex3, vertex3Normal, lights, camera, shading, color).green * 255;
  vertex3Color.blue = shadeVertex(vertex3, vertex3Normal, lights, camera, shading, color).blue * 255;
}

Color Triangle::shadeVertex(Eigen::Vector3d vectorArg, Eigen::Vector3d normal, std::vector<Light> lights, Camera camera, Shading shading, Color color) {

  Color localColor;
  int numLights = lights.size();
  localColor.red = 0;
  localColor.green = 0;
  localColor.blue = 0;
  for(int i=0; i < (int) lights.size(); i++) {
    Eigen::Vector3d light_direction = getLightDirection(lights[i], vectorArg);
    Eigen::Vector3d view_direction = getViewDirection(camera.from, vectorArg);
    Eigen::Vector3d h_vector = (light_direction + view_direction).normalized();
    double lightIntensity = 1/sqrt(numLights);
    double diffuse = max((double) 0, (double) normal.dot(light_direction));
    double specular = pow(max((double) 0, (double) normal.dot(h_vector)), shading.shine);
    localColor.red += (shading.Kd * color.red * diffuse + shading.Ks * specular) * lightIntensity;
    localColor.green += (shading.Kd * color.green * diffuse + shading.Ks * specular) * lightIntensity;
    localColor.blue += (shading.Kd * color.blue * diffuse + shading.Ks * specular) * lightIntensity;
  }
  localColor.red = fmax(min(localColor.red, (double) 1), (double) 0);
  localColor.green = fmax(min(localColor.green, (double) 1), (double) 0);
  localColor.blue = fmax(min(localColor.blue, (double) 1), (double) 0);
  
  return localColor;
}

void Triangle::transformVertices(Eigen::Matrix<double, 4, 4> M) {  
  Eigen::Vector4d p;
  Eigen::Vector4d q;
  Eigen::Vector4d r;
  
  Eigen::Vector4d a;
  Eigen::Vector4d b;
  Eigen::Vector4d c;
  
  // Apply transformation for vertex1 and vertex2
  a << (double) vertex1.x(), (double) vertex1.y(), (double) vertex1.z(), 1;
  b << (double) vertex2.x(), (double) vertex2.y(), (double) vertex2.z(), 1;
  c << (double) vertex3.x(), (double) vertex3.y(), (double) vertex3.z(), 1;
  
  p = M * a;
  q = M * b;
  r = M * c;
  
  transV1 << p.x()/p.w(), p.y()/p.w(), p.z()/p.w(); 
  transV2 << q.x()/q.w(), q.y()/q.w(), q.z()/q.w();
  transV3 << r.x()/r.w(), r.y()/r.w(), r.z()/r.w();
}


//rasterizes triangle image by calculating barycentric coordinates
vector<Fragment> Triangle::rasterize(vector<vector<array<double, 3>>> &pixels, std::vector<Light> lights, Camera camera, Shading shading, Color color) {
  //create bounding box
  double x0 = (double) transV1.x();
  double x1 = (double) transV2.x();
  double x2 = (double) transV3.x();
  double y0 = (double) transV1.y();
  double y1 = (double) transV2.y();
  double y2 = (double) transV3.y();
  
  double xmin = min({x0, x1, x2}); //min x coord
  double xmax = max({x0, x1, x2}); //max x coord
  double ymin = min({y0, y1, y2}); //min y coord
  double ymax = max({y0, y1, y2}); //max y coord
  
  xmin = max(xmin, 0.0);
  xmax = min(xmax, (double) camera.resx);
  ymin = max(ymin, 0.0);
  ymax = min(ymax, (double) camera.resx);
  
  double f01;
  double f12;
  double f20;
  
  double f12naught = ((y1 - y2) * x0) + ((x2 - x1) * y0) + (x1*y2) - (y1*x2);
  double f20naught = ((y2 - y0) * x1) + ((x0 - x2) * y1) + (x2*y0) - (y2*x0);
  double f01naught = ((y0 - y1) * x2) + ((x1 - x0) * y2) + (x0*y1) - (y0*x1);
  
  double alpha, beta, gamma;
  vector<Fragment> fragments;
  for (int i=ymin; i < ymax; i++) {
    for (int j=xmin; j < xmax; j++) {
      f01 = (transV1.y() - transV2.y()) * j + (transV2.x() - transV1.x()) * i + transV1.x()*transV2.y() - transV1.y()*transV2.x();
      f12 = (transV2.y() - transV3.y()) * j + (transV3.x() - transV2.x()) * i + transV2.x()*transV3.y() - transV2.y()*transV3.x();
      f20 = (transV3.y() - transV1.y()) * j + (transV1.x() - transV3.x()) * i + transV3.x()*transV1.y() - transV3.y()*transV1.x();
     
      alpha = f12/f12naught;
      beta = f20/f20naught;
      gamma = f01/f01naught;
     
      if (alpha > 0 && beta > 0 && gamma > 0) {
        std::array<double,3> pixelColor = {color.red, color.green, color.blue};
        Fragment fragment;
        fragment.z = transV1.z();
        fragment.i = camera.resy - 1 - i;
        fragment.j = j;
        
        Eigen::Vector3d currPoint = alpha * vertex1 + beta * vertex2 + gamma * vertex3;
        Eigen::Vector3d normal = alpha * vertex1Normal + beta * vertex2Normal + gamma * vertex3Normal;
        Color localColor;
        int numLights = lights.size();
     
        localColor.red = 0;
        localColor.green = 0;
        localColor.blue = 0;
        for(int i=0; i < (int) lights.size(); i++) {
          Eigen::Vector3d light_direction = getLightDirection(lights[i], currPoint);
          Eigen::Vector3d view_direction = getViewDirection(camera.from, currPoint);
          Eigen::Vector3d h_vector = (light_direction + view_direction).normalized();
          
          double lightIntensity = 1/sqrt(numLights);
          double diffuse = max((double) 0, (double) normal.dot(light_direction));
          double specular = pow(max((double) 0, (double) normal.dot(h_vector)), shading.shine);
          
          localColor.red += (shading.Kd * color.red * diffuse + shading.Ks * specular) * lightIntensity;
          localColor.green += (shading.Kd * color.green * diffuse + shading.Ks * specular) * lightIntensity;
          localColor.blue += (shading.Kd * color.blue * diffuse + shading.Ks * specular) * lightIntensity;
        }
        localColor.red = 255 * fmax(min(localColor.red, (double) 1), (double) 0);
        localColor.green = 255 * fmax(min(localColor.green, (double) 1), (double) 0);
        localColor.blue = 255 * fmax(min(localColor.blue, (double) 1), (double) 0);
        
        fragment.color.red = localColor.red;
        fragment.color.green = localColor.green;
        fragment.color.blue = localColor.blue;
        
        fragments.push_back(fragment);
      }
    }
  }
  
  return fragments;
}

double Triangle::getCameraDistance(Camera camera) {
  //get distance between vertex1 and camera.from; will probably update later but this is the current version
  double xMinusFrom = vertex1.x() - camera.from.x();
  double yMinusFrom = vertex1.y() - camera.from.y();
  double zMinusFrom = vertex1.z() - camera.from.z();
  double cameraDistance = sqrt(xMinusFrom * xMinusFrom + yMinusFrom * yMinusFrom + zMinusFrom * zMinusFrom);
  return cameraDistance;
}

//shades vertices for each triangle contained within polygon (ended up not being needed, as I implemented fragment shading later)
void Polygon::shadeVertices(vector<Light> lights, Camera camera) {
  for (int i=0; i < (int) triangles.size(); i++) {
    triangles[i].shadeVertices(lights, camera, color, shading);
  }
}

//transforms vertices on each of the triangles contained in polygon
void Polygon::transformVertices(Eigen::Matrix<double, 4, 4> M) {
  for (int i=0; i < (int) triangles.size(); i++) {
    triangles[i].transformVertices(M);
  }
}

//calls rasterize on each triangle contained within polygon
vector<Fragment> Polygon::rasterize(vector<vector<array<double, 3>>> &pixels, std::vector<Light> lights, Camera camera) {
  vector<Fragment> fragments;
  for (int i=0; i < (int) triangles.size(); i++) {
    //hopefully pixels updates? If not will pass it another way but fingers crossed  
    vector<Fragment> triangleFragments;
    triangleFragments = triangles[i].rasterize(pixels, lights, camera, shading, color);
    for (int j=0; j < (int) triangleFragments.size(); j++) {
      fragments.push_back(triangleFragments[j]);
    }
  }
  return fragments;
}

//ended up being unnecessary in my final implementation
void Polygon::getCameraDistance(Camera camera) {
  double totalDistance = 0;
  for (int i=0; i < (int) triangles.size(); i++) {
    totalDistance += triangles[i].getCameraDistance(camera);
  }
  cameraDistance = totalDistance / (int) triangles.size(); //average distance across all triangles
}

// constructor for Polygon derived class, populates vertices vector based on the input vector and 
// creates a fan of triangles based on that
Polygon::Polygon(std::vector<Eigen::Vector3d> inputVertices) {
  //populates vertices vector
    for (int i=0; i < (int) inputVertices.size(); i++) {
      vertices.push_back(inputVertices[i]);
    }
  //populates triangles vector
    Eigen::Vector3d baseVertex = vertices[0];
    for (int i=1; i < (int) vertices.size() - 1; i++) {
      Triangle newTriangle = Triangle(baseVertex, vertices[i], vertices[i+1], false);
      triangles.push_back(newTriangle);
    }
  }

// constructor for Polygon derived class, populates vertices vector based on the input vector and 
// creates a fan of triangles based on that
Polygon::Polygon(std::vector<Eigen::Vector3d> inputVertices, std::vector<Eigen::Vector3d> inputNormals) {
  //populates vertices vector
    for (int i=0; i < (int) inputVertices.size(); i++) {
      vertices.push_back(inputVertices[i]);
    }
  //populates triangles vector
    Eigen::Vector3d baseVertex = vertices[0];
    for (int i=1; i < (int) vertices.size() - 1; i++) {
      Triangle newTriangle = Triangle(baseVertex, vertices[i], vertices[i+1], true);
      newTriangle.vertex1Normal = inputNormals[0];
      newTriangle.vertex2Normal = inputNormals[i];
      newTriangle.vertex3Normal = inputNormals[i+1];
      triangles.push_back(newTriangle);
    }
  for (int i=0; i < (int) inputNormals.size(); i++) {
      normals.push_back(inputNormals[i]);
    }
  }


//triangle constructor, populates vertices based on the input arguments
Triangle::Triangle(Eigen::Vector3d firstVertex, Eigen::Vector3d secondVertex, Eigen::Vector3d thirdVertex, bool hasVertexNormals) {
  vertex1 = firstVertex;
  vertex2 = secondVertex;
  vertex3 = thirdVertex;
  surface_normal = ((secondVertex - firstVertex).cross(thirdVertex - firstVertex)).normalized();
  if (!hasVertexNormals) {
    vertex1Normal = surface_normal;
    vertex2Normal = surface_normal;
    vertex3Normal = surface_normal;
  }
}

//get light direction for vertex
Eigen::Vector3d Triangle::getLightDirection(Light light, Eigen::Vector3d vertex) {
  return (light.position - vertex).normalized();
}

//get view direction for vertex
Eigen::Vector3d Triangle::getViewDirection(Eigen::Vector3d from, Eigen::Vector3d vertex) {
  return (from - vertex).normalized(); //origin - intersection point
}


//overloaded operator (for debugging purposes)                                                                                                                               
bool operator < (const Polygon& polygon1, const Polygon& polygon2) {
  if (polygon1.cameraDistance < polygon2.cameraDistance)
    return true;
  return false;
}

//overloaded operator (for debugging purposes)                                                                                                                               
bool operator < (const Fragment& fragment1, const Fragment& fragment2) {
  if (fragment1.z < fragment2.z)
    return true;
  return false;
}

//overloaded operator (for debugging purposes)
ostream& operator<<(ostream& os, const Scene& scene) {
  for (int i=0; i < (int) scene.shapes.size(); i++) {
    Polygon currShape = (Polygon) scene.shapes[i];
    os << i << "\n";
    for (int i=0; i < (int) currShape.vertices.size(); i++) {
      os << currShape.vertices[i] << "\n";
    }
  }
  return os;
}

//overloaded operator (for debugging purposes)
ostream& operator<<(ostream& os, const Polygon& polygon) {
    for (int i=0; i < (int) polygon.vertices.size(); i++) {
      os << polygon.vertices[i] << "\n";
    }
  return os;
}

//overloaded operator (for debugging purposes)
ostream& operator<<(ostream& os, const Triangle& triangle) {
  os << triangle.vertex1 << "\n";
  os << triangle.vertex2 << "\n";
  os << triangle.vertex3 << "\n";
  return os;
}
