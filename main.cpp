#include <iostream>
#include <fstream>
#include "Scene.h"

using namespace std;
int main(int argc, char *argv[]) {
  Scene newScene(argv[1]);
  
  std::vector<Polygon> shapes;
  for (int i=0; i < (int) (newScene.shapes).size(); i++) {
    shapes.push_back(newScene.shapes[i]);
    //cout << shapes[i] << endl;
  }
  
  for (int i=0; i < newScene.camera.resy; i++) {
    for (int j=0; j < newScene.camera.resx; j++) {
      //cout << newScene.pixels[i][j][0] << " " << newScene.pixels[i][j][1] << newScene.pixels[i][j][2] << endl;
    }
  }
  newScene.rasterize();
  newScene.WriteToPPM();
  return 0;
}