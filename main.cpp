#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <imgui/imgui.h>
#include <iostream>
#include "scene.h"
#include <igl/copyleft/tetgen/tetrahedralize.h>

using namespace Eigen;
using namespace std;

MatrixXd V;
MatrixXi F;
MatrixXd P1, P2;   //for marking constraints
igl::opengl::glfw::Viewer mgpViewer[2];
int distConstV1, distConstV2;
Eigen::MatrixXi EConst;    //the user distance constraints in (V1,V2) format each row

MatrixXd platTriV;
MatrixXi platTriF;

float currTime = 0;
float timeStep = 0.02;
float CRCoeff= 1.0;
float friction = 0.0;

float power = 0.0;
float angleth = 0.0;
float anglephi = 360 * (3.14159265359 / 180.0); 

float airth = 0;
float airphi = 0;
float airforce = 0;

double tolerance=10e-6;
int maxIterations=100;
bool vertexChosen=false;
int currVertex=-1;
bool rcalcair = true;

Scene scene;

void CurrentBarDTWasModified() {
	// do whatever you need to do here to fire your signal
}


void createPlatform(MatrixXd& platV, MatrixXi& platF, RowVector3d& platCOM, RowVector4d& platOrientation)
{
  double platWidth=100.0;
  platCOM<<0.0,0.0,-0.0;
  platV.resize(8,3);
  platF.resize(12,3);
  platV<< -100, 0, -100,
	  -100, 10, -100,
	  100, 10, -100,
	  100, 0, -100,
	  -100, 0, 100,
	  -100, 10, 100,
	  100, 10, 100,
	  100, 0, 100;
  platF<< 0, 1, 3,
	  3, 1, 2,
	  0, 4, 1,
	  1, 4, 5,
	  3, 2, 7,
	  7, 2, 6,
	  4, 0, 3,
	  7, 4, 3,
	  6, 4, 7,
	  6, 5, 4,
	  1, 5, 6,
	  2, 1, 6;
  
  platOrientation<< 1.0, 0.0, 0.0, 1.0;
}


void update_mesh(igl::opengl::glfw::Viewer &viewer)
{
  viewer.data().clear();
  MatrixXi fullF(platTriF.rows()+F.rows(),3);
  fullF<<platTriF, F.array()+platTriV.rows();
  MatrixXd fullV(platTriV.rows()+V.rows(),3);
  fullV<<platTriV, V;
  viewer.data().set_mesh(fullV, fullF);

  Eigen::MatrixXd constV1(EConst.rows(),3), constV2(EConst.rows(),3);

  for (int i=0;i<EConst.rows();i++){
    constV1.row(i)=V.row(EConst(i,0));
    constV2.row(i)=V.row(EConst(i,1));
  }

  RowVector3d constColor; constColor << 0.4, 0.5, 0.3;

  //viewer.data().set_colors(constColor);


  viewer.data().set_face_based(true);
  viewer.data().add_edges(constV1, constV2, constColor);
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
	scene.shotmedown(power,angleth,anglephi, airth, airphi, airforce);
	rcalcair = true;
    return true;
  }

  if (key == 'P')
  {
	  viewer.core.is_animating = !viewer.core.is_animating;
	  return true;
  }

  if (key == 'R')
  {
	  currTime = 0;
	  scene.shenerestart();
	  return true;
  }

  if (key == 'S')
  {
    if (!viewer.core.is_animating){

      scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations,V,false, friction, currTime, power, angleth, anglephi, airth, airphi, airforce);
      currTime+=timeStep;
      update_mesh(viewer);
      std::cout <<"currTime: "<<currTime<<std::endl;
      return true;
    }
  }
  return false;
}


bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
  using namespace Eigen;
  using namespace std;
  
  if (viewer.core.is_animating){

	  if (rcalcair) {
		  airth = rand() % 361 * (3.14159265359 / 180.0);
		  airphi = rand() % 361 * (3.14159265359 / 180.0);
		  airforce = rand() % 5;
		  rcalcair = false;
	  }
    scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations, V,false, friction, currTime, power, angleth, anglephi, airth, airphi, airforce);

    update_mesh(viewer);
    currTime+=timeStep;
    //cout <<"currTime: "<<currTime<<endl;
  }
  
  return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{

  virtual void draw_viewer_menu() override
  {
    // Draw parent menu
    ImGuiMenu::draw_viewer_menu();
    
    // Add new group
    if (ImGui::CollapsingHeader("Game Options", ImGuiTreeNodeFlags_DefaultOpen))
    {  

	  ImGui::SliderFloat("Power", &power, 0, 450);
	  ImGui::SliderAngle("Angle thita", &angleth, 0, 90);
	  ImGui::SliderAngle("Angle phi", &anglephi, 270, 450);
	  
	  ImGui::Text("Press Space to Shot!");
    }
  }
};

bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifiers)
{
  if (!vertexChosen)
    return false;
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;

  Vector3f newPos=igl::unproject(Vector3f(x, y, viewer.down_mouse_z), (viewer.core.view * viewer.core.model).eval(), viewer.core.proj, viewer.core.viewport);
  

  if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
    //scene.createUserImpulse(currVertex, newPos);
    update_mesh(viewer);
    vertexChosen=false;
    return true;
  }
  
  if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Right){
    //creating new constraint
    int fid;
    Eigen::Vector3f bc;
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
                                viewer.core.proj, viewer.core.viewport, V, F, fid, bc));
    
    Eigen::MatrixXf::Index maxCol;
    bc.maxCoeff(&maxCol);
    int otherVertex=F(fid, maxCol);
    scene.addUserConstraint(currVertex, otherVertex, EConst);
    update_mesh(viewer);
    vertexChosen=false;
    return true;
  }
  
  return false;
}



int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
 
  
  //create platform

  RowVector3d platCOM;
  RowVector4d platOrientation;

  RowVector3d wallCOM;
  RowVector4d wallOrientation;

  createPlatform(platTriV, platTriF, platCOM, platOrientation);
  
  scene.loadScene(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]), F, EConst);

  scene.initScene(timeStep, 0.2, 0.2, V);
  scene.setPlatformBarriers(platTriV, CRCoeff);
  
  // Viewer Settings
  
  mgpViewer[0].callback_pre_draw = &pre_draw;
  mgpViewer[0].callback_key_down = &key_down;

  mgpViewer[0].core.animation_max_fps = 50.;
  
  CustomMenu menu;
  mgpViewer[0].plugins.push_back(&menu);

  update_mesh(mgpViewer[0]);

  mgpViewer[0].launch();
}
