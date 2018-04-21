#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <queue>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/edge_topology.h>
#include <igl/diag.h>
#include <igl/readMESH.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/igl_inline.h>
#include "constraints.h"
#include "auxfunctions.h"
#include "ccd.h"

using namespace Eigen;
using namespace std;


void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);

//the class the contains each individual rigid objects and their functionality
class Mesh{
public:
  
  //position
  VectorXd origPositions;     //3|V|x1 original vertex positions in xyzxyz format - never change this!
  VectorXd currPositions;     //3|V|x1 current vertex positions in xyzxyz format
  VectorXd currPositionsOld;    
  Vector3d ObjectCenter;
  int totaltime = 0; 
  int direction=-1;
  
  //kinematics
  bool isFixed;               //is the object immobile (infinite mass)
  VectorXd currImpulses;      //3|V|x1 correction impulses per coordinate
  VectorXd currVelocities;    //3|V|x1 velocities per coordinate in xyzxyz format.
  
  MatrixXi T;                 //|T|x4 tetrahdra
  VectorXd invMasses;         //|V|x1 inverse masses of vertices, computed in the beginning as 1.0/(density * vertex voronoi area)
  VectorXd voronoiVolumes;    //|V|x1 the voronoi volume of vertices
  VectorXd tetVolumes;        //|T|x1 tetrahedra volumes
  int globalOffset;           //the global index offset of the of opositions/velocities/impulses from the beginning of the global coordinates array in the containing scene class
  bool first_time = true;
  VectorXi boundTets;  //just the boundary tets, for collision
  
  double youngModulus, poissonRatio, density, alpha, beta;
  bool shot = false;

  SparseMatrix<double> A, K, M, D;   //The soft-body matrices
  
  SimplicialLLT<SparseMatrix<double>>* ASolver;   //the solver for the left-hand side matrix constructed for FEM
  
  ~Mesh(){if (ASolver!=NULL) delete ASolver;}
  
  //Quick-reject checking collision between mesh bounding boxes.
  bool isBoxCollide(const Mesh& m2){

    RowVector3d XMin1=RowVector3d::Constant(3276700.0);
    RowVector3d XMax1=RowVector3d::Constant(-3276700.0);
    RowVector3d XMin2=RowVector3d::Constant(3276700.0);
    RowVector3d XMax2=RowVector3d::Constant(-3276700.0);
    for (int i=0;i<origPositions.size();i+=3){
      XMin1=XMin1.array().min(currPositions.segment(i,3).array().transpose());
      XMax1=XMax1.array().max(currPositions.segment(i,3).array().transpose());
    }
    for (int i=0;i<m2.origPositions.size();i+=3){
      XMin2=XMin2.array().min(m2.currPositions.segment(i,3).array().transpose());
      XMax2=XMax2.array().max(m2.currPositions.segment(i,3).array().transpose());
    }
    
    /*double rmax1=vertexSphereRadii.maxCoeff();
     double rmax2=m2.vertexSphereRadii.maxCoeff();
     XMin1.array()-=rmax1;
     XMax1.array()+=rmax1;
     XMin2.array()-=rmax2;
     XMax2.array()+=rmax2;*/
    
    //checking all axes for non-intersection of the dimensional interval
    for (int i=0;i<3;i++)
      if ((XMax1(i)<XMin2(i))||(XMax2(i)<XMin1(i)))
        return false;
    
    return true;  //all dimensional intervals are overlapping = possible intersection
  }
  
  bool isNeighborTets(const RowVector4i& tet1, const RowVector4i& tet2){
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++)
        if (tet1(i)==tet2(j)) //shared vertex
          return true;
    
    return false;
  }
  
  
  //this function creates all collision constraints between vertices of the two meshes
  bool createCollisionConstraints(const Mesh& m, const bool sameMesh, const double timeStep, const double CRCoeff, vector<Constraint>& activeConstraints){
    

    //collision between bounding boxes
	if (!isBoxCollide(m)) {
		return false;
	}
    
	if ((isFixed && m.isFixed)) {//collision does nothing
		return false;
	}
    
    //creating tet spheres
    /*MatrixXd c1(T.rows(), 3);
     MatrixXd c2(m.T.rows(), 3);
     VectorXd r1(T.rows());
     VectorXd r2(m.T.rows());*/
    
    MatrixXd maxs1(boundTets.rows(),3);
    MatrixXd mins1(boundTets.rows(),3);
    MatrixXd maxs2(m.boundTets.rows(),3);
    MatrixXd mins2(m.boundTets.rows(),3);
    
    for (int i = 0; i < boundTets.size(); i++) {
      MatrixXd tet1(4, 3); tet1 << currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
      currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
      currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
      currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();
      
      //c1.row(i) = tet1.colwise().mean();
      //r1(i) = ((c1.row(i).replicate(4, 1) - tet1).rowwise().norm()).maxCoeff();
      mins1.row(i)=tet1.colwise().minCoeff();
      maxs1.row(i)=tet1.colwise().maxCoeff();
      
    }
    
    for (int i = 0; i < m.boundTets.size(); i++) {
      
      MatrixXd tet2(4, 3); tet2 << m.currPositions.segment(3 * m.T(m.boundTets(i), 0), 3).transpose(),
      m.currPositions.segment(3 * m.T(m.boundTets(i), 1), 3).transpose(),
      m.currPositions.segment(3 * m.T(m.boundTets(i), 2), 3).transpose(),
      m.currPositions.segment(3 * m.T(m.boundTets(i), 3), 3).transpose();
      
      //c2.row(i) = tet2.colwise().mean();
      //r2(i) = ((c2.row(i).replicate(4, 1) - tet2).rowwise().norm()).maxCoeff();
      mins2.row(i)=tet2.colwise().minCoeff();
      maxs2.row(i)=tet2.colwise().maxCoeff();
      
    }
    
    //checking collision between every tetrahedrons
    std::list<Constraint> collisionConstraints;
    for (int i=0;i<boundTets.size();i++){
      for (int j=0;j<m.boundTets.size();j++){
        
        //not checking for collisions between tetrahedra neighboring to the same vertices
        if (sameMesh)
          if (isNeighborTets(T.row(boundTets(i)), m.T.row(m.boundTets(j))))
            continue;  //not creating collisions between neighboring tets
        
        bool overlap=true;
        for (int k=0;k<3;k++)
          if ((maxs1(i,k)<mins2(j,k))||(maxs2(j,k)<mins1(i,k)))
            overlap=false;
        
        if (!overlap)
          continue;
        
        VectorXi globalCollisionIndices(24);
        VectorXd globalInvMasses(24);
        for (int t=0;t<4;t++){
          globalCollisionIndices.segment(3*t,3)<<globalOffset+3*(T(boundTets(i),t)), globalOffset+3*(T(boundTets(i),t))+1, globalOffset+3*(T(boundTets(i),t))+2;
          globalInvMasses.segment(3*t,3)<<invMasses(T(boundTets(i),t)), invMasses(T(boundTets(i),t)),invMasses(T(boundTets(i),t));
          globalCollisionIndices.segment(12+3*t,3)<<m.globalOffset+3*m.T(m.boundTets(j),t), m.globalOffset+3*m.T(m.boundTets(j),t)+1, m.globalOffset+3*m.T(m.boundTets(j),t)+2;
          globalInvMasses.segment(12+3*t,3)<<m.invMasses(m.T(m.boundTets(j),t)), m.invMasses(m.T(m.boundTets(j),t)),m.invMasses(m.T(m.boundTets(j),t));
        }
        
        ccd_t ccd;
        CCD_INIT(&ccd);
        ccd.support1       = support; // support function for first object
        ccd.support2       = support; // support function for second object
        ccd.center1         =center;
        ccd.center2         =center;
        
        ccd.first_dir       = stub_dir;
        ccd.max_iterations = 100;     // maximal number of iterations
        
        MatrixXd tet1(4, 3); tet1 << currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();
        
        MatrixXd tet2(4, 3); tet2 << m.currPositions.segment(3 * m.T(m.boundTets(j), 0), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(j), 1), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(j), 2), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(j), 3), 3).transpose();
        
        void* obj1=(void*)&tet1;
        void* obj2=(void*)&tet2;
        
        ccd_real_t _depth;
        ccd_vec3_t dir, pos;
        
        int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);
        
        if (nonintersect)
          continue;
        
        Vector3d intNormal, intPosition;
        double depth;
        for (int k=0;k<3;k++){
          intNormal(k)=dir.v[k];
          intPosition(k)=pos.v[k];
        }
        
        depth =_depth;
        intPosition-=depth*intNormal/2.0;
        
        Vector3d p1=intPosition+depth*intNormal;
        Vector3d p2=intPosition;
        
        //getting barycentric coordinates of each point
        
        MatrixXd PMat1(4,4); PMat1<<1.0,currPositions.segment(3*T(boundTets(i),0),3).transpose(),
        1.0,currPositions.segment(3*T(boundTets(i),1),3).transpose(),
        1.0,currPositions.segment(3*T(boundTets(i),2),3).transpose(),
        1.0,currPositions.segment(3*T(boundTets(i),3),3).transpose();
        PMat1.transposeInPlace();
        
        Vector4d rhs1; rhs1<<1.0,p1;
        
        Vector4d B1=PMat1.inverse()*rhs1;
        
        MatrixXd PMat2(4,4); PMat2<<1.0,m.currPositions.segment(3*m.T(m.boundTets(j),0),3).transpose(),
        1.0,m.currPositions.segment(3*m.T(m.boundTets(j),1),3).transpose(),
        1.0,m.currPositions.segment(3*m.T(m.boundTets(j),2),3).transpose(),
        1.0,m.currPositions.segment(3*m.T(m.boundTets(j),3),3).transpose();
        PMat2.transposeInPlace();
        
        Vector4d rhs2; rhs2<<1.0,p2;
        
        Vector4d B2=PMat2.inverse()*rhs2;
        
        //cout<<"B1: "<<B1<<endl;
        //cout<<"B2: "<<B2<<endl;
        
        //Matrix that encodes the vector between interpenetration points by the c
        MatrixXd v2cMat1(3,12); v2cMat1.setZero();
        for (int k=0;k<3;k++){
          v2cMat1(k,k)=B1(0);
          v2cMat1(k,3+k)=B1(1);
          v2cMat1(k,6+k)=B1(2);
          v2cMat1(k,9+k)=B1(3);
        }
        
        MatrixXd v2cMat2(3,12); v2cMat2.setZero();
        for (int k=0;k<3;k++){
          v2cMat2(k,k)=B2(0);
          v2cMat2(k,3+k)=B2(1);
          v2cMat2(k,6+k)=B2(2);
          v2cMat2(k,9+k)=B2(3);
        }
        
        MatrixXd v2dMat(3,24); v2dMat<<-v2cMat1, v2cMat2;
        VectorXd constVector=intNormal.transpose()*v2dMat;
        
        //cout<<"intNormal: "<<intNormal<<endl;
        //cout<<"n*(p2-p1): "<<intNormal.dot(p2-p1)<<endl;
        collisionConstraints.push_back(Constraint(COLLISION, INEQUALITY,  globalCollisionIndices, globalInvMasses, constVector, 0, CRCoeff));
        
        //i=10000000;
        //break;
        
      }
    }

    activeConstraints.insert(activeConstraints.end(), collisionConstraints.begin(), collisionConstraints.end());
	return true;
  }
  
  //where the matrices A,M,K,D are created and factorized to ASolver at each change of time step, or beginning of time.
  void createGlobalMatrices(const double timeStep, const double _alpha, const double _beta)
  {
    
    /***************************
     
     ***************************/

	  A = SparseMatrix<double>(currVelocities.size() , currVelocities.size() );
	  M = SparseMatrix<double>(currVelocities.size() , currVelocities.size() );
	  K = SparseMatrix<double>(currVelocities.size() , currVelocities.size() );
	  D = SparseMatrix<double>(currVelocities.size() , currVelocities.size() );
	  A.setZero();
	  M.setZero();
	  K.setZero();
	  D.setZero();
	  //youngModulus, poissonRatio, density
	  double m = youngModulus / (2 * (1 + poissonRatio));
	  double l = poissonRatio * youngModulus / ((1 + poissonRatio)*(1 - 2 * poissonRatio));

	  MatrixXd Pe(4, 4);
	  MatrixXd Ge(3, 4);
	  MatrixXd ZeroI(3, 4);
	  MatrixXd Je(9, 12);
	  MatrixXd Ddeform(6, 9);
	  MatrixXd Be(6, 12);
	  MatrixXd Ce(6, 6);
	  MatrixXd Ke(12, 12);
	  SparseMatrix<double> Kt(tetVolumes.size()*12 , tetVolumes.size()*12);
	  Kt.setZero();
	  Ce.setZero();

	  MatrixXd Output(3, 3);
	  vector<Triplet<double>> triplets;


	  SparseMatrix<double> Q(tetVolumes.size() * 12, currVelocities.size());
	  //Build the Q table
	  Q.setZero();

	  for (int t = 0; t < T.rows(); t++) {
		  for (int i = 0; i < 3; i++) {
			  for (int j = 0; j < 4; j++) {
				  triplets.push_back((Triplet<double>(t*12+i*4+j,T(t,j)*3+i, 1)));
			  }
		  }
	  }

	  Q.setFromTriplets(triplets.begin(), triplets.end());

	  triplets.clear();

	  //Build the ZeroI table
	  for (int i = 0; i < 3; i++) {
		  for (int j = 0; j < 4; j++) {
			  if (j == 0)
				  ZeroI(i, j) = 0;
			  else
				  if (i==j-1)
					  ZeroI(i, j) = 1;
				  else 
					  ZeroI(i, j) = 0;
		  }
	  }

	  //Build the Ce table
	  for (int i = 0; i < 3; i++) {
		  for (int j = 0; j < 3; j++) {
			  if (i == j) {
				  Ce(i, j) = l + 2 * m;
				  Ce(i + 3, j + 3) = l;
			  }
			  else {
				  Ce(i, j) = l;
			  }
		  }
	  }

	  //Build the Ddeform table
	  Ddeform.setZero();
	  Ddeform(0, 0) = 1;
	  Ddeform(1, 4) = 1;
	  Ddeform(2, 8) = 1;
	  Ddeform(3, 1) = 0.5;
	  Ddeform(3, 3) = 0.5;
	  Ddeform(4, 5) = 0.5;
	  Ddeform(4, 7) = 0.5;
	  Ddeform(5, 2) = 0.5;
	  Ddeform(5, 6) = 0.5;


	  for (int t = 0; t < T.rows(); t++) {

		  //Build the Pe table
		  for (int i = 0; i < 4; i++) {
			  for (int j = 0; j < 4; j++) {
				  if (j == 0)
					  Pe(i, j) = 1;
				  else
					  Pe(i, j) = origPositions.segment(3 * T(t, i), 3)[j-1];
					  
			  }
		  }

		  //Build the Ge table
		  Ge = ZeroI * Pe.inverse();

		  //Build the Je table
		  Je.setZero();
		  for (int k = 0; k < 3; k++) {
			  for (int i = 0; i < 3; i++) {
				  for (int j = 0; j < 4; j++) {
					  Je(i + k * 3, j + k * 4) = Ge(i, j);
				  }
			  }
		  }


		  //Build the Be table
		  Be = Ddeform * Je;

		  //Build the Ke table
		  Ke = 0.5 * tetVolumes(t) * Be.transpose() * Ce*Be;

		  for (int i = 0; i < 12; i++) {
			  for (int j = 0; j < 12; j++) {
				  triplets.push_back((Triplet<double>(i + t * 12, j + t * 12, Ke(i,j))));
			  }
		  }

	  }

	  Kt.setFromTriplets(triplets.begin(), triplets.end());

	  K = Q.transpose()*Kt*Q;


	  triplets.clear();

	  double Output2 = 0;

	  for (int i = 0; i < voronoiVolumes.size(); i++) {
		  Output2 = density* voronoiVolumes(i);
		  for (int row = 0; row < 3; row++) {
			  for (int col = 0; col < 3; col++) {
				  if (row==col)
					triplets.push_back((Triplet<double>(col + i * 3, row + i * 3, Output2)));
			  }
		  }
	  }

	  M.setFromTriplets(triplets.begin(), triplets.end());

	  D = _alpha * M + _beta * K;

	  A = M + timeStep * D + timeStep * timeStep*K;

    if (ASolver==NULL)
      ASolver=new SimplicialLLT<SparseMatrix<double>>();
    
	ASolver->compute(A);
	
  }

  
  //computes tet volumes, masses, and allocate voronoi areas and inverse masses to
  Vector3d initializeVolumesAndMasses()
  {

    tetVolumes.conservativeResize(T.rows());
    voronoiVolumes.conservativeResize(origPositions.size()/3);
    voronoiVolumes.setZero();
    invMasses.conservativeResize(origPositions.size()/3);
    Vector3d COM; COM.setZero();
    for (int i=0;i<T.rows();i++){
      Vector3d e01=origPositions.segment(3*T(i,1),3)-origPositions.segment(3*T(i,0),3);
      Vector3d e02=origPositions.segment(3*T(i,2),3)-origPositions.segment(3*T(i,0),3);
      Vector3d e03=origPositions.segment(3*T(i,3),3)-origPositions.segment(3*T(i,0),3);
      Vector3d tetCentroid=(origPositions.segment(3*T(i,0),3)+origPositions.segment(3*T(i,1),3)+origPositions.segment(3*T(i,2),3)+origPositions.segment(3*T(i,3),3))/4.0;
      tetVolumes(i)=std::abs(e01.dot(e02.cross(e03)))/6.0;
      for (int j=0;j<4;j++)
        voronoiVolumes(T(i,j))+=tetVolumes(i)/4.0;
      
      COM+=tetVolumes(i)*tetCentroid;
    }
    
    COM.array()/=tetVolumes.sum();

	

    for (int i=0;i<origPositions.size()/3;i++)
      invMasses(i)=1.0/(voronoiVolumes(i)*density);
	ObjectCenter = COM;
    return COM;
    
  }
  
  //performing the integration step of the soft body.
  void integrateVelocity(double timeStep,int index, double power, double angleth, double anglephi, double airth, double airphi, double airforce){
    
    /***************************
     
     ***************************/
	if (isFixed)
		return;

	double grav = 9.81;

	for (int i = 0; i < currVelocities.rows() / 3; i++)
	{
		if (shot) {
			currVelocities(i * 3 + 1) -= grav * timeStep;
		}
		if (index != 3 && index != 4) {

			currVelocities(i * 3 ) = currVelocities(i * 3 ) + direction * grav * timeStep;

			if (currVelocities(i * 3 ) < -10)
				direction = direction * -1;
			else if (currVelocities(i * 3 )>10)
				direction = direction * -1;
		}
	}

	if (index == 4) {
		if (!shot) {
			if (totaltime == 50) {
				for (int i = 0; i < currVelocities.rows() / 3; i++)
				{
					currVelocities(i * 3 + 2) = -sqrt(2 * power / 0.1)*cos(angleth)*cos(anglephi); //z
					currVelocities(i * 3 + 1) = sqrt(2 * power / 0.1)*sin(angleth)*cos(anglephi); //y
					currVelocities(i * 3) = sqrt(2 * power / 0.1)*sin(anglephi); //x

					currVelocities(i * 3 + 2) += -sqrt(2 * airforce / 0.1)*cos(airth)*cos(airphi); //z
					currVelocities(i * 3 + 1) += sqrt(2 * airforce / 0.1)*sin(airth)*cos(airphi); //y
					currVelocities(i * 3) += sqrt(2 * airforce / 0.1)*sin(airphi); //x
				}
				currPositions = origPositions;
				totaltime = 0;
			}
			else {
				for (int i = 0; i < currVelocities.rows() / 3; i++)
				{
					currVelocities(i * 3 + 1) -= grav * timeStep;
				}
				totaltime++;
			}
		}
		else 
			currPositions = origPositions;

	}

  }


  void setshotvel(double power,double angleth,double anglephi, double airth, double airphi, double airforce) {

	  shot = true;
	  double grav = 9.81;


	  for (int i = 0; i < currVelocities.rows() / 3; i++)
	  {
		  currVelocities(i * 3 +2) =    -sqrt(2 * power / 0.1)*cos(angleth)*cos(anglephi); //z
		  currVelocities(i * 3 + 1) = sqrt(2 * power / 0.1)*sin(angleth)*cos(anglephi); //y
		  currVelocities(i * 3 ) = sqrt(2 * power / 0.1)*sin(anglephi); //x

		  currVelocities(i * 3 + 2) += -sqrt(2 * airforce / 0.1)*cos(airth)*cos(airphi); //z
		  currVelocities(i * 3 + 1) += sqrt(2 * airforce / 0.1)*sin(airth)*cos(airphi); //y
		  currVelocities(i * 3) += sqrt(2 * airforce / 0.1)*sin(airphi); //x

		  currVelocities(i * 3 + 1) -= grav * 0.02;
	  }
  }


  void restartshot() {
	  shot = false;
	  currPositions = origPositions;
	  currVelocities = VectorXd::Zero(origPositions.rows());
	  currImpulses = VectorXd::Zero(origPositions.rows());

  }

  //performing the integration step of the soft body.
  void integrateVelocities(double timeStep,double currTime) {

	  /***************************

	  ***************************/
	  if (isFixed)
		  return;

	  VectorXd b= VectorXd::Zero(currVelocities.rows());;
	  SparseMatrix<double> psi(currVelocities.size(), currVelocities.size());
	  psi.setZero();

	  VectorXd Fgrav = VectorXd::Zero(currVelocities.rows());


	  vector<Triplet<double>> triplets;

	  double grav = -9.81;
	  for (int i = 0; i < currVelocities.size()/3; i++)
	  {
		  Fgrav(i * 3 + 1) = grav * currTime;
	  }
	  Fgrav = M * Fgrav;
	  if (shot)
		b = M * currVelocities - timeStep*((K*(currPositions - origPositions)));
	  else
		  b = M * currVelocities - timeStep * ((K*(currPositions - origPositions)));

	  currVelocities = ASolver->solve(b);

	  integratePosition(timeStep,false, currTime);
  }

  //Update the current position with the integrated velocity
  void integratePosition(double timeStep, bool squeezing,double currTime){
    if (isFixed)
      return;  //a fixed object is immobile
    
    currPositions+= currVelocities * timeStep;


	if (squeezing)
		currPositions+=VectorXd::Random(currPositions.size());

  }
  
  //the full integration for the time step (velocity + position)
  void integrate(double timeStep,bool squeezing,double currTime,int i, double power, double angleth, double anglephi, double airth, double airphi, double airforce){
	  
    integrateVelocity(timeStep,i, power, angleth, anglephi, airth, airphi, airforce);
    integratePosition(timeStep, squeezing, currTime);
  }
  
  
  Mesh(const VectorXd& _origPositions, const MatrixXi& boundF, const MatrixXi& _T, const int _globalOffset, const double _youngModulus, const double _poissonRatio, const double _density, const bool _isFixed, const RowVector3d& userCOM, const RowVector4d& userOrientation){
    origPositions=_origPositions;
    //cout<<"original origPositions: "<<origPositions<<endl;
    T=_T;
    isFixed=_isFixed;
    globalOffset=_globalOffset;
    density=_density;
    poissonRatio=_poissonRatio;
    youngModulus=_youngModulus;
    currVelocities=VectorXd::Zero(origPositions.rows());
    currImpulses=VectorXd::Zero(origPositions.rows());
    
    VectorXd naturalCOM=initializeVolumesAndMasses();
    //cout<<"naturalCOM: "<<naturalCOM<<endl;


    origPositions-= naturalCOM.replicate(origPositions.rows()/3,1);  //removing the natural COM of the OFF file (natural COM is never used again)
    //cout<<"after natrualCOM origPositions: "<<origPositions<<endl;
    
    for (int i=0;i<origPositions.size();i+=3)
      origPositions.segment(i,3)<<(QRot(origPositions.segment(i,3).transpose(), userOrientation)+userCOM).transpose();
    
    currPositions=origPositions; 
	currPositionsOld = origPositions;
    
    if (isFixed)
      invMasses.setZero();
    
    //finding boundary tets
    VectorXi boundVMask(origPositions.rows()/3);
    boundVMask.setZero();
    for (int i=0;i<boundF.rows();i++)
      for (int j=0;j<3;j++)
        boundVMask(boundF(i,j))=1;
    
    cout<<"boundVMask.sum(): "<<boundVMask.sum()<<endl;
    
    vector<int> boundTList;
    for (int i=0;i<T.rows();i++){
      int incidence=0;
      for (int j=0;j<4;j++)
        incidence+=boundVMask(T(i,j));
      if (incidence>2)
        boundTList.push_back(i);
    }
    
    boundTets.resize(boundTList.size());
    for (int i=0;i<boundTets.size();i++)
      boundTets(i)=boundTList[i];
    
    ASolver=NULL;
  }
  
};


//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime;
  
  VectorXd globalPositions;   //3*|V| all positions
  VectorXd globalVelocities;  //3*|V| all velocities
  VectorXd globalImpulses;
  VectorXd globalInvMasses;   //3*|V| all inverse masses  (NOTE: the invMasses in the Mesh class is |v| (one per vertex)!
  MatrixXi globalT;           //|T|x4 tetraheda in global index
  
  double oldtimeStep;
  vector<Mesh> meshes;
  
  vector<Constraint> userConstraints;   //provided from the scene
  vector<Constraint> barrierConstraints;  //provided by the platform
  
  //updates from global values back into mesh values
  void global2Mesh(){
    for (int i=0;i<meshes.size();i++){
      meshes[i].currPositions<<globalPositions.segment(meshes[i].globalOffset, meshes[i].currPositions.size());
      meshes[i].currVelocities<<globalVelocities.segment(meshes[i].globalOffset, meshes[i].currVelocities.size());
	  meshes[i].currImpulses << globalImpulses.segment(meshes[i].globalOffset, meshes[i].currImpulses.size());
    }
  }
  
  //update from mesh current values into global values
  void mesh2global(){
    for (int i=0;i<meshes.size();i++){
      globalPositions.segment(meshes[i].globalOffset, meshes[i].currPositions.size())<<meshes[i].currPositions;
      globalVelocities.segment(meshes[i].globalOffset, meshes[i].currVelocities.size())<< meshes[i].currVelocities;
	  globalImpulses.segment(meshes[i].globalOffset, meshes[i].currImpulses.size()) << meshes[i].currImpulses;
    }
  }
  
  
  //This should be called whenever the timestep changes
  void initScene(double timeStep, const double alpha, const double beta, MatrixXd& viewerV){
	
    for (int i=0;i<meshes.size();i++){
      if (!meshes[i].isFixed)
        meshes[i].createGlobalMatrices(timeStep, alpha, beta);
    }
	oldtimeStep = timeStep;
    mesh2global();
    
    //cout<<"globalPositions: "<<globalPositions<<endl;
    //updating viewer vertices
    viewerV.conservativeResize(globalPositions.size()/3,3);
    for (int i=0;i<globalPositions.size();i+=3)
      viewerV.row(i/3)<<globalPositions.segment(i,3).transpose();
  }

  void shotmedown(double power, double angleth,double anglephi,double airth, double airphi, double airforce) {
	  meshes[3].isFixed = false;
	  meshes[4].shot = true;
	  meshes[3].setshotvel(power, angleth, anglephi, airth, airphi, airforce);

  }

  void shenerestart() {

	  meshes[3].restartshot();

	  meshes[4].shot=false;
  }
  /*********************************************************************
   This function handles a single time step
   1. Integrating velocities and position from forces and previous impulses
   2. detecting collisions and generating collision constraints, alongside with given user constraints
   3. Resolving constraints iteratively by updating velocities until the system is valid (or maxIterations has passed)
   *********************************************************************/
  
  void updateScene(double timeStep, double CRCoeff, const double tolerance, const int maxIterations, MatrixXd& viewerV,bool squeezing,double friction,double currTime, double power, double angleth, double anglephi, double airth, double airphi, double airforce){
	 
    /*******************1. Integrating velocity and position from external and internal forces************************************/
	

	for (int i = 0; i < meshes.size(); i++)
		meshes[i].integrate(timeStep, squeezing, currTime, i, power, angleth, anglephi, airth, airphi, airforce);
		
    
    mesh2global();
    

    /*******************2. Creating and Aggregating constraints************************************/
    
    vector<Constraint> activeConstraints;

    //user constraints
    activeConstraints.insert(activeConstraints.end(), userConstraints.begin(), userConstraints.end());

    //barrier constraints
    activeConstraints.insert(activeConstraints.end(), barrierConstraints.begin(), barrierConstraints.end());


    //collision constraints
	int meshsize = meshes.size();
	for (int i = 0; i < meshsize; i++) {
		if (i == 4) continue;
		for (int j = i + 1; j < meshsize; j++)
		{
			if (j == 4) continue;
			if (meshes[i].createCollisionConstraints(meshes[j], i == j, timeStep, CRCoeff, activeConstraints)) {
				if (j != 4) {
					meshes[j].isFixed = true;

					//Mesh m = meshes[i];
					meshes[i].isFixed=true;

				}
				if (j == 9) {
					//meshes[i].isFixed = true;
				}
			}
		}
	}

	mesh2global();
    /*******************3. Resolving velocity constraints iteratively until the velocities are valid************************************/
    
    /***************************
     TODO
     ***************************/
    
	 //----------------------------Full Constraints Create-----------------------------------------------------------


	vector<VectorXd> ForPositions(activeConstraints.size());
	vector<VectorXd> ForImpules(activeConstraints.size());
	//activeConstraints.size()
	for (int i = 0; i < activeConstraints.size(); i++)
	{
		VectorXi particleIndices = activeConstraints.at(i).globalIndices;

		VectorXd currPosR(particleIndices.size());
		VectorXd currVelR(particleIndices.size());

		for (int j = 0; j < particleIndices.size(); j++)
		{
			currPosR(j) = globalPositions(activeConstraints[i].globalIndices(j));
			currVelR(j) = globalVelocities(activeConstraints[i].globalIndices(j));
		}

		activeConstraints.at(i).resolvePositionConstraint(currPosR, currVelR, ForPositions[i], 0);


		for (int j = 0; j < particleIndices.size(); j++)
		{
			globalPositions(activeConstraints[i].globalIndices(j)) += ForPositions[i](j);

		}

		if ((activeConstraints.at(i).constraintType != DISTANCE) && timeStep != 0)
		{
			activeConstraints.at(i).resolveVelocityConstraint(currPosR, currVelR, ForImpules[i], 0);
		}

		for (int j = 0; j < particleIndices.size(); j++)
		{

			if ((activeConstraints.at(i).constraintType != DISTANCE)&&(activeConstraints.at(i).constraintType != COLLISION) && timeStep != 0)
			{
				globalImpulses(activeConstraints[i].globalIndices(j)) += ((1 + CRCoeff)*ForImpules[i](j)/ timeStep) ;
			}
		}
	}

	for (int i = 0; i < globalVelocities.rows(); i++)
	{
		globalVelocities.row(i) += globalImpulses.row(i) ;
	}

	globalImpulses.setZero();

    global2Mesh();
    
    /*******************4. Solving for position drift************************************/
    
    mesh2global();

    /***************************
     TODO
     ***************************/

	//collision constraints
	for (int i = 0; i < meshes.size(); i++){
		meshes[i].integrateVelocities(timeStep,currTime);
	}

	if (oldtimeStep != timeStep) {
		for (int i = 0; i<meshes.size(); i++) {
			if (!meshes[i].isFixed)
				meshes[i].createGlobalMatrices(timeStep, 0.02, 0.02);
		}
		oldtimeStep = timeStep;
	}

    //global2Mesh();
    
    //updating viewer vertices
    viewerV.conservativeResize(globalPositions.size()/3,3);


    for (int i=0;i<globalPositions.size();i+=3)
      viewerV.row(i/3)<<globalPositions.segment(i,3).transpose();
  }
  
  //adding a constraint from the user
  void addUserConstraint(const int currVertex, const int otherVertex, MatrixXi& viewerEConst)
  {
    
    VectorXi coordIndices(6);
    coordIndices<<3*currVertex,3*currVertex+1,3*currVertex+2,3*otherVertex,3*otherVertex+1,3*otherVertex+2;
    
    VectorXd constraintInvMasses(6);
    constraintInvMasses<<globalInvMasses(currVertex),globalInvMasses(currVertex),globalInvMasses(currVertex),
    globalInvMasses(otherVertex),globalInvMasses(otherVertex),globalInvMasses(otherVertex);
    double refValue=(globalPositions.segment(3*currVertex,3)-globalPositions.segment(3*otherVertex,3)).norm();

    userConstraints.push_back(Constraint(DISTANCE, EQUALITY, coordIndices, constraintInvMasses, MatrixXd::Zero(1,1), refValue,0.0));
    
    viewerEConst.conservativeResize(viewerEConst.rows()+1,2);
    viewerEConst.row(viewerEConst.rows()-1)<<currVertex, otherVertex;
    
  }
  
  void setPlatformBarriers(const MatrixXd& platV, const double CRCoeff){
    
    RowVector3d minPlatform=platV.colwise().minCoeff();
    RowVector3d maxPlatform=platV.colwise().maxCoeff();

    //y value of maxPlatform is lower bound
    for (int i=1;i<globalPositions.size();i+=3){
      VectorXi coordIndices(1); coordIndices(0)=i;
      VectorXd constraintInvMasses(1); constraintInvMasses(0)=globalInvMasses(i);
      barrierConstraints.push_back(Constraint(BARRIER, INEQUALITY, coordIndices, constraintInvMasses, MatrixXd::Zero(1,1), maxPlatform(1),CRCoeff));

	}
    
  }
  
  
  //adding an object.
  void addMesh(const MatrixXd& V, const MatrixXi& boundF, const MatrixXi& T, const double youngModulus, const double PoissonRatio,  const double density, const bool isFixed, const RowVector3d& userCOM, const RowVector4d userOrientation){
    
    VectorXd Vxyz(3*V.rows());
    for (int i=0;i<V.rows();i++)
      Vxyz.segment(3*i,3)=V.row(i).transpose();
    
    //cout<<"Vxyz: "<<Vxyz<<endl;
    Mesh m(Vxyz,boundF, T, globalPositions.size(), youngModulus, PoissonRatio, density, isFixed, userCOM, userOrientation);
	
    meshes.push_back(m);
    int oldTsize=globalT.rows();
    globalT.conservativeResize(globalT.rows()+T.rows(),4);
    globalT.block(oldTsize,0,T.rows(),4)=T.array()+globalPositions.size()/3;  //to offset T to global index
    globalPositions.conservativeResize(globalPositions.size()+Vxyz.size());
    globalVelocities.conservativeResize(globalPositions.size());
	globalImpulses.conservativeResize(globalPositions.size());
    int oldIMsize=globalInvMasses.size();
    globalInvMasses.conservativeResize(globalPositions.size());
    for (int i=0;i<m.invMasses.size();i++)
      globalInvMasses.segment(oldIMsize+3*i,3)=Vector3d::Constant(m.invMasses(i));
    
    mesh2global();
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string dataFolder, const std::string sceneFileName, const std::string constraintFileName, MatrixXi& viewerF, MatrixXi& viewerEConst){
    
    ifstream sceneFileHandle;
    ifstream constraintFileHandle;
    sceneFileHandle.open(dataFolder+std::string("/")+sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;

    int numofObjects, numofConstraints;
    
    currTime=0;
    sceneFileHandle>>numofObjects;
    for (int i=0;i<numofObjects;i++){
      MatrixXi objT, objF;
      MatrixXd objV;
      std::string MESHFileName;
      bool isFixed;
      double youngModulus, poissonRatio, density;
      RowVector3d userCOM;
      RowVector4d userOrientation;
      sceneFileHandle>>MESHFileName>>density>>youngModulus>>poissonRatio>>isFixed;
      sceneFileHandle>>userCOM(0)>>userCOM(1)>>userCOM(2)>>userOrientation(0)>>userOrientation(1)>>userOrientation(2)>>userOrientation(3);
      userOrientation.normalize();
      //if the mesh is an OFF file, tetrahedralize it
      if (MESHFileName.find(".off") != std::string::npos){
        MatrixXd VOFF;
        MatrixXi FOFF;
        igl::readOFF(dataFolder+std::string("/")+MESHFileName,VOFF,FOFF);
        if (!isFixed)
          igl::copyleft::tetgen::tetrahedralize(VOFF,FOFF,"pq1.1", objV,objT,objF);
        else
          igl::copyleft::tetgen::tetrahedralize(VOFF,FOFF,"pq1.414Y", objV,objT,objF);
      } else {
        igl::readMESH(dataFolder+std::string("/")+MESHFileName,objV,objT, objF);
      }
      
      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;
      
      int oldFSize=viewerF.rows();
      viewerF.conservativeResize(viewerF.rows()+objF.rows(),3);
      viewerF.block(oldFSize,0,objF.rows(),3)=objF.array()+globalPositions.size()/3;
      //cout<<"objF: "<<objF<<endl;
      //cout<<"viewerF: "<<viewerF<<endl;

      addMesh(objV,objF, objT, youngModulus, poissonRatio,  density, isFixed, userCOM, userOrientation);
    }
    
    return true;
  }
  
  
  Scene(){}
  ~Scene(){}
};



/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p)
{
  // assume that obj_t is user-defined structure that holds info about
  // object (in this case box: x, y, z, pos, quat - dimensions of box,
  // position and rotation)
  //std::cout<<"calling support"<<std::endl;
  MatrixXd *obj = (MatrixXd *)_obj;
  RowVector3d p;
  RowVector3d d;
  for (int i=0;i<3;i++)
    d(i)=_d->v[i]; //p(i)=_p->v[i];
  
  
  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;
  
  RowVector3d objCOM=obj->colwise().mean();
  int maxVertex=-1;
  int maxDotProd=-32767.0;
  for (int i=0;i<obj->rows();i++){
    double currDotProd=d.dot(obj->row(i)-objCOM);
    if (maxDotProd < currDotProd){
      maxDotProd=currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex=i;
    }
    
  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;
  
  for (int i=0;i<3;i++)
    _p->v[i]=(*obj)(maxVertex,i);
  
  //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir)
{
  dir->v[0]=1.0;
  dir->v[1]=0.0;
  dir->v[2]=0.0;
}

void center(const void *_obj,ccd_vec3_t *center)
{
  MatrixXd *obj = (MatrixXd *)_obj;
  RowVector3d objCOM=obj->colwise().mean();
  for (int i=0;i<3;i++)
    center->v[i]=objCOM(i);
}




#endif
