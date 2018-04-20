#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION, BARRIER} ConstraintType;   //seems redundant, but you can expand it
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint{
public:
  
  VectorXi globalIndices;         //|V_c| list of participating indices (out of global indices of the scene)
  double currValue;               //The current value of the constraint
  VectorXd currGradient;       //Current gradient of the constraint
  MatrixXd invMassMatrix;         //M^{-1} matrix size |V_c| X |V_c| only over the participating  vertices
  double refValue;                //Reference values to use in the constraint, when needed
  VectorXd refVector;             //Reference vector when needed
  double CRCoeff;
  ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
  ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality

  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const VectorXi& _globalIndices, const VectorXd& invMasses, const VectorXd& _refVector, const double _refValue, const double _CRCoeff):constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), refVector(_refVector), refValue(_refValue), CRCoeff(_CRCoeff){
    currValue=0.0;
	CRCoeff = _CRCoeff;
	refValue = _refValue;
	refVector = _refVector;
    globalIndices=_globalIndices;
    currGradient=VectorXd::Zero(globalIndices.size());
    invMassMatrix=invMasses.asDiagonal();
  }
  
  ~Constraint(){}
  
  
  //updating the value and the gradient vector with given values
  void updateValueGradient(const VectorXd& currVars){
    switch (constraintType){
        
      case DISTANCE: {
        /***************************
         TODO
         ***************************/

		  Vector3d p0(3);
		  p0 << currVars(0), currVars(1), currVars(2);
		  Vector3d p1(3);
		  p1 << currVars(3), currVars(4), currVars(5);
		  
		  currValue = (p0 - p1).norm()- refValue;
		  currGradient.segment(0, 3) = ((p0 - p1)) / (p0 - p1).norm();

		  currGradient.segment(3, 3) = -((p0 - p1)) / (p0 - p1).norm();

        break;
      }
        
      case COLLISION:{
        /***************************
         TODO
         ***************************/
		  /*
		  Vector3d p0(3);
		  Vector3d p1(3);
		  Vector3d p2(3);
		  Vector3d p3(3);


		  p0 << currVars(0), currVars(1), currVars(2);

		  p1 << currVars(3), currVars(4), currVars(5);
		  p2 << currVars(6), currVars(7), currVars(8);
		  p3 << currVars(9), currVars(10), currVars(11);

		  //VectorXd normal = ((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm());

		  double grv1 = ((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).transpose()*(p0-p1);

		  currGradient.segment(3, 3) = ((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm());

		  currGradient.segment(6, 3)  = ((p3 - p0).cross(p1 - p0) + (((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).cross(p3 - p0))*(((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).dot(p1 - p0))) / ((p2 - p0).cross(p3 - p0)).norm();

		  currGradient.segment(9, 3) = -((p2 - p0).cross(p1 - p0) + (((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).cross(p2 - p0))*(((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).dot(p1 - p0))) / ((p2 - p0).cross(p3 - p0)).norm();

		  currGradient.segment(0, 3) = -currGradient.segment(3, 3) - currGradient.segment(6, 3)- currGradient.segment(9, 3);

		  Vector3d tmp(3);
		  tmp = p0;

		  p0 << currVars(12), currVars(13), currVars(14);

		  p1 << currVars(15), currVars(16), currVars(17);
		  p2 << currVars(18), currVars(19), currVars(20);
		  p3 << currVars(21), currVars(22), currVars(23);


		  //VectorXd normal = ((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm());

		  double grv2 = ((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).transpose()*(p0 - p1);

		  currGradient.segment(15, 3) = ((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm());

		  currGradient.segment(18, 3) = ((p3 - p0).cross(p1 - p0) + (((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).cross(p3 - p0))*(((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).dot(p1 - p0))) / ((p2 - p0).cross(p3 - p0)).norm();

		  currGradient.segment(21, 3) = -((p2 - p0).cross(p1 - p0) + (((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).cross(p2 - p0))*(((p2 - p0).cross(p3 - p0) / (p2 - p0).cross(p3 - p0).norm()).dot(p1 - p0))) / ((p2 - p0).cross(p3 - p0)).norm();

		  currGradient.segment(12, 3) = -currGradient.segment(15, 3) - currGradient.segment(18, 3) - currGradient.segment(21, 3);


		  currValue = refVector.transpose()*(p0-tmp);
		  */

		  currValue = refVector.transpose()*(currVars);
		  currGradient = refVector;
		  //currGradient.segment(0, 12) = refVector;

		  //currGradient.segment(12, 12) = -refVector;
		  break;
      }
        
      case BARRIER:{
        /***************************
         TODO
         ***************************/

		currValue = currVars(0) - refValue;
		currGradient(0) = 1;
		break;
      }
    }
  }
  
  
  //computes the impulse needed for all particles to resolve the velocity constraint
  //returns true if constraint was already good
  bool resolveVelocityConstraint(const VectorXd& currPositions, const VectorXd& currVelocities, VectorXd& newImpulses, double tolerance){
    
    updateValueGradient(currPositions);

    if ((constraintEqualityType==INEQUALITY)&&((currValue>-tolerance)||((currGradient*currVelocities)(0,0)>-tolerance))){
      //constraint is valid, or velocity is already solving it, so nothing happens
      newImpulses=VectorXd::Zero(globalIndices.size());
      return true;
    }
    if ((constraintEqualityType==EQUALITY)&&(abs((currGradient*currVelocities)(0,0))<tolerance)){
      newImpulses=VectorXd::Zero(globalIndices.size());
      return true;
    }

    /***************************
     TODO
     ***************************/

	double lamda = -currValue / ((currGradient.transpose()*invMassMatrix*currGradient));

	newImpulses = lamda * invMassMatrix*currGradient;

    return false;
  }
  
  //projects the position unto the constraint
  //returns true if constraint was already good
  bool resolvePositionConstraint(const VectorXd& currPositions, const VectorXd& currVelocities, VectorXd& newPosDiffs, double tolerance){
    
    updateValueGradient(currPositions);

    //cout<<"C(currPositions): "<<currValue<<endl;
    //cout<<"currPositions: "<<currPositions<<endl;
    if ((constraintEqualityType==INEQUALITY)&&(currValue>-tolerance)){
      //constraint is valid
      newPosDiffs=VectorXd::Zero(globalIndices.size());
      return true;
    }
    
    if ((constraintEqualityType==EQUALITY)&&(abs(currValue)<tolerance)){
      newPosDiffs=VectorXd::Zero(globalIndices.size());
      return true;
    }
   
    /***************************
     TODO
     ***************************/

	double lamda = -currValue / ((currGradient.transpose()*invMassMatrix*currGradient));

	newPosDiffs = lamda*invMassMatrix*currGradient;

    return false;
  }
};



#endif /* constraints_h */
