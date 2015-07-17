#include "MyWorld.h"
#include "dart/utils/Paths.h"
#include "dart/utils/SkelParser.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/Joint.h"
#include "dart/dynamics/Marker.h"
#include <string>
#include <iostream>

using namespace Eigen;
using namespace dart::dynamics;

using namespace std;

string nodeName[15];
MyWorld::MyWorld() {
    // Load a skeleton from file
    mSkel = dart::utils::SkelParser::readSkeleton(DART_DATA_PATH"skel/human.skel");
    // Assume that there is only one constraint
    mJ = MatrixXd::Zero(15 * 3, mSkel->getNumDofs());
    mC = VectorXd::Zero(15 * 3);
	
	
	

    mConstrainedMarker = -1;

    mTgtB = VectorXd::Zero(45);

    for (int i = 0; i < 15; i++)
    {
        aMrk[i] = false;
		BodyNode* node = mSkel->getMarker(i)->getBodyNode();
		nodeName[i] = node->getName();
    }


}

MyWorld::~MyWorld() {
    delete mSkel;
}

void MyWorld::solve() {
    if (mConstrainedMarker == -1)
        return; 
    int numIter = 300;
    double alpha = 0.01;
    int nDof = mSkel->getNumDofs();
    VectorXd gradients(nDof);
    VectorXd newPose(nDof);
    for (int i = 0; i < numIter; i++) {
        gradients = updateGradients();
        newPose = mSkel->getPositions() - alpha * gradients;
        mSkel->setPositions(newPose); 
        mSkel->computeForwardKinematics(true, false, false); // DART updates all the transformations based on newPose
    }
}

// Current code only works for the left ankle with only one constraint
VectorXd MyWorld::updateGradients()
{
    Vector3d mCTemp;
    mJ = MatrixXd::Zero(15 * 3, mSkel->getNumDofs());
	 
	
	//Intializing markers as false
    for (int i = 0; i < 15; i++)
    {

        if (aMrk[i] == false)
            continue;



        BodyNode* node = mSkel->getMarker(i)->getBodyNode();
		
		if (i ==0 || i==1)
        {
            // coded for dof = 3 as we know dof is 3 for root nodes
            mTgt[0] = mTgtB[3 * i];
            mTgt[1] = mTgtB[3 * i + 1];
            mTgt[2] = mTgtB[3 * i + 2];
			 // compute c(q)
            mCTemp = mSkel->getMarker(i)->getWorldPosition() - mTgt;
			mC[3 * i] = mCTemp[0];
            mC[3 * i + 1] = mCTemp[1];
            mC[3 * i + 2] = mCTemp[2];

            // compute J(q)
            Vector4d offset;
            offset << mSkel->getMarker(i)->getLocalPosition(), 1; // Create a vector in homogeneous coordinates

            // w.r.t ankle
            BodyNode *node = mSkel->getMarker(i)->getBodyNode();
            Joint *joint = node->getParentJoint();

            Matrix4d parentToJoint = joint->getTransformFromParentBodyNode().matrix();

            Matrix4d dR = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
            Matrix4d R = joint->getTransform(1).matrix();
            Matrix4d R1 = joint->getTransform(2).matrix();
            Matrix4d jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();

            Vector4d jCol = parentToJoint * dR * R * R1* jointToChild * offset;

            int colIndex = joint->getIndexInSkeleton(0);

            mJ.middleRows(3 * i, 3).col(colIndex) = jCol.head(3);

            dR = joint->getTransformDerivative(1);
            R = joint->getTransform(0).matrix();
            R1 = joint->getTransform(2).matrix();
            jCol = parentToJoint * R * dR * R1* jointToChild * offset;
            colIndex = joint->getIndexInSkeleton(1);

            mJ.middleRows(3 * i, 3).col(colIndex) = jCol.head(3);

            dR = joint->getTransformDerivative(2);
            R = joint->getTransform(0).matrix();
            R1 = joint->getTransform(1).matrix();
            jCol = parentToJoint * R * R1 * dR * jointToChild * offset;
            colIndex = joint->getIndexInSkeleton(2);

            mJ.middleRows(3 * i, 3).col(colIndex) = jCol.head(3);
		}

        else
        {
			//calculating index position for the parent node
			int p = i;
            while(node->getParentBodyNode()!= NULL)
            {
                int nDofs = node->getParentJoint()->getNumDofs();
                mTgt[0] = mTgtB[3 * i];
                mTgt[1] = mTgtB[3 * i + 1];
                mTgt[2] = mTgtB[3 * i + 2];

                // compute c(q)

                mCTemp = mSkel->getMarker(i)->getWorldPosition() - mTgt;

                mC[3 * i] = mCTemp[0];
                mC[3 * i + 1] = mCTemp[1];
                mC[3 * i + 2] = mCTemp[2];

              if (nDofs == 1)
              {
                    // compute J(q)
                    Vector4d offset;
                    offset << mSkel->getMarker(p)->getLocalPosition(), 1; // Create a vector in homogeneous coordinates

                    // w.r.t ankle
                    Joint *joint = node->getParentJoint();
                    Matrix4d worldToParent = node->getParentBodyNode()->getTransform().matrix();
                    Matrix4d parentToJoint = joint->getTransformFromParentBodyNode().matrix();
                    Matrix4d dR = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
                    Matrix4d jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
                    Vector4d jCol = worldToParent * parentToJoint * dR * jointToChild * offset;
                    int colIndex = joint->getIndexInSkeleton(0);

                    mJ.middleRows(3 * p, 3).col(colIndex) = jCol.head(3);
                }


                if (nDofs == 2)
                {
                    
                    // compute J(q)
                    Vector4d offset;
                    offset << mSkel->getMarker(p)->getLocalPosition(), 1; // Create a vector in homogeneous coordinates

                    // w.r.t ankle
                    Joint *joint = node->getParentJoint();
                    Matrix4d worldToParent = node->getParentBodyNode()->getTransform().matrix();
                    Matrix4d parentToJoint = joint->getTransformFromParentBodyNode().matrix();
			
                    Matrix4d dR = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
                    Matrix4d R = joint->getTransform(1).matrix();
                    Matrix4d jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
                    Vector4d jCol = worldToParent * parentToJoint * dR * R * jointToChild * offset;
                    int colIndex = joint->getIndexInSkeleton(0);

                    mJ.middleRows(3 * p, 3).col(colIndex) = jCol.head(3);

                    dR = joint->getTransformDerivative(1);
                    R = joint->getTransform(0).matrix();
                    jCol = worldToParent * parentToJoint * R * dR * jointToChild * offset;
                    colIndex = joint->getIndexInSkeleton(1);

                    mJ.middleRows(3 * p, 3).col(colIndex) = jCol.head(3);

                }


                if (nDofs == 3)
                {
                    // compute J(q)
                    Vector4d offset;
                    offset << mSkel->getMarker(p)->getLocalPosition(), 1; // Create a vector in homogeneous coordinates

                    // w.r.t ankle
                   
                    Joint *joint = node->getParentJoint();
                    Matrix4d worldToParent = node->getParentBodyNode()->getTransform().matrix();
                    Matrix4d parentToJoint = joint->getTransformFromParentBodyNode().matrix();
			
                    Matrix4d dR = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
                    Matrix4d R = joint->getTransform(1).matrix();
                    Matrix4d R1 = joint->getTransform(2).matrix();
                    Matrix4d jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
                    Vector4d jCol = worldToParent * parentToJoint * dR * R * R1* jointToChild * offset;
                    int colIndex = joint->getIndexInSkeleton(0);

                    mJ.middleRows(3 * p, 3).col(colIndex) = jCol.head(3);

                    dR = joint->getTransformDerivative(1);
                    R = joint->getTransform(0).matrix();
                    R1 = joint->getTransform(2).matrix();
                    jCol = worldToParent * parentToJoint * R * dR * R1* jointToChild * offset;
                    colIndex = joint->getIndexInSkeleton(1);

                    mJ.middleRows(3 * p, 3).col(colIndex) = jCol.head(3);

                    dR = joint->getTransformDerivative(2);
                    R = joint->getTransform(0).matrix();
                    R1 = joint->getTransform(1).matrix();
                    jCol = worldToParent * parentToJoint * R * R1 * dR * jointToChild * offset;
                    colIndex = joint->getIndexInSkeleton(2);

                    mJ.middleRows(3 * p, 3).col(colIndex) = jCol.head(3);


                }

                node = node->getParentBodyNode();
				string ndName = node->getName();

				// checking the node name of parent node and checking the specific index of that marker
				for (int k = 0; k < 15; k++){

					if (ndName == nodeName[k]){
						p = k;
					}

				}

				
            }

          

        }
        
    }

    // compute gradients
    VectorXd gradients = 2 * mJ.transpose() * mC;
    return gradients;
}

// Current code only handlse one constraint on the left foot.
void MyWorld::createConstraint(int _index) {

    aMrk[_index] = true;
    mConstrainedMarker = _index;
    
    mTgt = mSkel->getMarker(_index)->getWorldPosition();
    int iX = _index * 3;
    mTgtB[iX] = mTgt[0];
    mTgtB[iX + 1] = mTgt[1];
    mTgtB[iX + 2] = mTgt[2];
}

void MyWorld::modifyConstraint(Vector3d _delPos) {
    // update all targets with new delta
    mTgtB[3 * mConstrainedMarker] += _delPos[0];
    mTgtB[3 * mConstrainedMarker + 1] += _delPos[1];
    mTgtB[3 * mConstrainedMarker + 2] += _delPos[2];
    
}

void MyWorld::removeConstraint(int _index) {
    mConstrainedMarker = -1;
}


