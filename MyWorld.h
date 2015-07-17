#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>
#include <string>
#include "dart/dynamics/Skeleton.h"


class MyWorld {
 public:
    MyWorld();
    virtual ~MyWorld();
    dart::dynamics::Skeleton* getSkel() {
        return mSkel;
    }

    void solve();
    void createConstraint(int _index);
    void modifyConstraint(Eigen::Vector3d _delPos);
    void removeConstraint(int _index);

 protected:
    Eigen::VectorXd updateGradients();
    bool aMrk[15];
	
    dart::dynamics::Skeleton *mSkel;
    Eigen::VectorXd mC;
    Eigen::MatrixXd mJ;
    Eigen::Vector3d mTgt; // The target location of the constriant
    Eigen::VectorXd mTgtB;
    int mConstrainedMarker; // The index of the constrained marker
};

#endif
