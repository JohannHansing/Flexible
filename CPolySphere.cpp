#include "headers/CPolySphere.h"

using namespace std;

CPolySphere::CPolySphere(){
}

CPolySphere::CPolySphere( Eigen::Vector3d PolySpherePos){
    pos = PolySpherePos;
    f_mob = Eigen::Vector3d::Zero();
    f_sto = Eigen::Vector3d::Zero();
}
