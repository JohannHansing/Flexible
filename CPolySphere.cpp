#include "headers/CPolySphere.h"

using namespace std;

CPolySphere::CPolySphere(){
}

CPolySphere::CPolySphere( Eigen::Vector3d PolySpherePos){
    pos_abs = PolySpherePos;
    pos_pbc = PolySpherePos;
    f_mob = Eigen::Vector3d::Zero();
    f_sto = Eigen::Vector3d::Zero();
    image_corr.resize(3, Eigen::Vector3d::Zero()); // In the case of n_edge =1,i.e. SingleSpheres, there can be three image_corr vectors
    n_rn = 0;
}
