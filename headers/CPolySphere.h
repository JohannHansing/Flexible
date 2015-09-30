/*
 * CPolySphere.h
 *
 *  Created on: Mar 11, 2015
 *      Author: johann hansing
 */

#ifndef CPOLYSPHERE_H_
#define CPOLYSPHERE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>
#include <Eigen/Dense>

using namespace std;

class CPolySphere {
public:
    CPolySphere();
    CPolySphere( Eigen::Vector3d PolySpherePos);
        
    // Public variables
    Eigen::Vector3d pos;  //fixed length array for position, this allocates less memory than a (dynamic) vector 
    Eigen::Vector3d f_mob; // Vector to store the mobility force on the polysphere
    Eigen::Vector3d f_sto; // Stores random vector for displacement of polySPheres

};



#endif /* CPOLYSPHERE_H_ */
