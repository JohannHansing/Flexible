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
    std::vector<Eigen::Vector3d> image_corr; // Correction vector for polySpheres that have a right neighbor in another simulation box (periodic b.c.)
    int n_rn; // number of right neigbors for spring interaction. Initialized to zero and incremented by 1 for each rightneigbor added.
    bool image[3];
    int i_rn[3]; // indexes of right neighbors are stored here. 
    
    
    void addRightNeighbor(int rightneighborIndex, int myindex, Eigen::Vector3d pbc_shift = Eigen::Vector3d::Zero()){
        if ((pbc_shift != Eigen::Vector3d::Zero()) && (rightneighborIndex < myindex)){
            image_corr[n_rn] = pbc_shift;
            image[n_rn] = true;
            //cout << "shift for particle " << myindex << " with neighbor " << rightneighborindex  << endl;
        }
        else image[n_rn] = false;
        
        i_rn[n_rn] = rightneighborIndex;
        n_rn++;
    }

};



#endif /* CPOLYSPHERE_H_ */
