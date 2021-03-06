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
    Eigen::Vector3d pos_abs;  // pos_abs is the absolute position of the network bead relative to each other. It does not change, when the beads leave the simulation box
    Eigen::Vector3d pos_pbc;  // pos_pbc is the position of network bead inside the simulation box. It is shifted, if the sphere crosses the border of the sim box.
    Eigen::Vector3d f_mob; // Vector to store the mobility force on the polysphere
    Eigen::Vector3d f_sto; // Stores random vector for displacement of polySPheres
    std::vector<Eigen::Vector3d> image_corr; // Correction vector for polySpheres that have a right neighbor in another simulation box (periodic b.c.)
    int n_rn; // number of right neigbors for spring interaction. Initialized to zero and incremented by 1 for each rightneigbor added.
    bool image[3];
    int i_rn[3]; // indexes of right neighbors are stored here. 
    std::vector<int> i_LJ; // indexes of extra LJ neigbors are stored here.
    
    
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
    
    void addextraLJneighbor(int LJneighbor){
        i_LJ.push_back(LJneighbor);
    }

};



#endif /* CPOLYSPHERE_H_ */
