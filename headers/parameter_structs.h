//
//  paramater_structs.h
//
//
//  Created by Richard Schwarzl on 28/04/15.
//
//

#ifndef ____paramater_structs__ //checks if already defined by other header file
#define ____paramater_structs__


#include <string.h>
#include <Eigen/Dense>
#include "xdrfile.h"
#include "xdrfile_xtc.h"

struct sim_param_desc {
    int runs;
    double timestep;                      // Internal time step
    unsigned int simtime;
    int instantvalues;
    unsigned int steps;
    unsigned int  saveInt;
};

struct sim_triggers {
    // TEST CUE to modify the directory the output data is written to!!
    string _testcue = "FillEdge";
    //bool bendPot;
};

struct file_desc {
    // TODO MAYBE
    string folder;
    XDRFILE *xd;                    // handle to trajectory file
    std::string xtc_filename;       // trajectory filename
};


struct model_param_desc {
    double polymersize;
    double particlesize;
    double boxsize;
    double urange;
    double ustrength;
    int n_cells;  // Number of cells along one edge of the simulation box
    double kspring;
    double kbend;
};



 struct simbox_desc {
     // TODO MAYBE
//     unsigned int numberOfAtoms;              // Number of atoms in the box
//     vec box;                        // box lengths x,y,z
//     double minHalfBoxLength;        // maximum distance between to particles due to minimal image convention
//     matrix box_matrix;              // box vectors for xtc writing
//     double volume;                  // Volume of the box
//     double density;                 // Density of the System
 };

// struct state_desc {
//     std::vector<Eigen::Vector3d> allPositions;                              // Position of atoms in the box
//     std::vector<std::vector<Eigen::Vector3d> > allDistanceUnitV;                // Matrix of the distance vector between all pairs of atoms.
//     std::vector<std::vector<double> > allScalarDistances;       // Matrix of the scalar distances between all pairs of atoms.
//     std::vector<std::vector<double> > allScalarForces;          // Tensor of the forces between all pairs of atoms.
//   //  std::vector<std::vector<Eigen::Vector3d> > allForces;                   // Tensor of the forcevectors between all pairs of atoms.
//     //std::vector<Eigen::Vector3d> allSummedForces;                           // Array that contains the forces on every atom respectively.
//     //std::vector<std::vector<Eigen::Vector3d> > allLJForces;                 // Tensor of the Lennard Jones forcevectors between all pairs of atoms.
//     // The following will be needed for addPullingForce():
//     std::vector<std::vector<double> > allInitScalarDistances;   // Initial scalar distances from init.gro
//     double scalarPullForce;                                     // scalar amount of pulling force
//     double pullExtension;                                       // pull extension
// };
//
// struct ener_desc {
//     double internalEnergy;
//     double pairVirial;
//     double pressure;
//
//     std::vector<std::vector<double> > virialTensor;
//     std::vector<std::vector<double> > pressureTensor;
// };

#endif /* defined(____paramater_structs__) */ //end check
