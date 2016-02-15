/*
 * FlexibleSim.h
 *
 *  Created on: November 3, 2015
 *      Author: Johann Hansing
 */

#ifndef FLEXIBLE_H_
#define FLEXIBLE_H_



#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <math.h>
#include <boost/filesystem.hpp>
#include "xdrfile.h"
#include "xdrfile_xtc.h"

#include "CAverage.h"
#include "CConfiguration.h"
#include "parameter_structs.h"


// Create global instances of structs
struct sim_param_desc _simpar;
struct model_param_desc _modelpar;
struct file_desc _files;
struct sim_triggers _triggers;



// void read_run_params(const std::string& filename) {
//     CSimpleIniA ini;
//     ini.SetUnicode();
//     if (ini.LoadFile(filename.c_str()) < 0) {
//       std::cout << "Could not open run parameters file " << filename << std::endl;
//       abort();
//     }
//     run_params.dt = atof(ini.GetValue("", "dt", "0.002"));
//     run_params.numberOfSteps = atoi(ini.GetValue("", "nsteps", "1000"));
//     run_params.nPrint = atoi(ini.GetValue("", "nstlog", "10"));
//     run_params.nSave = atoi(ini.GetValue("", "nstxtcout", "1"));
//     run_params.prec = atof(ini.GetValue("", "xtc-precision", "1000"));
//     run_params.T = atof(ini.GetValue("", "ref-t", "300"));
//     run_params.gamma = atof(ini.GetValue("", "bd-fric", "0"));
//     run_params.cutOff = atof(ini.GetValue("", "rvdw", "0.9"));
//     run_params.pull_k = atof(ini.GetValue("", "pull_k", "1."));
//     run_params.pull_vel = atof(ini.GetValue("", "pull_vel", "1."));
//     run_params.pull_bead_l = atoi(ini.GetValue("", "pull_bead_l", "0"));
//     run_params.pull_bead_r = atoi(ini.GetValue("", "pull_bead_r", "1"));
//     run_params.bead_r = atof(ini.GetValue("", "bead_r", "40."));
//
//     std::cout << "-------------- Read run file ----------------------" << std::endl;
//
//     simbox.micStatus = atoi(ini.GetValue("", "mic_status", "1"));
//     if (simbox.micStatus == 0) {
//       std::cout << "minimal image convention is: \t\t\t DISABLED" << std::endl;
//     } else {
//       if (simbox.micStatus == 1) {
//         std::cout << "minimal image convention is: \t\t\t ACTIVE" << std::endl;
//       } else {
//         std::cout << simbox.micStatus << " is not a valid value for mic_status." << std::endl;
//         abort();
//       }
//     }
//
//     run_params.noise_status = atoi(ini.GetValue("", "noise_status", "1"));
//     if (run_params.noise_status == 0) {
//       std::cout << "noise term (in brownian dynamics integrator): \t DISABLED" << std::endl;
//     } else {
//       if (run_params.noise_status == 1) {
//         std::cout << "noise term (in brownian dynamics integrator): \t ACTIVE" << std::endl;
//       } else {
//         std::cout << run_params.noise_status << " is not a valid value for noise_status." << std::endl;
//         abort();
//       }
//     }
//
//     run_params.shear_rate = atof(ini.GetValue("", "shear_rate", "0."));
//     if (run_params.shear_rate == 0.) {
//       std::cout << "shear is: \t\t\t\t\t DISABLED" << std::endl;
//     } else {
//       std::cout << "shear is set to: \t\t\t\t " << run_params.shear_rate << " [1/ps]" << std::endl;
//     }
//
//     if (run_params.shear_rate>0 && simbox.micStatus==1) {
//       std::cout << "minimal image is not yet implemented for shear flow conditions." << std::endl;
//       abort();
//     }
//
//     /*   currently not used: calculate gamma from mass
//     run_params.tau_t = atof(ini.GetValue("", "tau-t", "0"));
//     if (run_params.gamma == 0) {
//         if (run_params.tau_t != 0) {
//             run_params.gamma = forcefield.mass/run_params.tau_t;
//         } else {
//             std::cout << "State either bd-fric oder tau-t" << std::endl;
//             abort();
//         }
//     } else {
//         if (run_params.tau_t == 0) {
//             std::cout << "bd-fric and tau-t given, tau-t will be irgnored" << std::endl;
//         }
//     }
//     */
//
//     run_params.kT = run_params.T * 8.31e-3;
//     //forcefield.internalEnergyShift = (forcefield.c12 * std::pow(run_params.cutOff,-12)) - (forcefield.c6 * std::pow(run_params.cutOff,-6));
//
//     if (run_params.numberOfSteps < 10) {
//         run_params.nPrint = 1;
//     }
// }


template<typename T>
string toString(const T& value){
    ostringstream oss;
    oss << value;
    return oss.str();
}

template <typename T, size_t N>
inline
size_t sizeOfArray( const T(&)[ N ] )
{
  return N;
}



void createDataFolder(string testcue){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[5];
    //sprintf(range, "%.3f", _modelpar.potRange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data";
    if (!testcue.empty()) folder += testcue;
    folder += "/n_cells" + toString(_modelpar.n_cells);
    folder = folder
            + "/dt" + toString(_simpar.timestep)
            + "/t" + toString(_simpar.simtime)
            + "/kb" + toString(_modelpar.kbend)
            + "/ks" + toString(_modelpar.kspring)
            + "/a" + toString(_modelpar.polymersize)
            + "/b" + toString(_modelpar.boxsize)
            + "/p" + toString(_modelpar.particlesize);
        //    + "/k" + range
        //    + "/u" + toString(potStrength);
    boost::filesystem::create_directories(folder);
    boost::filesystem::create_directory(folder + "/InstantValues");
    boost::filesystem::create_directory(folder + "/Coordinates");
	cout << "Writing data to folder:\n" << folder << endl;
    _files.folder = folder;
}


void parameterFile(string testcue){
    //Creates a file where the simulation settings are stored
    ofstream parameterFile;
    parameterFile.open((_files.folder + "/parameters.txt").c_str());

    // Print time and date
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    parameterFile << "date " << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday
         << endl;
    parameterFile << "starttime " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << endl;
    parameterFile << "Sim_dir " << _files.folder << endl;
    if (!testcue.empty()) parameterFile << "Test cue " << testcue << endl;
    parameterFile << "p " << _modelpar.particlesize << endl;
    parameterFile << "n " << _modelpar.n_cells << endl;
    parameterFile << "b " << _modelpar.boxsize << endl;
    parameterFile << "dt " << _simpar.timestep << endl;
    parameterFile << "runs " << _simpar.runs << endl;
    parameterFile << "steps " << _simpar.steps << endl;
    parameterFile << "time " << _simpar.timestep*_simpar.steps << endl;
    parameterFile << "a " << _modelpar.polymersize << endl;
    parameterFile << "ks " << _modelpar.kspring << endl;
    parameterFile << "kb " << _modelpar.kbend << endl;
    parameterFile << "n_cells " << _modelpar.n_cells << endl;
    parameterFile << "n_edge " << _modelpar.n_edge << endl;


    parameterFile.close();
}

void parameterFileAppend(double executiontime){
    // Appends parameters to parameters file
    ofstream parameterFile;
    parameterFile.open((_files.folder + "/parameters.txt").c_str(), std::ios_base::app);
    parameterFile << "ExecutionTime " << executiontime << " s" << endl;
    parameterFile.close();
}



#endif /* FLEXIBLE_H_ */
