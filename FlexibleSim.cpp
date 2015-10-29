/*
 * FlexibleSim.cpp
 *
 *  Created on: June 3, 2015
 *      Author: Johann Hansing
 */

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

#include "headers/CAverage.h"
#include "headers/CConfiguration.h"
#include "headers/parameter_structs.h"


using namespace std;


// Create global instances of structs
struct sim_param_desc _simpar;
struct model_param_desc _modelpar;
struct file_desc _files;
struct sim_triggers _triggers;


//Function declarations //TODO put all this in header Flexible.h!
void createDataFolder(string cue);
void parameterFile(string cue);
void parameterFileAppend(double executiontime);

template<typename T>
string toString(const T& value);

template <typename T, size_t N>
inline
size_t sizeOfArray( const T(&)[ N ] );

int main(int argc, const char* argv[]){
	// measure runtime of the simulation
	clock_t start, end;

    // INPUT PARAMETERS:
	int boolpar = 0;

    //TRIGGERS:
    //_triggers.bendPot = (strcmp(argv[1] , "true") == 0 ) ;
	// Checking for correct structure of input arguments
	for (int k= 0; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
	for (int b_i=1; b_i<=boolpar; b_i++){
		if ((strcmp(argv[b_i] , "true") == 1 )  && (strcmp(argv[b_i] , "false") == 1 )){
			cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
			exit(1);
		}
	}

    //NUMBERS
    _simpar.runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    _simpar.timestep = atof( argv[boolpar+2] );
    _simpar.simtime = atoi( argv[boolpar+3] );                   // simulation time
    _simpar.instantvalues = 200;
    _simpar.steps = _simpar.simtime/_simpar.timestep;
    _simpar.saveInt = _simpar.steps/_simpar.instantvalues;

    _modelpar.polymersize = atof( argv[boolpar+4] );
    _modelpar.particlesize = atof( argv[boolpar+5] );
    _modelpar.n_cells = atof( argv[boolpar+6] );
    _modelpar.urange = atof( argv[boolpar+7] );
    _modelpar.ustrength = atof( argv[boolpar+8] );
    _modelpar.kspring = atof( argv[boolpar+9] );
    _modelpar.kbend = atof( argv[boolpar+10] );
    _modelpar.boxsize = _modelpar.n_cells * 10.;


    int instValIndex;                             //Counter for addInstantValue during runs-loop
    const int trajout = 10/_simpar.timestep;    // trajout * timestep = 10!


    //initialize instance of configuration
	CConfiguration conf;
	try{
        conf = CConfiguration(_simpar.timestep, _modelpar, _triggers, _files);
    }
	catch(int e){
		cout << "An exception occurred. Exception Nr. " << e << '\n';
		return 1;
	}
    //Create data folders and print location as string to string "folder"
    createDataFolder(conf.getTestCue());
    if (!conf.printGroFile(_files.folder)){   //This also throws an exception. Hence, it can be moved back into conf without problems
		return 1;
	}


    //create file to save the trajectory
//    string traj_file = _files.folder + "/Coordinates/single_traj.xyz";
//    if (writeTrajectory) conf.saveXYZTraj(traj_file,0,"w");

    // Initialize counters and stuff
    start = clock();
    cout << "Starting Simulation!" << endl;

    unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((_files.folder + "/Coordinates/trajectory.txt").c_str());
	_files.xtc_filename = "TEST.xtc";    // TODO xtc
    _files.xd = xdrfile_open((_files.folder + "/Coordinates/" + _files.xtc_filename).c_str(), "w");
	if (!_files.xd) {
      std::cout << "Error: Could not open trajectory file " << _files.xtc_filename << " for writing." << std::endl;
      return 1;
    }

    //create .xyz file to save the trajectory for VMD
    string traj_file = _files.folder + "/Coordinates/single_traj.xyz";
    conf.saveXYZTraj(traj_file,0,"w");

    // Write parameter file parameters.txt
    parameterFile(conf.getTestCue());




// ************** START OF SIMULATION RUNS-LOOP *****************
    for (int l = 0; l<_simpar.runs; l++){
        instValIndex = 0;
        //TODO conf.updateStartpos();

// -------------- iteration step-loop ------------
        for (int i = 0; i < _simpar.steps; i++){

            conf.calcStochasticForces();
            conf.calcMobilityForces();

            stepcount++;
            conf.makeStep();    //move particle at the end of iteration

            // steric hard-sphere interaction
            // while (includeSteric && conf.testOverlap()){
//                 conf.moveBack();
//                 conf.calcStochasticForces();
//                 conf.makeStep();
//             }


            // Write trajectory to trajectoryfile
            if (stepcount%trajout == 0) {
                conf.saveCoordinates(trajectoryfile, stepcount);
				conf.save_traj_step(_files.xd,i);  //TODO xtc
            }
            if (((i+1)%100 == 0) && (l == 0)){       //Save the first trajectory to file
                conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }
        }
        if (l==0) conf.saveXYZTraj(traj_file, _simpar.steps, "c"); // Close XYZ traj_file

    }//----------END OF RUNS-LOOP ----------------



	cout << "Simulation Finished" << endl;
	end = clock();
	double runtime = (double)((end-start)/(CLOCKS_PER_SEC));
	cout << runtime << " seconds runtime." << endl;

    parameterFileAppend(runtime);

	trajectoryfile.close();
	xdrfile_close(_files.xd); //TODO xtc


    return 0;
}









//--------------------------------------------------------------------------
//**************************************************************************
//--------------------------------------------------------------------------




void createDataFolder(string testcue){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[5];
    //sprintf(range, "%.3f", _modelpar.potRange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data";
    if (!testcue.empty()) folder = folder + "/test/" + testcue;
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
 //   parameterFile << "k " << _modelpar.potRange << endl;
 //   parameterFile << "U_0 " << _modelpar.potStrength << endl;

    parameterFile.close();
}

void parameterFileAppend(double executiontime){
    // Appends parameters to parameters file
    ofstream parameterFile;
    parameterFile.open((_files.folder + "/parameters.txt").c_str(), std::ios_base::app);
    parameterFile << "ExecutionTime " << executiontime << " s" << endl;
    parameterFile.close();
}


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
