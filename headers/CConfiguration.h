#ifndef CCONFIGURATION_H_
#define CCONFIGURATION_H_

#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string>
#include <vector>
#include <array>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include "CPolySphere.h"
#include "parameter_structs.h"


class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
private:
    //SCALING
    double _timestep;         //This is the RESCALED timestep! timestep = dt * kT / (frictionCoeffcient * particlesize)
    double _mu_sto;
    double _mu_sto_poly;
    double _epsilon;


    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    //LJ Potential
    double _epsilonLJ;
    double _uLJ;  // LJ potential of tracer
    //SPRING Potential
    double _kappaSP;
    double _r0SP;
    double _uspring;
    //BEND Potential
    bool _bendPot;
    double _kappaBend;
    double _ubend; 
    
    


    //COUNTERS AND INIT VALUES
    int _boxnumberXYZ[3];           //counter to calculate the actual position of the particle
    Eigen::Vector3d _prevpos;           //Stores previous particle position before particle is moved.
    string _testcue;

    //Particle parameters
    Eigen::Vector3d _ppos;    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _pradius;     //particle size is most relevant for scaling! (default = 1)
    double _upot;
    Eigen::Vector3d _f_mob;   //store mobility and stochastic force
    Eigen::Vector3d _f_sto;
    Eigen::Vector3d _startpos;
    
    //Lattice parameters
    std::vector<CPolySphere> _polySpheres;
    int _n_cellsAlongb;
    int _N_polySpheres;
    int _edgeParticles;
    int _N_cellParticles;
    double _boxsize;          // ALWAYS define boxsize through particlesize due to scaling!
    double _polyrad;
    
    //Storage Matrices for distances and distance vectors
    std::vector<std::vector<double> > _Mrabs; 
    std::vector<std::vector<Eigen::Vector3d> > _Mrvec;
        
	
    // Bead-Spring interaction
    std::vector<std::vector<bool> > _springMatrix; // Bool parameter is true for spring-connected neighbors
    int _N_springInteractions;
    //SemiFlexiblePolymer
    std::vector<std::array<int, 2> > _MspringTupel;
    std::vector<std::array<int, 3> > _MbendTupel;

    //MISC
    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".
    FILE* m_traj_file;



public:
    CConfiguration();
    CConfiguration(double timestep, model_param_desc modelpar, sim_triggers triggers);
    void updateStartpos();
    void makeStep();
    void checkBoxCrossing();
    void calcStochasticForces();
    void calcMobilityForces();
    void saveXYZTraj(string name, const int& move, string flag);
    void saveGRO(string foldername, const int& move);
    void saveCoordinates(std::ostream& trajectoryfile, unsigned int stepcount);
    Eigen::Vector3d minImage(Eigen::Vector3d rij);
	std::vector<double> getppos();
    double getUpot(){ return _upot; }
    string getTestCue(){ return _testcue; };



private:
    void setRanNumberGen(double seed);
    void initConstants(){
        // Function to init constants, so this doesn't clutter the cpp file
        _edgeParticles =  (int) ( _boxsize/_n_cellsAlongb/(2*_polyrad) + 0.0001);
        _N_cellParticles = 3 * _edgeParticles - 2;
        _N_polySpheres = _N_cellParticles * pow(_n_cellsAlongb,3);
        _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
        _mu_sto_poly = sqrt( 2 * _timestep * _pradius/_polyrad );
        _epsilonLJ = 1;  // in kT 
        _upot = 0;
        _r0SP = _boxsize/_n_cellsAlongb/_edgeParticles; //equilibrium distance for spring potential
        cout << "TODO: Adjust bending and spring parameters! Consider Schlagberger2006 and Metzler paper." << endl;
    }
    
    
    int makeIndex(int i, int nx, int ny, int nz){ // Transforms a NxKxM matrix index (i,j,k) to an array index ind = i + j*K + l*M
        assert(i < _N_cellParticles  && "**** Index i out of range in makeIndex()!");
        if (nx == -1) nx += _n_cellsAlongb;
        if (ny == -1) ny += _n_cellsAlongb;
        if (nz == -1) nz += _n_cellsAlongb;
        if (nx == _n_cellsAlongb) nx = 0;
        if (ny == _n_cellsAlongb) ny = 0;
        if (nz == _n_cellsAlongb) nz = 0;
        return i + nx * _N_cellParticles + ny * _n_cellsAlongb * _N_cellParticles + nz * pow(_n_cellsAlongb,2) * _N_cellParticles;
    }
    
    //POTENTIALS
    void addLJPot(const double r, double& U, double& Fr, const double r_steric){
        //Function to calculate the Lennard-Jones Potential
        if ( r < 1.122462 * r_steric ){ // steric parameter r_steric is the added radii of both Lennard-Jones particles
            double  por6 = pow((r_steric / r ), 6);      //por6 stands for "p over r to the power of 6" . r_steric is the total steric parameter
            U += 4 * _epsilonLJ * ( por6*por6 - por6 + 0.25 );
            Fr +=  24  * _epsilonLJ / ( r * r ) * ( 2 * por6*por6 - por6 );
        }
    }
    
    void addSpringPot(const double& r, double &U, double &Fr) { 
        U += (_kappaSP * pow(r - _r0SP, 2)); 
        Fr += (_kappaSP * 2 * (_r0SP/r - 1));
    }
    
    
    void calcBendPot(Eigen::Vector3d vec_r12, Eigen::Vector3d vec_r32, double r12, double r32, double& U, Eigen::Vector3d& fvec1, Eigen::Vector3d& fvec3){
        // Function that returns the bending potential between three particles, where partice 2 is in the middle 
        Eigen::Vector3d vec_n = vec_r12.cross(vec_r32);
        double r12r32 = r12*r32;
       // cout << "theta = " << 180 / 3.14 * acos(vec_r12.dot(vec_r32)/r12r32) << endl;
        fvec1 =  _kappaBend *  vec_r12.cross(vec_n) / (r12*r12*r12r32);
        fvec3 =  _kappaBend *  vec_n.cross(vec_r32) / (r32*r32*r12r32);
        U     =  _kappaBend *  ( 1 + vec_r12.dot(vec_r32)/r12r32 );
    }
    
    
    
    
    
    // LATTICE INITIALIZATION FUNCTIONS
    
    /******* TEMPLATE
        - overwrite neccessary parameters like _r0SP
        - Fill edgeParticles Vector
        - resize storage matrices/vectors
        - init 2p and 3p interaction vectors / matrices
    */
    
    void initSemiFlexibleLattice(){
        assert(_edgeParticles > 1 && "**** Edgeparticles need to be more than 1 for semiflexible lattice!");
        _polySpheres.clear();
        Eigen::Vector3d nvec;
        Eigen::Vector3d ijk;
    
        // init PolySpheres
        std::vector<Eigen::Vector3d> zeroPos( _N_cellParticles , Eigen::Vector3d::Zero() );
    	// store the edgeParticle positions in first cell in zeroPos
        assert(((_boxsize/_n_cellsAlongb) / _edgeParticles > 2*_pradius)  && "**** Error: Overlap between polySpheres!");
    	for (int i = 1; i < _edgeParticles; i++){
    		double tmp = i * (_boxsize/_n_cellsAlongb) / _edgeParticles;
    		zeroPos[i](0) = tmp;
    		zeroPos[i + (_edgeParticles - 1)](1) = tmp;  
    		zeroPos[i + 2 * (_edgeParticles - 1)](2) = tmp;
    	}
    	for (int nz=0; nz < _n_cellsAlongb; nz++){
    	    for (int ny=0; ny < _n_cellsAlongb; ny++){
                for (int nx=0; nx < _n_cellsAlongb; nx++){
                    nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                    for (unsigned int i = 0; i < _N_cellParticles; i++){
                        _polySpheres.push_back( CPolySphere( zeroPos[i] + nvec * _boxsize/_n_cellsAlongb ) );
                	}
                }
            }
    	}
        assert(_N_polySpheres == _polySpheres.size()  && "**** Size of _polySpheres vector is incorrect, i.e. not equal to _N_polySpheres!");
        
        //resize storage containers for distances and vectors
        _Mrvec.resize(_N_polySpheres+1);
        _Mrabs.resize(_N_polySpheres+1);
        for (int i=0; i<_N_polySpheres+1; i++){
            _Mrvec[i].resize(_N_polySpheres+1);
            _Mrabs[i].resize(_N_polySpheres+1);
        }
    
        // Init two-particle spring and three-particle bending interaction Vectors
        _N_springInteractions = (3 * _edgeParticles) * pow(_n_cellsAlongb,3);
        _MspringTupel.resize(_N_springInteractions);
        _MbendTupel.resize(_N_springInteractions);
        // Init two-partice spring Matris to all false
        _springMatrix.resize(_N_polySpheres);
        for (int i=0;i<_N_polySpheres;i++){
            _springMatrix[i].resize(_N_polySpheres,false);
        }
    
        int counter = 0;
        int ne = _edgeParticles - 1; // "extra particles" between 0 particles. This I use in my paper version for the code
    	for (int nz=0; nz < _n_cellsAlongb; nz++){
    	    for (int ny=0; ny < _n_cellsAlongb; ny++){
                for (int nx=0; nx < _n_cellsAlongb; nx++){
                    nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                    for (int i = 0; i < _N_cellParticles; i++){
                        // particle at 0 position of cell - 3 right neighbors and 3 bending interactions
                        if (i == 0){
                            for (int t=0; t<3; t++){
                                int a = 0, b = 0, c = 0;
                                if (t==0) a=-1;
                                else if (t==1) b=-1;
                                else c=-1;
                                _springMatrix[makeIndex(i,nx,ny,nz)][makeIndex(1+ne*t,nx,ny,nz)] = true;
                                _MspringTupel[counter][0] = makeIndex(i,nx,ny,nz);
                                _MspringTupel[counter][1] = makeIndex(1+ne*t,nx,ny,nz);
                                _MbendTupel[counter][0] = makeIndex(ne+ne*t,nx+a,ny+b,nz+c);
                                _MbendTupel[counter][1] = _MspringTupel[counter][0];
                                _MbendTupel[counter][2] = _MspringTupel[counter][1];
                                counter++;
                            }
                        }
                        else if (i % ne == 0 && _edgeParticles == 2){ // Special case for _edgeParticles == 2
                            int a = 0, b = 0, c = 0;
                            if (i == ne ) a = 1;
                            else if (i == 2*ne ) b = 1;
                            else if (i == 3*ne ) c = 1;
                            else {
                                cout << "Error in CConfiguration initSemiFlexLattice() !!!!" << endl;
                                abort();
                            }
                            _springMatrix[makeIndex(i,nx,ny,nz)][makeIndex(0,nx+a,ny+b,nz+c)] = true;
                            _MspringTupel[counter][0] = makeIndex(i,nx,ny,nz);
                            _MspringTupel[counter][1] = makeIndex(0,nx+a,ny+b,nz+c);
                            _MbendTupel[counter][0] = makeIndex(0,nx,ny,nz);
                            _MbendTupel[counter][1] = _MspringTupel[counter][0];
                            _MbendTupel[counter][2] = _MspringTupel[counter][1];
                            counter++;
                        }
                        else if (i % ne == 0){ // if i = ne, 2ne or 3ne
                            int a = 0, b = 0, c = 0;
                            if (i == ne ) a = 1;
                            else if (i == 2*ne ) b = 1;
                            else if (i == 3*ne ) c = 1;
                            else {
                                cout << "Error in CConfiguration initSemiFlexLattice() !!!!" << endl;
                                abort();
                            }
                            _springMatrix[makeIndex(i,nx,ny,nz)][makeIndex(0,nx+a,ny+b,nz+c)] = true;
                            _MspringTupel[counter][0] = makeIndex(i,nx,ny,nz);
                            _MspringTupel[counter][1] = makeIndex(0,nx+a,ny+b,nz+c);
                            _MbendTupel[counter][0] = makeIndex(i-1,nx,ny,nz);
                            _MbendTupel[counter][1] = _MspringTupel[counter][0];
                            _MbendTupel[counter][2] = _MspringTupel[counter][1];
                            counter++;
                        }
                        else  {
                            int a = i-1;
                            if ((i == ne+1) || (i == 2*ne+1)) a = 0;
                            _springMatrix[makeIndex(i,nx,ny,nz)][makeIndex(i+1,nx,ny,nz)] = true;
                            _MspringTupel[counter][0] = makeIndex(i,nx,ny,nz);
                            _MspringTupel[counter][1] = makeIndex(i+1,nx,ny,nz);
                            _MbendTupel[counter][0] = makeIndex(a,nx,ny,nz);
                            _MbendTupel[counter][1] = _MspringTupel[counter][0];
                            _MbendTupel[counter][2] = _MspringTupel[counter][1];
                            counter++;
                        }
                    }
                }
            }
        }
        assert((counter == _N_springInteractions) && "**** final counter Not equal size as _N_springInteractions");
        //printSpringMatrix();
    }
    

    void initZhouPolyspheres(){
        // initialize the PolymerSphere positions
        _polySpheres.clear();
        Eigen::Vector3d nvec;
    	for (int nz=0; nz < _n_cellsAlongb; nz++){
    	    for (int ny=0; ny < _n_cellsAlongb; ny++){
                for (int nx=0; nx < _n_cellsAlongb; nx++){
                    nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                    _polySpheres.push_back( CPolySphere(  nvec * _boxsize/_n_cellsAlongb ) );
                }
            }
    	}
    }
    
    void initZhouInteractionMatrix(){
        assert(_N_polySpheres != 0 && "**** _N_polySpheres is zero. Can't initialize interaction Matrix!");
        Eigen::Vector3i nvec;
        Eigen::Vector3i ijk;
        // Init two-partice spring Matris to all false
        _springMatrix.resize(_N_polySpheres);
        for (int i=0;i<_N_polySpheres;i++){
            _springMatrix[i].resize(_N_polySpheres,false);
        }
    
        // Change appropriate spring Matrix parts to true
    
        // init cubic lattice w spheres in corners
        _r0SP = _boxsize/_n_cellsAlongb; //equilibrium distance for spring potential Leave this here since it belongs to cubic lattice interactionMatrix
    	for (int nz=0; nz < _n_cellsAlongb; nz++){
    	    for (int ny=0; ny < _n_cellsAlongb; ny++){
                for (int nx=0; nx < _n_cellsAlongb; nx++){
                    nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                    int nvec_to_m = nvec(0)+nvec(1)*_n_cellsAlongb+nvec(2)*pow(_n_cellsAlongb,2);
                    int ijk_to_m = 0;
                    for (int n=0;n<3;n++){
                        ijk = nvec;
                        ijk(n) = nvec(n) - 1;  //neighbor in -n direction
                        if (ijk(n) == -1) ijk(n) += _n_cellsAlongb;  //periodic BC
                        ijk_to_m = ijk(0) + ijk(1) * _n_cellsAlongb + ijk(2)*pow(_n_cellsAlongb,2);
                        _springMatrix[nvec_to_m][ijk_to_m] = true;
        
                        ijk = nvec;
                        ijk(n) = nvec(n) + 1;  //neighbor in +n direction
                        if (ijk(n) == _n_cellsAlongb ) ijk(n) -= _n_cellsAlongb;  //periodic BC
                        ijk_to_m = ijk(0) + ijk(1) * _n_cellsAlongb + ijk(2)*pow(_n_cellsAlongb,2);
                        _springMatrix[nvec_to_m][ijk_to_m] = true;
                    }
                }
            }
    	}
    }
    
    
    //MISC
    template<typename T>
    string toString(const T& value){
        ostringstream oss;
        oss << value;
        return oss.str();
    }
    void printSpringMatrix(){
        cout << "*********** SPRING MATRIX ************" << endl;
        for (unsigned int i = 0; i < _MspringTupel.size() ; i++) {
            cout << "[ ";
            for (unsigned int j = 0; j < 2 ; j++) {
                cout << _MspringTupel[i][j] << " , ";
            }
            cout <<" ]"<< endl;
        }
        cout << "\n\n\n*********** BENDING MATRIX ************" << endl;
        for (unsigned int i = 0; i < _MbendTupel.size() ; i++) {
            cout << "[ ";
            for (unsigned int j = 0; j < 3 ; j++) {
                cout << _MbendTupel[i][j] << " , ";
            }
            cout <<" ]"<< endl;
        }
    }
    

};



#endif /* CCONFIGURATION_H_ */
