/*
 * CConfiguration.cpp
 *
 *  Created on: June 4, 2015
 *      Author: Johann Hansing
 */


#include "headers/CConfiguration.h"


using namespace Eigen;
using namespace std;

CConfiguration::CConfiguration(){
}


CConfiguration::CConfiguration(double timestep, model_param_desc modelpar, sim_triggers triggers, file_desc files){
    _pradius = modelpar.particlesize/2;   //_pradius is now the actual radius of the particle. hence, I need to change the definition of the LJ potential to include (_pradius + _polyrad)   -  or maybe leave LJ pot out
    _polyrad = modelpar.polymersize / 2;   //This is needed for testOverlap for steric and HI stuff !!
	_boxsize = modelpar.boxsize;
    _n_cellsAlongb = modelpar.n_cells;
    _kappaSP = modelpar.kspring; // in kT/a^2
    _kappaBend = modelpar.kbend; // in kT
    _timestep = timestep;

    for (int i = 0; i < 3; i++){
        _ppos(i) = _boxsize/_n_cellsAlongb/2;
        _boxnumberXYZ[i] = 0;
        _startpos(i) = _ppos(i);
    }

    initConstants();

    // init random number generator
    setRanNumberGen(0);
    //initPolyspheres();
    //initInteractionMatrix(); // This comes only AFTER initPolyspheres
    initSemiFlexibleLattice();
    cout << "NOTE: Implement periodic boundary conditions for the polySpheres. I.e. checkBoxCrossing function! " << endl;
    // TEST CUE to modify the directory the output data is written to!!
    _testcue = "FillEdge";
    //if ( _n_cellsAlongb != 1 ) _testcue = "n" + toString(_n_cellsAlongb);

}




void CConfiguration::calcMobilityForces(){
    Vector3d vec_rij(3);
    double rij = 0;
    double utmp = 0, frtmp = 0;
    _f_mob = Vector3d::Zero();
    double LJ_cut = 0;
    _ubend=0;
    _uspring=0;
    // Calculate two-particle interactions between tracer and edgeparticles
    for (unsigned int i = 0; i < _N_polySpheres ; i++) {
        _polySpheres[i].f_mob  = Vector3d::Zero();
        frtmp = 0;
        vec_rij = minImage(_polySpheres[i].pos - _ppos);
        rij = vec_rij.norm();
        _Mrabs[i+1][0] = rij;
        _Mrabs[0][i+1] = rij;
        _Mrvec[i+1][0] = vec_rij; // stores vector going from tracer to _polySphere[i]
        _Mrvec[0][i+1] = -vec_rij;

        addLJPot(rij, _uLJ, frtmp, _pradius+_polyrad);

        // add total directional forces
        _f_mob += - frtmp * vec_rij;
        _polySpheres[i].f_mob += frtmp * vec_rij;
    }
    utmp = _uLJ; // store total LJ potential of tracer

    // Calculate 2p interaction for edgeParticles
    LJ_cut = 1.122462 * 2 * _polyrad;
    for (unsigned int i = 0; i < _N_polySpheres ; i++) {
        //nSpringPot=0;
        for (unsigned int j = i+1; j < _N_polySpheres ; j++) {
            frtmp = 0;
            vec_rij = minImage(_polySpheres[j].pos - _polySpheres[i].pos);
            rij = vec_rij.norm();
            _Mrabs[j+1][i+1] = rij;
            _Mrabs[i+1][j+1] = rij;
            _Mrvec[j+1][i+1] = vec_rij; // stores vector going from _polySphere[i] to _polySphere[j]
            _Mrvec[i+1][j+1] = -vec_rij;

            // LJ Interaction
            addLJPot(rij, utmp, frtmp, 2 * _polyrad);

            // Spring interaction
            if (_springMatrix[i][j] == true){
                addSpringPot(rij, _uspring, frtmp);
            }

            // add total directional forces to instances of polyspheres
            _polySpheres[i].f_mob += - frtmp * vec_rij;
            _polySpheres[j].f_mob += frtmp * vec_rij;
        }
    }
    utmp+=_uspring;

    // Loop over matrix of 3 particle bending interaction tupel _MBend[N_bendInteractions][3]
    if (_kappaBend != 0){
        Vector3d fvec1, fvec3;
        double utmp3p=0;
        int i1, i2, i3;

        for (int i=0; i<_N_springInteractions; i++){
            i1 = _MbendTupel[i][0] + 1;
            i2 = _MbendTupel[i][1] + 1;
            i3 = _MbendTupel[i][2] + 1;
            calcBendPot(_Mrvec[i1][i2],_Mrvec[i3][i2],_Mrabs[i1][i2],_Mrabs[i3][i2],utmp3p,fvec1,fvec3);
            _polySpheres[i1-1].f_mob += fvec1;
            _polySpheres[i3-1].f_mob += fvec3;
            _polySpheres[i2-1].f_mob += -(fvec1+fvec3);
            utmp+=utmp3p;
            _ubend+=utmp3p; // store bend potential
        }
    }



    _upot = utmp;
}

Vector3d CConfiguration::minImage(Vector3d rij){
    // returns disctance vector with minimal image convention.
    // For info - Check wikipedia
    double bhalf = _boxsize/2.;
    for (int i=0;i<3;i++){
        if (rij(i) > bhalf)          rij(i) -= _boxsize;
        else if (rij(i) <= -bhalf)   rij(i) += _boxsize;
    }
    return rij;
}

void CConfiguration::updateStartpos(){
    //This function is used if the particle should keep moving after each run, and not start at _resetpos again, like in the first run
    //This (hopefully) will give better averages without having to spend a lot of steps in the beginning of each run to get away from _resetpos
    for (int i = 0; i < 3; i++){
        _startpos(i) = _ppos(i) + _boxsize * _boxnumberXYZ[i];
    }
}


void CConfiguration::makeStep(){
    //move the tracer particle according to the forces and record trajectory like watched by outsider
    _prevpos = _ppos;
    _ppos += _f_mob * _timestep + _f_sto * _mu_sto;

    // move polymerSperes
    double mu_correction = _pradius/_polyrad;
    for (unsigned int i = 0; i < _polySpheres.size() ; i++) {
        _polySpheres[i].pos += _polySpheres[i].f_mob * _timestep * mu_correction + _polySpheres[i].f_sto * _mu_sto_poly;
    }
    if ((_prevpos-_ppos).squaredNorm() > 1 ){
        cout <<"\nCConfiguration.cpp ERROR: Way too big jump!!\nprevpos:\n" << _prevpos
            << "\nppos:\n" << _ppos
                << "\n_upot = " << _upot
                    << "\nubend = " << _ubend << " -- uspring = " << _uspring << " -- uLJ = " << _uLJ <<  endl;
        throw 2;
    }
    // Check if particle has crossed the confinenment of the box
    checkBoxCrossing();
}


void CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    for (int i = 0; i < 3; i++){
        if (_ppos(i) < 0){
            _ppos(i) += _boxsize;
            _boxnumberXYZ[i] -= 1;
        }
        else if (_ppos(i) > _boxsize){
            _ppos(i) -= _boxsize;
            _boxnumberXYZ[i] += 1;
        }
    }
}


void CConfiguration::calcStochasticForces(){
    // the variate generator uses m_igen (int rand number generator),
    // samples from normal distribution with variance 1 (later sqrt(2) is multiplied)
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
            *m_igen, boost::normal_distribution<double>(0, 1));

	Vector3d ran_v = Vector3d::Zero();

    ran_v(0) = ran_gen();
	ran_v(1) = ran_gen();
	ran_v(2) = ran_gen();

    _f_sto = ran_v;

    for (unsigned int i = 0; i < _polySpheres.size() ; i++) {
        _polySpheres[i].f_sto(0) = ran_gen();
        _polySpheres[i].f_sto(1) = ran_gen();
        _polySpheres[i].f_sto(2) = ran_gen();
    }
}





//***************************** MISC **************************


std::vector<double> CConfiguration::getppos(){ // returns pointer to current particle position array
	std::vector<double> pos (3);
	for (int i = 0; i < 3; i++){
		pos[i] = _ppos(i) + _boxsize * _boxnumberXYZ[i];
	}
	return pos;
}



void CConfiguration::setRanNumberGen(double seed){
    if (seed == 0) {
        m_igen = new boost::mt19937(static_cast<unsigned int>(time(NULL)));
        cout << "random seed is time!" << endl;
    } else {
        m_igen = new boost::mt19937(static_cast<unsigned int>(seed));
        cout << "random seed is " << seed << endl;
    }
}


/*****************************************************************************
 * save the coordinates to a file named _name UPDATED VERSION by Matthias
 */
void CConfiguration::saveXYZTraj(string name, const int& move, string flag) {
    Vector3d boxCoordinates;
    boxCoordinates << _boxsize *_boxnumberXYZ[0], _boxsize *_boxnumberXYZ[1], _boxsize *_boxnumberXYZ[2];
    Vector3d rtmp;
    if(flag=="w") {    //write to new file
        /*if(m_traj_file!=NULL) {
            fclose(m_traj_file);
        }*/
        m_traj_file = fopen(name.c_str(), flag.c_str());
        if(m_traj_file==NULL) {
            cout << "error creating trajfile" << endl;
        }
    }

    fprintf(m_traj_file, "%d\n%s (%8.3f %8.3f %8.3f) t=%d \n", _N_polySpheres + 1, "sim_name", _boxsize, _boxsize, _boxsize, move);


    // Tracer
    rtmp = _ppos;//+boxCoordinates;
    fprintf(m_traj_file, "%3s%9.3f%9.3f%9.3f \n","H", rtmp(0), rtmp(1),  rtmp(2));
    // polymer particles
    for (unsigned int i = 0; i < _N_polySpheres; i++) {
        rtmp = _polySpheres[i].pos;//+boxCoordinates;
        fprintf(m_traj_file, "%3s%9.3f%9.3f%9.3f \n","O", rtmp(0), rtmp(1),  rtmp(2));
    }

    //fflush(m_traj_file);

    if(flag=="c") {    //close file
        if(m_traj_file!=NULL) { fclose(m_traj_file); }
    }
}

void CConfiguration::save_traj_step(XDRFILE *xd, const int stepcount) {
    // copy necessary to convert double to float
    std::vector<std::array<float, 3> > rvecs(_N_polySpheres+1);  //TODO xtc TAKE CARE OF THIS
    // std::copy(
    //     reinterpret_cast<const double *>(&state.allPositions.front()),
    //     reinterpret_cast<const double *>(&state.allPositions.front()) + state.allPositions.size()*3,
    //     reinterpret_cast<float *>(&rvecs.front())
    // );
    //TODO xtc inefficient copying (?) check benchmark
    rvecs[0][0] = (float)_ppos(0);
    rvecs[0][1] = (float)_ppos(1);
    rvecs[0][2] = (float)_ppos(2);
    for (int i=0; i<_N_polySpheres; i++){
        for (int k=0; k <3; k++){
            rvecs[i+1][k] = (float)_polySpheres[i].pos(k);
        }
    }

    //TODO xtc   put his matrix thing in header, dont define it here everytime!
    matrix box_matrix;
    box_matrix[0][0] = _boxsize;
    box_matrix[1][1] = _boxsize;
    box_matrix[2][2] = _boxsize;
    write_xtc(xd, (_N_polySpheres+1), (stepcount), (stepcount * _timestep), box_matrix, (rvec*) &rvecs[0], 1000);  //TODO xtc redefine last parameter precision
}

void CConfiguration::saveCoordinates(ostream& trajectoryfile, unsigned int stepcount) {
    // So far this only writes the tracer particle position
	trajectoryfile << fixed << stepcount * _timestep << "\t"
        << _ppos(0)+ _boxsize *_boxnumberXYZ[0] << " " << _ppos(1) +_boxsize*_boxnumberXYZ[1]
            << " " << _ppos(2) +_boxsize*_boxnumberXYZ[2] << endl;

	// for (unsigned int i = 0; i < m_NumberParticles; i++) {
//         myfile << m_rx[i] << ", " << m_ry[i] << ", " << m_rz[i] << endl;
//     }
}