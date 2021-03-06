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
    _polyrad = modelpar.polymersize/2;   //This is needed for testOverlap for steric and HI stuff !!
	_boxsize = modelpar.boxsize;
    _n_cellsAlongb = modelpar.n_cells;
    _kappaSP = modelpar.kspring; // in kT/a^2 NO! I this would only be true, if all lengths were rescaled by the particle size. I, in contrast, rescale with 0.1b!
    _kappaBend = modelpar.kbend; // in kT
    _timestep = timestep;

    for (int i = 0; i < 3; i++){
        _ppos(i) = _boxsize/_n_cellsAlongb/2.0;
        _boxnumberXYZ[i] = 0;
        _startpos(i) = _ppos(i);
    }
    cout << "TODO IMPLEMENT OVERLAP TEST FUNCTION BETWEEN PARTICLE AND POLYMER LATTICE AT START OF SIMULATION!!!" << endl;

    initConstants(modelpar);

    // init random number generator
    setRanNumberGen(0);
    //initPolyspheres();
    //initInteractionMatrix(); // This comes only AFTER initPolyspheres
    initSemiFlexibleLattice();
    //updateLJlist();
    cout << "NOTE: Implement periodic boundary conditions for the polySpheres. I.e. checkBoxCrossing function! " << endl;
    // TEST CUE to modify the directory the output data is written to!!
    _testcue = "";
    _testcue += "/n_edge" + toString(_edgeParticles);

}

void CConfiguration::checkDisplacementforLJlist(){
    double movedsq = 0;
    for (int i=0; i<3; i++){
        movedsq += pow(_ppos[i] + _boxsize *  _boxnumberXYZ[i] - _lastcheck[i] , 2);
    }
    if ( movedsq > _cutoffUpdateLJlistSq ){
        updateLJlist();
        // new start point for displacement check
        _lastcheck[0] = _ppos[0] + _boxsize *  _boxnumberXYZ[0];
        _lastcheck[1] = _ppos[1] + _boxsize *  _boxnumberXYZ[1];
        _lastcheck[2] = _ppos[2] + _boxsize *  _boxnumberXYZ[2];
    }
}


void CConfiguration::updateLJlist(){
    Vector3d vec_rij;
    _LJlist.resize(0);
    for (unsigned int i = 0; i < _N_polySpheres ; i++) {
        vec_rij = minImage(_polySpheres[i].pos_pbc - _ppos);
        if ( vec_rij.squaredNorm() < _cutoffAddToLJlistSq ) _LJlist.push_back(i);
    }
}

void CConfiguration::calcMobilityForces(){
    Vector3d vec_rij, vec_rij_min;
    double rij = 0, rijSq=0, rij_min=0;
    double utmp = 0, frtmp = 0;
    Vector3d faddtmp;
    _f_mob = Vector3d::Zero();
    _ubend=0;
    _uspring=0;
    // Calculate two-particle interactions between tracer and edgeparticles
    // LENNARD JONES INTERACTION TRACER
    const double tracerLJSq = pow(_pradius+_polyrad, 2);
    for (unsigned int i = 0; i < _N_polySpheres ; i++) {
    //for (unsigned int l = 0; l < _LJlist.size() ; l++) {
        //const int i = _LJlist[l];
        _polySpheres[i].f_mob  = Vector3d::Zero();
        frtmp = 0;
        vec_rij = minImage(_polySpheres[i].pos_pbc - _ppos);
        rijSq = vec_rij.squaredNorm();
        if ( rijSq <  1.25992 * tracerLJSq ){ // steric parameter d_steric is the added radii of both Lennard-Jones particles
            addLJPot(rijSq, _uLJ, frtmp,tracerLJSq);

            // add total directional forces
            faddtmp = frtmp * vec_rij;
            _f_mob += - faddtmp;
            _polySpheres[i].f_mob += faddtmp;
        }
    }
    utmp = _uLJ; // store total LJ potential of tracer

    //INTERACTION POLYSPHERES
    // here; we calculate the spring and LJ interaction between the polyspheres, if the network consists of connected springs with _edgeParticles>1
    int rn;
    const double polysLJSq = pow(2.*_polyrad, 2);
    if (_edgeParticles != 1 && _kappaSP > 0){// _kappaSP > 0 since -1 means fixed spheres (no interaction) and 0 is taken care of in the else if statement below.
        for (unsigned int i = 0; i < _N_polySpheres ; i++) {
            // SPRING INTERACTION AND LENNARD JONES ONLY FOR CLOSE PARTICLES
            //TODO is LJ really included?
            //TODO is spring really included?
            for (int k=0; k<_polySpheres[i].n_rn; k++){
                frtmp = 0;
                rn = _polySpheres[i].i_rn[k]; // index of right neighbor
                vec_rij = _polySpheres[rn].pos_abs - _polySpheres[i].pos_abs;
                if (_polySpheres[i].image[k]) {
                    vec_rij += _polySpheres[i].image_corr[k];
                }
                rijSq = vec_rij.squaredNorm();  

                // SPRING
                rij = sqrt(rijSq);
                addSpringPot(rij, _uspring, frtmp);
                // LENNARD-JONES
                if ( rijSq <  1.25992 * polysLJSq )  addLJPot(rijSq, utmp, frtmp,polysLJSq);
                //if (rij > 15) cout << "Sphere " << i << " ---- rightneighbor index = " << rn << endl;
                faddtmp = frtmp * vec_rij;
                _polySpheres[i].f_mob += - faddtmp ;
                _polySpheres[rn].f_mob += faddtmp;
                if (_polySpheres[rn].n_rn == 3){
                    int j = 0;
                    //LJ Interact with these three, too!
                    for (int f=0;f<3;f++){
                        j=_polySpheres[rn].i_rn[f];
                        frtmp =0;
                        vec_rij = minImage(_polySpheres[j].pos_abs - _polySpheres[i].pos_abs);
                        rijSq = vec_rij.squaredNorm();
                        if (rijSq < 1.25992 * polysLJSq ){
                            addLJPot(rijSq, utmp, frtmp, polysLJSq);

                            faddtmp = frtmp * vec_rij;
                            _polySpheres[i].f_mob += - faddtmp;
                            _polySpheres[j].f_mob += faddtmp;
                        }
                    }

                    // BUT ALSO WITH THE OTHER NEIGHBORS OF THIS GUY!
                    // Maybe I need to store leftneighbors too in the case of a sphere at a vertex.
                    for (int f=0;f<_polySpheres[i].i_LJ.size();f++){
                        j=_polySpheres[i].i_LJ[f];
                        frtmp =0;
                        vec_rij = minImage(_polySpheres[j].pos_abs - _polySpheres[i].pos_abs);
                        rijSq = vec_rij.squaredNorm();
                        if (rijSq < 1.25992 * polysLJSq ){
                            addLJPot(rijSq, utmp, frtmp, polysLJSq);

                            faddtmp = frtmp * vec_rij;
                            _polySpheres[i].f_mob += - faddtmp;
                            _polySpheres[j].f_mob += faddtmp;
                        }
                    }
                }

            }
        }
    }
    else if(_kappaSP != -1) {//If the spheres are fixed, no interaction needs to be calculated
        // ++++++++++++++++++++++++++++++ FULL CALCULATION (slow) NEED TO USE THIS FOR _kappaSP = 0 ++++++++++++++++
        for (unsigned int i = 0; i < _N_polySpheres ; i++) {
            // LENNARD JONES -- TODO this is inefficient for edgeparticles==1
            for (unsigned int j = i+1; j < _N_polySpheres ; j++) {
                frtmp = 0;
                if (_kappaSP != 0.) vec_rij = _polySpheres[j].pos_abs - _polySpheres[i].pos_abs;
                else vec_rij = _polySpheres[j].pos_pbc - _polySpheres[i].pos_pbc;
                _Mrvec[j][i] = vec_rij; // stores vector going from _polySphere[i] to _polySphere[j]
                _Mrvec[i][j] = -vec_rij;

                vec_rij_min = minImage(vec_rij);  // For LJ interaction minImage needs to be employed
                rijSq = vec_rij_min.squaredNorm();
                if (rijSq < 1.25992 * polysLJSq ){
                    addLJPot(rijSq, utmp, frtmp, polysLJSq);

                    faddtmp = frtmp * vec_rij_min;
                    _polySpheres[i].f_mob += - faddtmp;
                    _polySpheres[j].f_mob += faddtmp;
                }
            }
            if (_kappaSP > 0.){//only calc spring pot if spheres are not fixed (ks=-1) or unconnected (ks=0)
                // SPRING INTERACTION
                for (int k=0; k<_polySpheres[i].n_rn; k++){
                    frtmp = 0;
                    rn = _polySpheres[i].i_rn[k]; // index of right neighbor
                    if (_polySpheres[i].image[k]) {
                        //cout << " **********\n" << _Mrvec[rn][i].norm() << endl;
                        _Mrvec[rn][i] += _polySpheres[i].image_corr[k];
                        _Mrvec[i][rn] = -_Mrvec[rn][i];
                        //cout << "\n--\n" << _Mrvec[rn][i].norm() << endl;
                    }
                    rij = _Mrvec[rn][i].norm();
                    _Mrabs[rn][i] = rij;
                    _Mrabs[i][rn] = rij;

                    addSpringPot(rij, _uspring, frtmp);
                    //if (rij > 15) cout << "Sphere " << i << " ---- rightneighbor index = " << rn << endl;
                    faddtmp = frtmp * _Mrvec[rn][i];
                    _polySpheres[i].f_mob += - faddtmp;
                    _polySpheres[rn].f_mob += faddtmp;
                }
            }
        }
    }

    utmp+=_uspring;

    // Loop over matrix of 3 particle bending interaction tupel _MBend[N_bendInteractions][3]
    if (_kappaBend != 0){
        Vector3d fvec1, fvec3;
        double utmp3p=0;
        int i1, i2, i3;

        for (int i=0; i<_N_springInteractions; i++){
            i1 = _MbendTupel[i][0];
            i2 = _MbendTupel[i][1];
            i3 = _MbendTupel[i][2];
            calcBendPot(_Mrvec[i1][i2],_Mrvec[i3][i2],_Mrabs[i1][i2],_Mrabs[i3][i2],utmp3p,fvec1,fvec3);
            _polySpheres[i1].f_mob += fvec1;
            _polySpheres[i3].f_mob += fvec3;
            _polySpheres[i2].f_mob += -(fvec1+fvec3);
            utmp+=utmp3p;
            _ubend+=utmp3p; // store bend potential
        }
    }



    _upot = utmp;
}


Vector3d CConfiguration::minImage(Vector3d rij){
    // returns disctance vector with minimal image convention.
    // For info - Check wikipedia
    int abc[3];
    Vector3d minvec = rij;
        for (int p=0;p<3;p++){
            abc[p]= minvec(p)*(_b2inv);
            minvec(p) -= abc[p] * _boxsize;
            abc[p]= minvec(p)*(_b2inv);
            minvec(p) -= abc[p] * _boxsize;
        }
    // Vector3d rij_rem;
//     rij_rem(0) = remainder(rij(0),_boxsize);
//     rij_rem(1) = remainder(rij(1),_boxsize);
//     rij_rem(2) = remainder(rij(2),_boxsize);
//     if (rij_rem - minvec != Vector3d::Zero()) cout << "===============\n" << rij_rem- minvec << "\nrij_rem:\n" << rij_rem <<  "\nminvec:\n" << minvec << endl;
        return minvec;
    // This was used to check, if rij_rem minimage verison produces the correct result.
    //if (  rij_rem - rij != Vector3d::Zero() )cout << "=========\n" << rij_rem << "\n--\n" << rij << "\n****\n" << endl;
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
    Vector3d dr;
    _prevpos = _ppos;
    _ppos += _f_mob * _timestep + _f_sto * _mu_sto;

    // move polymerSperes
    if (_kappaSP != -1.){
        double mu_correction = _pradius/_polyrad;
        for (unsigned int i = 0; i < _N_polySpheres ; i++) {
            dr = _polySpheres[i].f_mob * _timestep * mu_correction + _polySpheres[i].f_sto * _mu_sto_poly;
            //shift both the absolute position of the polyspheres (never shifted during sim) as well as the pbc corrected position by the same amount dr.
           _polySpheres[i].pos_abs += dr;
           _polySpheres[i].pos_pbc += dr;
        }
    }
    if ((_prevpos-_ppos).squaredNorm() > 1 ){
        cout <<"\nCConfiguration.cpp ERROR: Way too big jump!!\nprevpos:\n" << _prevpos
            << "\nppos:\n" << _ppos
                << "\n_upot = " << _upot
                    << "\nubend = " << _ubend << " -- uspring = " << _uspring << " -- uLJ = " << _uLJ <<  endl;
        throw 2;
    }
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
        for (int k=0;k<_N_polySpheres;k++){
            if (_polySpheres[k].pos_pbc(i) < 0){
                _polySpheres[k].pos_pbc(i) += _boxsize;
            }
            else if (_polySpheres[k].pos_pbc(i) > _boxsize){
                _polySpheres[k].pos_pbc(i) -= _boxsize;
            }
        }
    }
    //enable this for use of LJlist (no real speedup actually
    //checkDisplacementforLJlist();
}


void CConfiguration::calcStochasticForces(){
    // the variate generator uses m_igen (int rand number generator),
    // samples from normal distribution with variance 1 (later sqrt(2) is multiplied)
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
            *m_igen, boost::normal_distribution<double>(0, 1));


    _f_sto(0) = ran_gen();
	_f_sto(1) = ran_gen();
	_f_sto(2) = ran_gen();


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
    fprintf(m_traj_file, "%3s%9.3f%9.3f%9.3f \n","O", rtmp(0), rtmp(1),  rtmp(2));
    // polymer particles
    for (unsigned int i = 0; i < _N_polySpheres; i++) {
        rtmp = _polySpheres[i].pos_abs;//+boxCoordinates;
        fprintf(m_traj_file, "%3s%9.3f%9.3f%9.3f \n","H", rtmp(0), rtmp(1),  rtmp(2));
    }

    //fflush(m_traj_file);

    if(flag=="c") {    //close file
        if(m_traj_file!=NULL) { fclose(m_traj_file); }
    }
}

void CConfiguration::save_traj_step(XDRFILE *xd, const int stepcount) {
    cout << "Warning! save_traj_step: So far saves the position of the particles relative to simulation box! Need to add _boxsize." << endl;
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
            rvecs[i+1][k] = (float)_polySpheres[i].pos_pbc(k);
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
