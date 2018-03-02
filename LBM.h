/* This code accompanies
 *   The Lattice Boltzmann Method: Principles and Practice
 *   T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
 *   ISBN 978-3-319-44649-3 (Electronic) 
 *        978-3-319-44647-9 (Print)
 *   http://www.springer.com/978-3-319-44647-9
 *
 * This code is provided under the MIT license. See LICENSE.txt.
 *
 * Author: Orest Shardt
 *
 */
#ifndef __LBM_H
#define __LBM_H

const unsigned int scale = 1;
const unsigned int NX = 200;
const unsigned int NY = 64;

const unsigned int ndir = 9;
const size_t mem_size_0dir   = sizeof(double)*NX*NY;
const size_t mem_size_n0dir  = sizeof(double)*NX*NY*(ndir-1);
const size_t mem_size_scalar = sizeof(double)*NX*NY;
const size_t mem_size_bool = sizeof(bool)*NX*NY;

const double w0 = 4.0/9.0;  // zero weight
const double ws = 1.0/9.0;  // adjacent weight
const double wd = 1.0/36.0; // diagonal weight



// Simulation Parameters
const double u_max = 0.05;
const double rho0 = 1.0;
const double CS2 = 1.0/3.0; 
const double nu = 0.01; //viscosity
const double Lambda_nu = 1.0/(nu/CS2+0.5);  // Lamba_nu

const unsigned int NSTEPS = 200000;
const unsigned int NSAVE  =  50000;
const unsigned int NMSG   =  10000;

// compute L2 error and energy?
// disable for speed testing
const bool computeFlowProperties = true;

// suppress verbose output
const bool quiet = true;

void taylor_green(unsigned int,unsigned int,unsigned int,double*,double*,double*);
void taylor_green(unsigned int,double*,double*,double*);
void stream_collide_save(double*,double*,double*,double*,double*,double*,bool);
void init_equilibrium(double*,double*,double*,double*,double*);
void compute_flow_properties_taylor(unsigned int,double*,double*,double*,double*);
void report_flow_properties_taylor(unsigned int,double*,double*,double*);
void save_scalar(const char*,double*,unsigned int);
void save_vector(const char*, double *,double *, unsigned int );
void init_obs(bool*);
void init_rest(double*, double*, double*, bool*);
void obs_collideBGK_stream_save(double*,double*,double*,double*,double*,double*,bool, bool*);
void compute_flow_properties(unsigned int,double*,double*,double*,double*, bool*);
void report_flow_properties(unsigned int,double*,double*,double*,bool*);
void compute_curl(double *, double *, double *, bool *);
void obs_collideMRT_stream_save(double*,double*,double*,double*,double*,double*,bool, bool*, double*);


void save_obs(bool *);

inline size_t field0_index(unsigned int x, unsigned int y)
{
    return NX*y+x;
}

inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return NX*y+x;
}

inline size_t fieldn_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (ndir-1)*(NX*y+x)+(d-1);
}

#endif /* __LBM_H */

