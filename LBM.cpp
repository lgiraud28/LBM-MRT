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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "LBM.h"

void taylor_green(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v)
{
    double kx = 2.0*M_PI/NX;
    double ky = 2.0*M_PI/NY;
    double td = 1.0/(nu*(kx*kx+ky*ky));
    
    double X = x+0.5;
    double Y = y+0.5;
    double ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
    double uy =  u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);
    double P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2.0*kx*X)+(kx/ky)*cos(2.0*ky*Y))*exp(-2.0*t/td);
    double rho = rho0+3.0*P;
    
    *r = rho;
    *u = ux;
    *v = uy;
}

void taylor_green(unsigned int t, double *r, double *u, double *v)
{
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            size_t sidx = scalar_index(x,y);

            taylor_green(t,x,y,&r[sidx],&u[sidx],&v[sidx]);
        }
    }
}

void init_rest(double *r, double *u, double *v, bool *obs)
{
	
    double ux = 0.1;
    double uy =  0.0;

    /* initialize random seed: */
    srand (time(NULL));
    
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            size_t sidx = scalar_index(x,y);
            if (!obs[sidx]) 
            {
           		 r[sidx]=rho0;
           		 u[sidx]=ux+(rand() %100)/1000000.0;  //small perturbation
            	 v[sidx]=uy+(rand() %100)/1000000.0;
            } else
            {
            	r[sidx]=rho0;
            	u[sidx]=0.0;
           		v[sidx]=0.0;
            }
        }
    }
    
}

void init_equilibrium( double *f0, double *f1,double *r, double *u, double *v)
{
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x,y)];
            double ux  = u[scalar_index(x,y)];
            double uy  = v[scalar_index(x,y)];
            
            // load equilibrium
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]
            
            // temporary variables
            double w0r = w0*rho;
            double wsr = ws*rho;
            double wdr = wd*rho;
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy);
            
            double tux = 3.0*ux;
            double tuy = 3.0*uy;
            
            f0[field0_index(x,y)]    = w0r*(omusq);
            
            double cidot3u = tux;
            f1[fieldn_index(x,y,1)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f1[fieldn_index(x,y,2)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f1[fieldn_index(x,y,3)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f1[fieldn_index(x,y,4)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            
            cidot3u = tux+tuy;
            f1[fieldn_index(x,y,5)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f1[fieldn_index(x,y,6)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f1[fieldn_index(x,y,7)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f1[fieldn_index(x,y,8)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}

void stream_collide_save(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save)
{
    // useful constants
    const double tauinv = 2.0/(6.0*nu+1.0); // 1/Lamba_nu
    const double omtauinv = 1.0-tauinv;     // 1 - 1/Lamba_nu

    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x+1)%NX;
            unsigned int yp1 = (y+1)%NY;
            unsigned int xm1 = (NX+x-1)%NX;
            unsigned int ym1 = (NY+y-1)%NY;
            
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8
            
            double ft0 = f0[field0_index(x,y)];
            
            // load populations from adjacent nodes
            double ft1 = f1[fieldn_index(xm1,y,  1)];
            double ft2 = f1[fieldn_index(x,  ym1,2)];
            double ft3 = f1[fieldn_index(xp1,y,  3)];
            double ft4 = f1[fieldn_index(x,  yp1,4)];
            double ft5 = f1[fieldn_index(xm1,ym1,5)];
            double ft6 = f1[fieldn_index(xp1,ym1,6)];
            double ft7 = f1[fieldn_index(xp1,yp1,7)];
            double ft8 = f1[fieldn_index(xm1,yp1,8)];
            
            // compute moments
            double rho = ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8;
            double rhoinv = 1.0/rho;
            
            double ux = rhoinv*(ft1+ft5+ft8-(ft3+ft6+ft7));
            double uy = rhoinv*(ft2+ft5+ft6-(ft4+ft7+ft8));
            
            // only write to memory when needed
            if(save)
            {
                r[scalar_index(x,y)] = rho;
                u[scalar_index(x,y)] = ux;
                v[scalar_index(x,y)] = uy;
            }
            
            // now compute and relax to equilibrium
            // note that
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]
            
            // temporary variables
            double tw0r = tauinv*w0*rho; //   w[0]*rho/tau 
            double twsr = tauinv*ws*rho; // w[1-4]*rho/tau
            double twdr = tauinv*wd*rho; // w[5-8]*rho/tau
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy); // 1-(3/2)u.u
            
            double tux = 3.0*ux;
            double tuy = 3.0*uy;
            
            
            f0[field0_index(x,y)]    = omtauinv*ft0  + tw0r*(omusq);
            
            double cidot3u = tux;
            f2[fieldn_index(x,y,1)]  = omtauinv*ft1  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f2[fieldn_index(x,y,2)]  = omtauinv*ft2  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f2[fieldn_index(x,y,3)]  = omtauinv*ft3  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f2[fieldn_index(x,y,4)]  = omtauinv*ft4  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            
            cidot3u = tux+tuy;
            f2[fieldn_index(x,y,5)]  = omtauinv*ft5  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f2[fieldn_index(x,y,6)]  = omtauinv*ft6  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f2[fieldn_index(x,y,7)]  = omtauinv*ft7  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f2[fieldn_index(x,y,8)]  = omtauinv*ft8  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}

void obs_collideBGK_stream_save(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save,bool* obs)
{
    // useful constants
    const double tauinv = 2.0/(6.0*nu+1.0); // 1/tau
    const double omtauinv = 1.0-tauinv;     // 1 - 1/tau
	double rho, ux, uy;


    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x+1)%NX;
            unsigned int yp1 = (y+1)%NY;
            unsigned int xm1 = (NX+x-1)%NX;
            unsigned int ym1 = (NY+y-1)%NY; 
            
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8
            
            double ft0 = f0[field0_index(x,y)];
            // load populations from adjacent nodes
            double ft1 = f1[fieldn_index(x,  y,1)];
            double ft2 = f1[fieldn_index(x,  y,2)];
            double ft3 = f1[fieldn_index(x,  y,3)];
            double ft4 = f1[fieldn_index(x,  y,4)];
            double ft5 = f1[fieldn_index(x,  y,5)];
            double ft6 = f1[fieldn_index(x,  y,6)];
            double ft7 = f1[fieldn_index(x,  y,7)];
            double ft8 = f1[fieldn_index(x,  y,8)];

            if (!obs[scalar_index(x,y)])
            {
            	// compute moments
            	rho = ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8;
            	double rhoinv = 1.0/rho;
            
            	ux = rhoinv*(ft1+ft5+ft8-(ft3+ft6+ft7));
            	// forcing term
            	ux=ux+0.0001;
            	uy = rhoinv*(ft2+ft5+ft6-(ft4+ft7+ft8));
            // now compute and relax to equilibrium
            // note that
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]
            
            // temporary variables
            	double tw0r = tauinv*w0*rho; //   w[0]*rho/tau 
            	double twsr = tauinv*ws*rho; // w[1-4]*rho/tau
            	double twdr = tauinv*wd*rho; // w[5-8]*rho/tau
            	double omusq = 1.0 - 1.5*(ux*ux+uy*uy); // 1-(3/2)u.u
            
            	double tux = 3.0*ux;
            	double tuy = 3.0*uy;
            
            
            	f0[field0_index(x,y)]    = omtauinv*ft0  + tw0r*(omusq);
            
            	double cidot3u = tux;
            	f2[fieldn_index(xp1,y,1)]  = omtauinv*ft1  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            	cidot3u = tuy;
            	f2[fieldn_index(x,yp1,2)]  = omtauinv*ft2  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            	cidot3u = -tux;
            	f2[fieldn_index(xm1,y,3)]  = omtauinv*ft3  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            	cidot3u = -tuy;
            	f2[fieldn_index(x,ym1,4)]  = omtauinv*ft4  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            
            	cidot3u = tux+tuy;
            	f2[fieldn_index(xp1,yp1,5)]  = omtauinv*ft5  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            	cidot3u = tuy-tux;
            	f2[fieldn_index(xm1,yp1,6)]  = omtauinv*ft6  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            	cidot3u = -(tux+tuy);
            	f2[fieldn_index(xm1,ym1,7)]  = omtauinv*ft7  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            	cidot3u = tux-tuy;
            	f2[fieldn_index(xp1,ym1,8)]  = omtauinv*ft8  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            } else
            {    //this is an obstable bounce back
            	f2[fieldn_index(xp1,y,  1)] = ft3;
             	f2[fieldn_index(x,yp1,  2)] = ft4;
            	f2[fieldn_index(xm1,y,  3)] = ft1;
            	f2[fieldn_index(x,ym1,  4)] = ft2;
            	f2[fieldn_index(xp1,yp1,  5)] = ft7;
            	f2[fieldn_index(xm1,yp1,  6)] = ft8;
            	f2[fieldn_index(xm1,ym1,  7)] = ft5;
            	f2[fieldn_index(xp1,ym1,  8)] = ft6;

           		rho=0.0;
           		ux=0.0;
           		uy=0.0;
            }
            
            
            // only write to memory when needed
            if(save)
            {
                r[scalar_index(x,y)] = rho;
                u[scalar_index(x,y)] = ux;
                v[scalar_index(x,y)] = uy;
            }
            
            
        }
    }
}
  
void compute_flow_properties_taylor(unsigned int t, double *r, double *u, double *v, double *prop)
{
    // prop must point to space for 4 doubles:
    // 0: energy
    // 1: L2 error in rho
    // 2: L2 error in ux
    // 3: L2 error in uy
    
    double E = 0.0; // kinetic energy
    
    double sumrhoe2 = 0.0; // sum of error squared in rho
    double sumuxe2 = 0.0;  //                         ux
    double sumuye2 = 0.0;  //                         uy
    
    double sumrhoa2 = 0.0; // sum of analytical rho squared
    double sumuxa2 = 0.0;  //                   ux
    double sumuya2 = 0.0;  //                   uy
    
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x,y)];
            double ux  = u[scalar_index(x,y)];
            double uy  = v[scalar_index(x,y)];
            E += rho*(ux*ux + uy*uy);
            
            double rhoa, uxa, uya;
            taylor_green(t,x,y,&rhoa,&uxa,&uya);
            
            sumrhoe2 += (rho-rhoa)*(rho-rhoa);
            sumuxe2  += (ux-uxa)*(ux-uxa);
            sumuye2  += (uy-uya)*(uy-uya);

            sumrhoa2 += (rhoa-rho0)*(rhoa-rho0);
            sumuxa2  += uxa*uxa;
            sumuya2  += uya*uya;
        }
    }
    
    prop[0] = E;
    prop[1] = sqrt(sumrhoe2/sumrhoa2);
    prop[2] = sqrt(sumuxe2/sumuxa2);
    prop[3] = sqrt(sumuye2/sumuya2);
}
void compute_flow_properties(unsigned int t, double *r, double *u, double *v, double *prop, bool *obs)
{
    
    
    int i=0;
    double sumrho = 0.0; // sum of error squared in rho
    double sumux = 0.0;
    double sumuy = 0.0;

    
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            if (!obs[scalar_index(x,y)])
            	{ 
            		i++;
            		sumrho += r[scalar_index(x,y)]-rho0;
            		sumux+=u[scalar_index(x,y)];
            		sumuy+=v[scalar_index(x,y)];
            	}
            
            
        }
    }
    
    prop[0] = sumrho;
    prop[1] = sumux/i;
    prop[2] = sumuy/i;
    prop[3] = (double)i;
}

void report_flow_properties(unsigned int t, double *rho, double *ux, double *uy, bool *obs)
{
    double prop[4];
    compute_flow_properties(t,rho,ux,uy,prop, obs);
    printf("time step %u, mass: %g av Ux:%g av Uy:%g active nodes: %g  Obs nodes: %g\n",t,prop[0],prop[1],prop[2],prop[3], NX*NY-prop[3]);
}
void report_flow_properties_taylor(unsigned int t, double *rho, double *ux, double *uy)
{
    double prop[4];
    compute_flow_properties_taylor(t,rho,ux,uy,prop);
    printf("%u,%g,%g,%g,%g\n",t,prop[0],prop[1],prop[2],prop[3]);
}

void save_scalar(const char* name, double *scalar, unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[16];
    unsigned int i,j;
    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS)+1.0);
    
    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format,"%%s%%0%dd.vtk",ndigits);
    sprintf(filename,format,name,n);
    
    // open file for writing
    FILE *fout = fopen(filename,"w+");
    
    // write data
    //fwrite(scalar,1,mem_size_scalar,fout);
    fprintf(fout,"# vtk DataFile Version 2.0 \n");
    fprintf(fout,"Lattice_Boltzmann_fluid_flow\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET RECTILINEAR_GRID\n");
    fprintf(fout,"DIMENSIONS %d %d 1\n",NX,NY); 
    fprintf(fout,"X_COORDINATES %d float\n",NX); 
    for (i=0;i<NX;i++)
        fprintf(fout,"%d ", i);
    fprintf(fout," \n");
    fprintf(fout,"Y_COORDINATES %d float\n",NY);
    for (i=0;i<NY;i++)
        fprintf(fout,"%d ", i);
    fprintf(fout," \n");
    fprintf(fout,"Z_COORDINATES 1 float\n");
    fprintf(fout,"0\n");
    fprintf(fout,"POINT_DATA %d\n",NX*NY);
    fprintf(fout," \n");
    //fprintf(fout," \n");
    fprintf(fout,"SCALARS %s float 1 \n",name);
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i=0;i<NY;i++) 
    {   for (j=0;j<NX;j++)
           fprintf(fout,"%f ", scalar[scalar_index(j,i)]);
        fprintf(fout," \n");
    }
    


    
    // close file
    fclose(fout);
    
    if(ferror(fout))
    {
        fprintf(stderr,"Error saving to %s\n",filename);
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to %s\n",filename);
    }
}

void save_vector(const char* name, double *ux,double *uy, unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[16];
    unsigned int i,j;
    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS)+1.0);
    
    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format,"%%s%%0%dd.vtk",ndigits);
    sprintf(filename,format,name,n);
    
    // open file for writing
    FILE *fout = fopen(filename,"w+");
    
    // write data
    //fwrite(scalar,1,mem_size_scalar,fout);
    fprintf(fout,"# vtk DataFile Version 2.0 \n");
    fprintf(fout,"Lattice_Boltzmann_fluid_flow\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET RECTILINEAR_GRID\n");
    fprintf(fout,"DIMENSIONS %d %d 1\n",NX,NY); 
    fprintf(fout,"X_COORDINATES %d float\n",NX); 
    for (i=0;i<NX;i++)
        fprintf(fout,"%d ", i);
    fprintf(fout," \n");
    fprintf(fout,"Y_COORDINATES %d float\n",NY);
    for (i=0;i<NY;i++)
        fprintf(fout,"%d ", i);
    fprintf(fout," \n");
    fprintf(fout,"Z_COORDINATES 1 float\n");
    fprintf(fout,"0\n");
    fprintf(fout,"POINT_DATA %d\n",NX*NY);
    fprintf(fout," \n");
    fprintf(fout," \n");
    fprintf(fout,"VECTORS %s float \n",name);
    for (i=0;i<NY;i++) 
    {   for (j=0;j<NX;j++)
           fprintf(fout,"%f %f 0.00000 ", ux[scalar_index(j,i)], uy[scalar_index(j,i)]);
        fprintf(fout," \n");
    }
    


    
    // close file
    fclose(fout);
    
    if(ferror(fout))
    {
        fprintf(stderr,"Error saving to %s\n",filename);
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to %s\n",filename);
    }
}
void init_obs(bool* obs )
{
	unsigned int y,x;

   for(y = 0; y < NY; ++y)
    {
        for(x = 0; x < NX; ++x)
        {
            obs[scalar_index(x,y)]=false;
        }
    }
    // init a plane obstacle
    x=32;
    for (y=NY/2-8; y<NY/2+8;y++)
    	obs[scalar_index(32,y)]=true; 
}
void save_obs(bool *obs)
{
	// assume reasonably-sized file names
    
    unsigned int i,j;
    float value;
    
       
    // open file for writing
    FILE *fout = fopen("obs.vtk","w+");
    
    // write data
    //fwrite(scalar,1,mem_size_scalar,fout);
    fprintf(fout,"# vtk DataFile Version 2.0 \n");
    fprintf(fout,"Lattice_Boltzmann_fluid_flow : Obstacles\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET RECTILINEAR_GRID\n");
    fprintf(fout,"DIMENSIONS %d %d 1\n",NX,NY); 
    fprintf(fout,"X_COORDINATES %d float\n",NX); 
    for (i=0;i<NX;i++)
        fprintf(fout,"%d ", i);
    fprintf(fout," \n");
    fprintf(fout,"Y_COORDINATES %d float\n",NY);
    for (i=0;i<NY;i++)
        fprintf(fout,"%d ", i);
    fprintf(fout," \n");
    fprintf(fout,"Z_COORDINATES 1 float\n");
    fprintf(fout,"0\n");
    fprintf(fout,"POINT_DATA %d\n",NX*NY);
    fprintf(fout," \n");
    //fprintf(fout," \n");
    fprintf(fout,"SCALARS Obstacles float 1 \n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (j=0;j<NY;j++) 
    {   for (i=0;i<NX;i++)
    	{
    	   if (obs[scalar_index(i,j)]==true) value=1.0; else value=0.0;
           fprintf(fout,"%f ", value);
        }
        fprintf(fout," \n");
    }
    


    
    // close file
    fclose(fout);
    
    if(ferror(fout))
    {
        fprintf(stderr,"Error saving to obs.vtk\n");
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to obs.vtk\n");
    }
}
void compute_curl(double *u, double *v, double *c,bool *obs)
{
	for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            if (!obs[scalar_index(x,y)])
            {
               unsigned int xp1 = (x+1)%NX;
               unsigned int yp1 = (y+1)%NY;
               unsigned int xm1 = (NX+x-1)%NX;
               unsigned int ym1 = (NY+y-1)%NY;

               c[scalar_index(x,y)]=v[scalar_index(xm1,y)]-v[scalar_index(xp1,y)]-u[scalar_index(x,ym1)]+u[scalar_index(x,yp1)];
            } else
             c[scalar_index(x,y)]=0.0;
        }
    }

}
void obs_collideMRT_stream_save(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save,bool* obs, double* Phy)
{
   // useful constants
    
	double rho, ux, uy;


    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x+1)%NX;
            unsigned int yp1 = (y+1)%NY;
            unsigned int xm1 = (NX+x-1)%NX;
            unsigned int ym1 = (NY+y-1)%NY; 
            
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8
            
            double ft0 = f0[field0_index(x,y)];
            // load populations from adjacent nodes
            double ft1 = f1[fieldn_index(x,  y,1)];
            double ft2 = f1[fieldn_index(x,  y,2)];
            double ft3 = f1[fieldn_index(x,  y,3)];
            double ft4 = f1[fieldn_index(x,  y,4)];
            double ft5 = f1[fieldn_index(x,  y,5)];
            double ft6 = f1[fieldn_index(x,  y,6)];
            double ft7 = f1[fieldn_index(x,  y,7)];
            double ft8 = f1[fieldn_index(x,  y,8)];

            if (!obs[scalar_index(x,y)])
            {
            	// compute moments
            	double aa= ft1+ft2+ft3+ft4;
            	double bb=ft5+ft6+ft7+ft8;
            	rho = ft0+aa+bb;
            	double rhoinv=1.0/rho;
            
            	ux = (ft1+ft5+ft8-(ft3+ft6+ft7))*rhoinv;
            	// forcing term
            	ux=ux+0.005;
            	uy = (ft2+ft5+ft6-(ft4+ft7+ft8))*rhoinv;
            	
            	double cc=ft5-ft6;
            	double dd=ft7-ft8;
            	double ux2 = ux * ux;				// pre-compute terms used repeatedly...
				double uy2 = uy * uy;
				double u2 = ux2 + uy2;
				double ss=4.0*ft0+Phy[4]*(u2);
				double E=ss+ (aa) -2.0* (bb)-Phy[2]*rho;
				double H=ss-2.0* (aa) +(bb)-Phy[3]*rho;
				double XX=2.0*(ft1-ft3)-cc+ dd -ux;
				double XY=2.0*(ft2-ft4)-(ft5+ft6)+ft7+ft8-uy;
				double SX=ft1-ft2+ft3-ft4-Phy[5]*(ux2-uy2);
				double SY=cc+dd-Phy[5]*ux*uy;
				double ee= Phy[8]*E;
				double ff= Phy[12]*E;
				double ii= Phy[7]*H;
				double jj= Phy[9]*H;
				double kk= Phy[14]*XX;
				double ll= Phy[14]*XY;
				double mm= Phy[11]*SX;
				double nn= Phy[11]*SY;
				double oo= Phy[10]*(XX+XY);
				double pp= Phy[10]*(XY-XX);
				double qq= -ee+ii;
				double rr=ff-jj;
            
            	f0[field0_index(x,y)]    = ft0 -Phy[13]*E-Phy[6]*H;
            
            
            	f2[fieldn_index(xp1,y,1)]  = ft1+qq-kk-mm;
            	f2[fieldn_index(x,yp1,2)]  = ft2+qq-ll+mm;
            	f2[fieldn_index(xm1,y,3)]  = ft3+qq+kk-mm;
            	f2[fieldn_index(x,ym1,4)]  = ft4+qq+ll+mm;
            
            	f2[fieldn_index(xp1,yp1,5)]  = ft5+rr+oo-nn;
            	f2[fieldn_index(xm1,yp1,6)]  = ft6+rr+pp+nn;
            	f2[fieldn_index(xm1,ym1,7)]  = ft7+rr-oo-nn;
            	f2[fieldn_index(xp1,ym1,8)]  = ft8+rr-pp+nn;
            } else
            {    //this is an obstable bounce back
            	f2[fieldn_index(xp1,y,  1)] = ft3;
             	f2[fieldn_index(x,yp1,  2)] = ft4;
            	f2[fieldn_index(xm1,y,  3)] = ft1;
            	f2[fieldn_index(x,ym1,  4)] = ft2;
            	f2[fieldn_index(xp1,yp1,  5)] = ft7;
            	f2[fieldn_index(xm1,yp1,  6)] = ft8;
            	f2[fieldn_index(xm1,ym1,  7)] = ft5;
            	f2[fieldn_index(xp1,ym1,  8)] = ft6;

           		rho=ft0+ft1+ft2+ft2+ft3*ft4+ft5+ft6+ft7+ft8;
           		ux=0.0;
           		uy=0.0;
            }
            
            
            // only write to memory when needed
            if(save)
            {
                r[scalar_index(x,y)] = rho;
                u[scalar_index(x,y)] = ux;
                v[scalar_index(x,y)] = uy;
            }
            
            
        }
    }
}
