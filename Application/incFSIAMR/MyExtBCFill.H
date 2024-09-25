#ifndef _incFSI_H_
#define _incFSI_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Mon Dec 21 16:40:52 2020


************************************************************************/

#include <amrexMesh.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LO_BCTYPES.H>
#include <functional>
#include <defs.H>

#include <Interface.H>
#include <Cylinder.H>
#include <Ellipse.H>
#include <Foil.H>
#include <Bubble.H>
#include <InterceptData.H>

#include <Mask.H>

#include <AMReX_AmrCore.H>
#include <AMReX_PhysBCFunct.H>

namespace mycode
{

    struct MyExtBCFill {
        int velflag;// 1 = u, 2 = v, 3 = w
        amrex::BoxArray ba;
        BCArray lobc;
        BCArray hibc;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> bcf;
        bool isAxisymmetric;
        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                         const int dcomp, const int numcomp,
                         amrex::GeometryData const& geom, const amrex::Real time,
                         const amrex::BCRec* bcr, const int bcomp,
                         const int orig_comp) const
        {
            using namespace amrex;
    
            // do something for external Dirichlet (BCType::ext_dir)
            const int i = iv[0];
            const int j = iv[1];
            const int k = 0;
            const amrex::Box &domain = geom.Domain();
            const amrex::Real *prob_lo = geom.ProbLo();
            const amrex::Real* prob_hi = geom.ProbHi();
            const amrex::Real* dx = geom.CellSize();
    
            if(velflag == 1)//Boundaru condition for u velocity
            {
                amrex::Real x = prob_lo[0] + dx[0] * i;
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                if(i == domain.smallEnd(0))//Left boundary
                {
                    //int shift = i - domain.smallEnd(0);
                    //if(shift == 0)
                    {
                        if (lobc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = bcf[0](x, y, time);
                        else if (lobc[0] == amrex::LinOpBCType::Neumann)
                        {
                            //dest(i, j, k) = dest(i + 1, j, k) - dx[0] * bcf[0](x, y, time);
                            dest(i - 1, j, k) = dest(i + 1, j, k) - 2.0*dx[0] * bcf[0](x, y, time);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.smallEnd(1) - 1)//bottom boundary
                {
                    int shift = j - domain.smallEnd(1);
                    if(shift == -1)
                    {
                        if (lobc[1] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[1](x, y, time) - dest(i , j + 1, k);
                        else if (lobc[1] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = -1.0* dx[1] * bcf[1](x, y, time) + dest(i , j + 1, k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(i == domain.bigEnd(0) + 1)//right boundary
                {
                    //int shift = i - (domain.bigEnd(0) + 1);
                    //if(shift == 0)
                    {
                        if (hibc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = bcf[2](x, y, time) ;
                        else if (hibc[0] == amrex::LinOpBCType::Neumann)
                        { 
                            //dest(i, j, k) = dx[0] * bcf[2](x, y, time) + dest(i - 1, j, k); 
                            dest(i + 1, j, k) = 2.0 * dx[0] * bcf[2](x, y, time) + dest(i - 1, j, k);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.bigEnd(1) + 1)//top boundary
                {
                    //int shift = j - domain.bigEnd(1);
                    //if(shift == 1)
                    {
                        if (hibc[1] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[3](x, y, time) + dest(i , j - 1 ,k);
                        else if (hibc[1] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = dx[1] * bcf[3](x, y, time) + dest(i , j - 1 ,k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else
                    dest(i, j, k) = 0.0;
               //dest(i, j, k) = 0.0;
            }
            if(velflag == 2)//Boundary condition for V velocity
            {   
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5); 
                amrex::Real y = prob_lo[1] + dx[1] * j;
                if(i == domain.smallEnd(0) - 1)//Left boundary
                {
                    //int shift = i - domain.smallEnd(0);
                    //if(shift == -1)
                    {
                        if (lobc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[0](x, y, time) - dest(i + 1, j, k);
                        else if (lobc[0] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = -1.0* dx[0] * bcf[0](x, y, time) + dest(i + 1, j, k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.smallEnd(1))//bottom boundary
                {
                    //int shift = j - domain.smallEnd(1);
                    //if(shift == 0)
                    {
                        if (lobc[1] == amrex::LinOpBCType::Dirichlet) 
                            dest(i , j, k) = bcf[1](x, y, time);
                        else if (lobc[1] == amrex::LinOpBCType::Neumann)
                        {   
                            if(isAxisymmetric)
                            {
                                dest(i, j, k) = 0.0;
                            }
                            //else
                            //{
                            //    dest(i , j, k) =  -1.0 * dx[1] * bcf[1](x, y, time) + dest(i, j + 1, k);
                            //}
                            dest(i , j - 1, k) = -2.0 * dx[1] * bcf[1](x, y, time) + dest(i, j + 1, k);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(i == domain.bigEnd(0) + 1)//right boundary
                {
                    //int shift = i - domain.bigEnd(0) ;
                    //if(shift == 1)
                    {
                        if (hibc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[2](x, y, time) - dest(i - 1, j, k);
                        else if (hibc[0] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = dx[0] * bcf[2](x, y, time) + dest(i - 1, j, k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.bigEnd(1) + 1)//top boundary
                {
                    //int shift = j - (domain.bigEnd(1) + 1);
                    //if(shift == 0)
                    {
                        if (hibc[1] == amrex::LinOpBCType::Dirichlet) 
                            dest(i , j, k) = bcf[3](x, y, time);
                        else if (hibc[1] == amrex::LinOpBCType::Neumann)
                        {
                            //dest(i , j, k) = dx[1] * bcf[3](x, y, time) + dest(i, j - 1, k);
                            dest(i , j + 1, k) = 2.0 * dx[1] * bcf[3](x, y, time) + dest(i, j - 1, k);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else
                    dest(i, j, k) = 0.0;
            }
            //dest(i, j, k) = 0.0;
        }
    };

} /*End namespace mycode */

#endif
