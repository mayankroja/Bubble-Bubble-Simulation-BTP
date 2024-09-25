#include "incFSI.H"

namespace mycode
{
    
void incFSI::XVelBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf)
{
    /// fill internal ghost cells
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &u = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (lobc_u[0] == amrex::LinOpBCType::Dirichlet)
                {
                    u(i, j, k) = u_bcf[0](x, y, time);
                }
                else if (lobc_u[0] == amrex::LinOpBCType::Neumann)
                {
                    //u(i, j, k) = u(i + 1, j, k) - dx[0] * u_bcf[0](x, y, time);
                    u(i - 1, j, k) = u(i + 1, j, k) - 2.0 * dx[0] * u_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for u on xlo face");
                }
                /*if( i == 0 && j == 10 && lev == 1)
                {
                    amrex::Print()<<"u Boundary condition "<<'\n';
                    amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
                    amrex::Print()<<"u("<<i<<","<< j<<", "<<k<<") = "<<u(i, j, k)<<'\n';
                    amrex::Print()<<"u("<<i - 1<<", "<<j <<", "<<k<<") = "<<u(i - 1, j, k)<<'\n';
                    amrex::Print()<<"*****************************"<<'\n';
                }*/

            }
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            //if(isAxisymmetric)//Assuming axis of symmetry is along y = 0
            //if ( isAxisymmetric && lobc_u[1] == amrex::LinOpBCType::Neumann)
            //{
            //    for (int shift = 1; shift <= Nghost; shift++)
            //    {
            //        for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            //            u(i, j - shift, k) = u(i, j + shift - 1, k);
            //    }
            //}
            //else
            {
                for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                {
                    x = prob_lo[0] + dx[0] * i;
                    if (lobc_u[1] == amrex::LinOpBCType::Dirichlet)
                    {
                        u(i, j - 1, k) = 2.0 * u_bcf[1](x, y, time) - u(i, j, k);
                    }
                    else if (lobc_u[1] == amrex::LinOpBCType::Neumann)
                    {
                        u(i, j - 1, k) = u(i, j, k) - dx[1] * u_bcf[1](x, y, time);
                    }
                    else
                    {
                        amrex::Abort("invalid bc for u on ylo face");
                    }                
                }
            }
        }
        
        /// right
        if (bx.bigEnd(0) > domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (hibc_u[0] == amrex::LinOpBCType::Dirichlet)
                {
                    u(i, j, k) = u_bcf[2](x, y, time);
                }
                else if (hibc_u[0] == amrex::LinOpBCType::Neumann)
                {
                    //u(i, j, k) = u(i - 1, j, k) + dx[0] * u_bcf[2](x, y, time);
                    u(i + 1, j, k) = u(i - 1, j, k) + 2.0 * dx[0] * u_bcf[2](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for u on xhi face");
                }
            }
        }
 
        /// top
        //std::cout<<"bx = "<<bx<<'\n';
        //std::cout<<"domain = "<<domain<<'\n';
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            {
                x = prob_lo[0] + dx[0] * i;

               //if((x < 0.0 && x > -0.05) || (x > 0.0 && x <0.05))
               //{
               //     std::cout<<"Location :"<<x<<" , "<<y<<'\n';
               //     std::cout<<"hibc_u = "<<hibc_u<<'\n';
               //}

                if (hibc_u[1] == amrex::LinOpBCType::Dirichlet)
                {
                    //if((x < 0.0 && x > -0.05) || (x > 0.0 && x <0.05))
                    //{
                    //    amrex::Print()<<"in BC Location :"<<x<<" , "<<y<<'\n';
                    //}

                    u(i, j + 1, k) = 2.0 * u_bcf[3](x, y, time) - u(i, j, k);
                }
                else if (hibc_u[1] == amrex::LinOpBCType::Neumann)
                {
                    u(i, j + 1, k) = u(i, j, k) + dx[1] * u_bcf[3](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for u on yhi face");
                }                
                /*if( i == 350 && j > 500)
                {
                    amrex::Print()<<"u Boundary condition "<<'\n';
                    amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
                    amrex::Print()<<"u("<<i<<","<< j<<", "<<k<<") = "<<u(i, j, k)<<'\n';
                    amrex::Print()<<"u("<<i<<", "<<j + 1<<", "<<k<<") = "<<u(i, j + 1, k)<<'\n';
                    amrex::Print()<<"*****************************"<<'\n';
                }*/

            }
        }
    }
}

void incFSI::YVelBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf)
{
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &v = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * j;
                if (lobc_v[0] == amrex::LinOpBCType::Dirichlet)
                {
                    v(i - 1, j, k) = 2.0 * v_bcf[0](x, y, time) - v(i, j, k);
                }
                else if (lobc_v[0] == amrex::LinOpBCType::Neumann)
                {
                    v(i - 1, j, k) = v(i, j, k) - dx[0] * v_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for v on xlo face");
                }
                /*if( i == 0 && j == 10 && lev == 1)
                {
                    amrex::Print()<<"v Boundary condition "<<'\n';
                    amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
                    amrex::Print()<<"v("<<i<<","<< j<<", "<<k<<") = "<<v(i, j, k)<<'\n';
                    amrex::Print()<<"v("<<i - 1<<", "<<j <<", "<<k<<") = "<<v(i - 1, j, k)<<'\n';
                    amrex::Print()<<"*****************************"<<'\n';
                }*/

            }
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            //if(isAxisymmetric)//Assuming axis of symmetry is along y = 0
            //if ( isAxisymmetric && lobc_v[1] == amrex::LinOpBCType::Neumann)
            //{
            //    for (int shift = 1; shift <= Nghost; shift++)
            //    {
            //        for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            //            v(i, j - shift, k) = -1.0*v(i, j + shift , k);
            //    }
            //}
            //else
            {
                for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                {
                    x = prob_lo[0] + dx[0] * (i + 0.5);
                    if (lobc_v[1] == amrex::LinOpBCType::Dirichlet)
                    {
                        v(i, j, k) = v_bcf[1](x, y, time);
                    }
                    else if (lobc_v[1] == amrex::LinOpBCType::Neumann)
                    {
                        //v(i, j, k) = v(i, j + 1, k) - dx[1] * v_bcf[1](x, y, time);
                        v(i, j - 1, k) = v(i, j + 1, k) - 2.0 * dx[1] * v_bcf[1](x, y, time);
                    }
                    else
                    {
                        amrex::Abort("invalid bc for v on ylo face");
                    }
                }
            }
        }

        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * j;
                if (hibc_v[0] == amrex::LinOpBCType::Dirichlet)
                {
                    v(i + 1, j, k) = 2.0 * v_bcf[2](x, y, time) - v(i, j, k);
                }
                else if (hibc_v[0] == amrex::LinOpBCType::Neumann)
                {
                    v(i + 1, j, k) = v(i, j, k) + dx[0] * v_bcf[2](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for v on xhi face");
                }
            }
        }

        /// top
        if (bx.bigEnd(1) > domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (hibc_v[1] == amrex::LinOpBCType::Dirichlet)
                {
                    v(i, j, k) = v_bcf[3](x, y, time);
                }
                else if (hibc_v[1] == amrex::LinOpBCType::Neumann)
                {
                    //v(i, j, k) = v(i, j - 1, k) + dx[1] * v_bcf[3](x, y, time);
                    v(i, j + 1, k) = v(i, j - 1, k) + 2.0 * dx[1] * v_bcf[3](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for v on yhi face");
                }
                /*if( i == 350 && j > 500)
                {
                    amrex::Print()<<"v Boundary condition "<<'\n';
                    amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
                    amrex::Print()<<"v("<<i<<","<< j<<", "<<k<<") = "<<v(i, j, k)<<'\n';
                    amrex::Print()<<"v("<<i<<", "<<j + 1<<", "<<k<<") = "<<v(i, j + 1, k)<<'\n';
                    amrex::Print()<<"*****************************"<<'\n';
                }*/
            }
        }
    }
}

void incFSI::PressureBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf)
{
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &P = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (lobc_p[0] == amrex::LinOpBCType::Dirichlet)
                {
                    P(i - 1, j, k) = 2.0 * p_bcf[0](x, y, time) - P(i, j, k);
                }
                else if (lobc_p[0] == amrex::LinOpBCType::Neumann)
                {
                    P(i - 1, j, k) = P(i, j, k) - dx[0] * p_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on xlo face");
                }
            }
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            //if(isAxisymmetric)//Assuming axis of symmetry is along y = 0
            //if ( isAxisymmetric && lobc_p[1] == amrex::LinOpBCType::Neumann)
            //{
            //    for (int shift = 1; shift <= Nghost; shift++)
            //    {
            //        for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            //            P(i, j - shift, k) = P(i, j + shift - 1, k);
            //    }
            //}
            //else
            {
                for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                {
                    x = prob_lo[0] + dx[0] * (i + 0.5);
                    if (lobc_p[1] == amrex::LinOpBCType::Dirichlet)
                    {
                        P(i, j - 1, k) = 2.0 * p_bcf[1](x, y, time) - P(i, j, k);
                    }
                    else if (lobc_p[1] == amrex::LinOpBCType::Neumann)
                    {
                        P(i, j - 1, k) = P(i, j, k) - dx[1] * p_bcf[1](x, y, time);
                    }
                    else
                    {
                        amrex::Print()<<" bc = "<<lobc_p[1]<<'\n';
                        amrex::Abort("invalid bc for pressure on ylo face");
                    }
                }
            }
        }

        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (hibc_p[0] == amrex::LinOpBCType::Dirichlet)
                {
                    P(i + 1, j, k) = 2.0 * p_bcf[2](x, y, time) - P(i, j, k);
                }
                else if (hibc_p[0] == amrex::LinOpBCType::Neumann)
                {
                    P(i + 1, j, k) = dx[0] * p_bcf[2](x, y, time) + P(i, j, k);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on xhi face");
                }
            }
        }

        /// top
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (hibc_p[1] == amrex::LinOpBCType::Dirichlet)
                {
                    P(i, j + 1, k) = 2.0 * p_bcf[3](x, y, time) - P(i, j, k);
                }
                else if (hibc_p[1] == amrex::LinOpBCType::Neumann)
                {
                    P(i, j + 1, k) = dx[1] * p_bcf[3](x, y, time) + P(i, j, k);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on yhi face");
                }
            }        
        }
    }
}

void incFSI::TemperatureBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf)
{
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &T = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (lobc_T[0] == amrex::LinOpBCType::Dirichlet)
                {
                    T(i - 1, j, k) = 2.0 * T_bcf[0](x, y, time) - T(i, j, k);
                }
                else if (lobc_T[0] == amrex::LinOpBCType::Neumann)
                {
                    T(i - 1, j, k) = T(i, j, k) - dx[0] * T_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for temperature on xlo face");
                }
            }
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            //if(isAxisymmetric)//Assuming axis of symmetry is along y = 0
            //if ( isAxisymmetric && lobc_T[1] == amrex::LinOpBCType::Neumann)
            //{
            //    for (int shift = 1; shift <= Nghost; shift++)
            //    {
            //        for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            //            T(i, j - shift, k) = T(i, j + shift - 1, k);
            //    }
            //}
            //else
            {
                for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                {
                    x = prob_lo[0] + dx[0] * (i + 0.5);
                    if (lobc_T[1] == amrex::LinOpBCType::Dirichlet)
                    {
                        T(i, j - 1, k) = 2.0 * T_bcf[1](x, y, time) - T(i, j, k);
                    }
                    else if (lobc_T[1] == amrex::LinOpBCType::Neumann)
                    {
                        T(i, j - 1, k) = T(i, j, k) - dx[1] * T_bcf[1](x, y, time);
                    }
                    else
                    {
                        amrex::Print()<<" bc = "<<lobc_p[1]<<'\n';
                        amrex::Abort("invalid bc for temperature on ylo face");
                    }
                }
            }
        }

        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (hibc_T[0] == amrex::LinOpBCType::Dirichlet)
                {
                    T(i + 1, j, k) = 2.0 * T_bcf[2](x, y, time) - T(i, j, k);
                }
                else if (hibc_T[0] == amrex::LinOpBCType::Neumann)
                {
                    T(i + 1, j, k) = dx[0] * T_bcf[2](x, y, time) + T(i, j, k);
                }
                else
                {
                    amrex::Abort("invalid bc for temperature on xhi face");
                }
            }
        }

        /// top
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (hibc_T[1] == amrex::LinOpBCType::Dirichlet)
                {
                    T(i, j + 1, k) = 2.0 * T_bcf[3](x, y, time) - T(i, j, k);
                }
                else if (hibc_T[1] == amrex::LinOpBCType::Neumann)
                {
                    T(i, j + 1, k) = dx[1] * T_bcf[3](x, y, time) + T(i, j, k);
                }
                else
                {
                    amrex::Abort("invalid bc for temperature on yhi face");
                }
            }        
        }
    }
}

void incFSI::PhiBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf)
{
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &phi = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (lobc_phi[0] == amrex::LinOpBCType::Dirichlet)
                {
                    phi(i - 1, j, k) = 2.0 * Phi_bcf[0](x, y, time) - phi(i, j, k);
		    phi(i - 2, j, k) = 2.0 * Phi_bcf[0](x, y, time) - phi(i + 1, j, k);
		    phi(i - 3, j, k) = 2.0 * Phi_bcf[0](x, y, time) - phi(i + 2, j, k);
		    phi(i - 4, j, k) = 2.0 * Phi_bcf[0](x, y, time) - phi(i + 3, j, k);
                }
                else if (lobc_phi[0] == amrex::LinOpBCType::Neumann)
                {
                    phi(i - 1, j, k) = phi(i, j, k) - dx[0] * Phi_bcf[0](x, y, time);
		    phi(i - 2, j, k) = phi(i + 1, j, k) - dx[0] * Phi_bcf[0](x, y, time);
		    phi(i - 3, j, k) = phi(i + 2, j, k) - dx[0] * Phi_bcf[0](x, y, time);
		    phi(i - 4, j, k) = phi(i + 3, j, k) - dx[0] * Phi_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for phi on xlo face");
                }
            }
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            //if(isAxisymmetric)//Assuming axis of symmetry is along y = 0
            //if ( isAxisymmetric && lobc_p[1] == amrex::LinOpBCType::Neumann)
            //{
            //    for (int shift = 1; shift <= Nghost; shift++)
            //    {
            //        for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            //            phi(i, j - shift, k) = phi(i, j + shift - 1, k);
            //    }
            //}
            //else
            {
                for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                {
                    x = prob_lo[0] + dx[0] * (i + 0.5);
                    if (lobc_phi[1] == amrex::LinOpBCType::Dirichlet)
                    {
                        phi(i, j - 1, k) = 2.0 * Phi_bcf[1](x, y, time) - phi(i, j, k);
			phi(i, j - 2, k) = 2.0 * Phi_bcf[1](x, y, time) - phi(i, j + 1, k);
			phi(i, j - 3, k) = 2.0 * Phi_bcf[1](x, y, time) - phi(i, j + 2, k);
			phi(i, j - 4, k) = 2.0 * Phi_bcf[1](x, y, time) - phi(i, j + 3, k);
                    }
                    else if (lobc_phi[1] == amrex::LinOpBCType::Neumann)
                    {
                        phi(i, j - 1, k) = phi(i, j, k) - dx[1] * Phi_bcf[1](x, y, time);
			phi(i, j - 2, k) = phi(i, j + 1, k) - dx[1] * Phi_bcf[1](x, y, time);
			phi(i, j - 3, k) = phi(i, j + 2, k) - dx[1] * Phi_bcf[1](x, y, time);
			phi(i, j - 4, k) = phi(i, j + 3, k) - dx[1] * Phi_bcf[1](x, y, time);
                    }
                    else
                    {
                        amrex::Print()<<" bc = "<<lobc_phi[1]<<'\n';
                        amrex::Abort("invalid bc for phi on ylo face");
                    }
                }
            }
        }

        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (hibc_phi[0] == amrex::LinOpBCType::Dirichlet)
                {
                    phi(i + 1, j, k) = 2.0 * Phi_bcf[2](x, y, time) - phi(i, j, k);
		    phi(i + 2, j, k) = 2.0 * Phi_bcf[2](x, y, time) - phi(i - 1, j, k);
		    phi(i + 3, j, k) = 2.0 * Phi_bcf[2](x, y, time) - phi(i - 2, j, k);
		    phi(i + 4, j, k) = 2.0 * Phi_bcf[2](x, y, time) - phi(i - 3, j, k);
                }
                else if (hibc_phi[0] == amrex::LinOpBCType::Neumann)
                {
                    phi(i + 1, j, k) = dx[0] * Phi_bcf[2](x, y, time) + phi(i, j, k);
		    phi(i + 2, j, k) = dx[0] * Phi_bcf[2](x, y, time) + phi(i - 1, j, k);
		    phi(i + 3, j, k) = dx[0] * Phi_bcf[2](x, y, time) + phi(i - 2, j, k);
		    phi(i + 4, j, k) = dx[0] * Phi_bcf[2](x, y, time) + phi(i - 3, j, k);
                }
                else
                {
                    amrex::Abort("invalid bc for phi on xhi face");
                }
            }
        }

        /// top
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (hibc_phi[1] == amrex::LinOpBCType::Dirichlet)
                {
                    phi(i, j + 1, k) = 2.0 * Phi_bcf[3](x, y, time) - phi(i, j, k);
		    phi(i, j + 2, k) = 2.0 * Phi_bcf[3](x, y, time) - phi(i, j - 1, k);
		    phi(i, j + 3, k) = 2.0 * Phi_bcf[3](x, y, time) - phi(i, j - 2, k);
		    phi(i, j + 4, k) = 2.0 * Phi_bcf[3](x, y, time) - phi(i, j - 3, k);
                }
                else if (hibc_phi[1] == amrex::LinOpBCType::Neumann)
                {
                    phi(i, j + 1, k) = dx[1] * Phi_bcf[3](x, y, time) + phi(i, j, k);
		    phi(i, j + 2, k) = dx[1] * Phi_bcf[3](x, y, time) + phi(i, j - 1, k);
		    phi(i, j + 3, k) = dx[1] * Phi_bcf[3](x, y, time) + phi(i, j - 2, k);
		    phi(i, j + 4, k) = dx[1] * Phi_bcf[3](x, y, time) + phi(i, j - 3, k);
                }
                else
                {
                    amrex::Abort("invalid bc for Phi on yhi face");
                }
            }        
        }
    }
}


void incFSI::ScalarBoundaryConditions(int lev,int iscalar, amrex::Real time, amrex::MultiFab& mf)
{
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &scalar = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (lobc_Scalars[iscalar][0] == amrex::LinOpBCType::Dirichlet)
                {
                    scalar(i - 1, j, k) = 2.0 * Scalars_bcf[iscalar][0](x, y, time) - scalar(i, j, k);
                }
                else if (lobc_Scalars[iscalar][0] == amrex::LinOpBCType::Neumann)
                {
                    scalar(i - 1, j, k) = scalar(i, j, k) - dx[0] * Scalars_bcf[iscalar][0](x, y, time);
                }
                else
                {
	            amrex::Print()<<"iscalar = "<<iscalar<<'\n';
                    amrex::Abort("invalid bc for scalar on xlo face");
                }
            }
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            if(isAxisymmetric)//Assuming axis of symmetry is along y = 0
            {
                for (int shift = 1; shift <= Nghost; shift++)
                {
                    for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                    {
			if(iscalar == F11_num || iscalar == F22_num || iscalar == F33_num)
                            scalar(i, j - shift, k) = scalar(i, j + shift - 1, k);
			else if(iscalar == F12_num || iscalar == F21_num)
                            scalar(i, j - shift, k) = -1.0*scalar(i, j + shift - 1, k);
		    }
                }
            }
            else
            {
                for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
                {
                    x = prob_lo[0] + dx[0] * (i + 0.5);
                    if (lobc_Scalars[iscalar][1] == amrex::LinOpBCType::Dirichlet)
                    {
                        scalar(i, j - 1, k) = 2.0 * Scalars_bcf[iscalar][1](x, y, time) - scalar(i, j, k);
                    }
                    else if (lobc_Scalars[iscalar][1] == amrex::LinOpBCType::Neumann)
                    {
                        scalar(i, j - 1, k) = scalar(i, j, k) - dx[1] * Scalars_bcf[iscalar][1](x, y, time);
                    }
                    else
                    {
                        amrex::Print()<<" bc = "<<lobc_Scalars[1]<<'\n';
                        amrex::Abort("invalid bc for scalar on ylo face");
                    }
                }
            }
        }

        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost; j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (hibc_Scalars[iscalar][0] == amrex::LinOpBCType::Dirichlet)
                {
                    scalar(i + 1, j, k) = 2.0 * Scalars_bcf[iscalar][2](x, y, time) - scalar(i, j, k);
                }
                else if (hibc_Scalars[iscalar][0] == amrex::LinOpBCType::Neumann)
                {
                    scalar(i + 1, j, k) = scalar(i, j, k) + dx[0] * Scalars_bcf[iscalar][2](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for scalar on xhi face");
                }
            }
        }

        /// top
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (hibc_Scalars[iscalar][1] == amrex::LinOpBCType::Dirichlet)
                {
                    scalar(i, j + 1, k) = 2.0 * Scalars_bcf[iscalar][3](x, y, time) - scalar(i, j, k);
                }
                else if (hibc_Scalars[iscalar][1] == amrex::LinOpBCType::Neumann)
                {
                    scalar(i, j + 1, k) = scalar(i, j, k) + dx[1] * Scalars_bcf[iscalar][3](x, y, time);
                }
                else
                {
                    amrex::Print()<<"iscalar = "<<iscalar<<'\n';
                    amrex::Abort("invalid bc for scalar on yhi face");
                }
            }        
        }
    }
}

void incFSI::CollocatedVelocityBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf)
{
    mf.FillBoundary();

    const amrex::Box &domain = geom[lev].Domain();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real* prob_hi = geom[lev].ProbHi();
    const amrex::Real* dx = geom[lev].CellSize();

    amrex::Real x, y;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &Vel = mf.array(mfi);
        int k = 0; /// as 2d

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            x = prob_lo[0];
            int i = bx.smallEnd(0);
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (lobc_u[0] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i - 1, j, k, 0) = 2.0 * u_bcf[0](x, y, time) - Vel(i, j, k, 0);
                }
                else if (lobc_u[0] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i - 1, j, k, 0) = Vel(i, j, k, 0) - dx[0] * u_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on xlo face");
                }
                if (lobc_v[0] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i - 1, j, k, 1) = 2.0 * v_bcf[0](x, y, time) - Vel(i, j, k, 1);
                }
                else if (lobc_v[0] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i - 1, j, k, 1) = Vel(i, j, k, 1) - dx[0] * v_bcf[0](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on xlo face");
                }
            }
    
        }

        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            y = prob_lo[1];
            int j = bx.smallEnd(1);
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (lobc_u[1] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i, j - 1, k, 0) = 2.0 * u_bcf[1](x, y, time) - Vel(i, j, k, 0);
                }
                else if (lobc_u[1] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i, j - 1, k, 0) = Vel(i, j, k, 0) - dx[1] * u_bcf[1](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on ylo face");
                }
                
                if (lobc_v[1] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i, j - 1, k, 1) = 2.0 * v_bcf[1](x, y, time) - Vel(i, j, k, 1);
                }
                else if (lobc_v[1] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i, j - 1, k, 1) = Vel(i, j, k, 1) - dx[1] * v_bcf[1](x, y, time);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on ylo face");
                }               
            }
            //if ( isAxisymmetric && lobc_v[1] == amrex::LinOpBCType::Neumann)
            //{
            //    for (int shift = 1; shift <= Nghost; shift++)
            //    {
            //        for (int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost; i++)
            //        {
            //            Vel(i, j - shift, k, 0) =  Vel(i, j + shift - 1, k, 0);
            //            Vel(i, j - shift, k, 1) =  -1.0*Vel(i, j + shift - 1, k, 1);
            //        }
            //    }
            //}
        }

        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            x = prob_hi[0];
            int i = bx.bigEnd(0);
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
            {
                y = prob_lo[1] + dx[1] * (j + 0.5);
                if (hibc_u[0] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i + 1, j, k, 0) = 2.0 * u_bcf[2](x, y, time) - Vel(i, j, k, 0);
                }
                else if (hibc_u[0] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i + 1, j, k, 0) = dx[0] * u_bcf[2](x, y, time) + Vel(i, j, k, 0);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on xhi face");
                }
                
                if (hibc_v[0] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i + 1, j, k, 1) = 2.0 * v_bcf[2](x, y, time) - Vel(i, j, k, 1);
                }
                else if (hibc_v[0] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i + 1, j, k, 1) = dx[0] * v_bcf[2](x, y, time) + Vel(i, j, k, 1);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on xhi face");
                }
            }
        }

        /// top
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            y = prob_hi[1];
            int j = bx.bigEnd(1);
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
            {
                x = prob_lo[0] + dx[0] * (i + 0.5);
                if (hibc_u[1] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i, j + 1, k, 0) = 2.0 * u_bcf[3](x, y, time) - Vel(i, j, k, 0);
                }
                else if (hibc_u[1] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i, j + 1, k, 0) = dx[1] * u_bcf[1](x, y, time) + Vel(i, j, k, 0);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on yhi face");
                }
                
                if (hibc_v[1] == amrex::LinOpBCType::Dirichlet)
                {
                    Vel(i, j + 1, k, 1) = 2.0 * v_bcf[3](x, y, time) - Vel(i, j, k, 1);
                }
                else if (hibc_v[1] == amrex::LinOpBCType::Neumann)
                {
                    Vel(i, j + 1, k, 1) = dx[1] * v_bcf[1](x, y, time) + Vel(i, j, k, 1);
                }
                else
                {
                    amrex::Abort("invalid bc for pressure on yhi face");
                }
            }        
        }
    }
}

} // namespace mycode
