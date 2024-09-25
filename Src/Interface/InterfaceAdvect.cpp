#include "InterfaceAdvect.H"
#include <AMReX_ParmParse.H>
#include <AdvectLS_helper.H>

namespace mycode
{

InterfaceAdvect::InterfaceAdvect
(
    const amrexMesh &mesh,
    const std::string& name
) 
: 
Interface(mesh, name)
{
    Mask_.define(mesh.grid(), mesh.dmap(), 1, 0);
    Index_.define(mesh.grid(), mesh.dmap());
    Source.define(mesh.grid(), mesh.dmap());

    const amrex::Real *dx = mesh.geometry().CellSize();
    L_BETA = 4.0 * dx[0];
    L_GAMMA = 9.0 * dx[0];

    //is_advect_ls = true;
    //Initialize Mask_
    Mask_.setVal(0);
}

InterfaceAdvect::~InterfaceAdvect() {}

void InterfaceAdvect::TubeIdentification()
{
    PhaseFieldBC();
    Mask_.setVal(0);
    NTube_ = 0;
    TTube_ = 0;

    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    const amrex::Real *prob_hi = mesh_.geometry().ProbHi();
    //const amrex::Box &domain =mesh_.geometry().Domain();

    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<int> const &mask = Mask_.array(mfi);
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);

        auto &index = Index_[mfi];

        //! NOTE : before adding cell index to vector
        //!        we should clear already present data if any
        //!        otherwise in each iteration new cell index will be added to vector
        //!        and it will keep on growing
        if (!index.empty())
        {
            index.clear();
            index.shrink_to_fit();
        }
        index.reserve(100);

        //amrex::Box gbx(bx);
        //gbx.grow(1);

        /// take intersection with domain box as we do not want to compute for cells outside domain
        //amrex::Box isect = gbx & domain;


        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if ( fabs(Psi(i,j,k)) < L_BETA)
            {
                mask(i, j, k) = 3;
                index.emplace_back(AMREX_D_DECL(i, j, k));
            }
            else if (fabs(Psi(i,j,k)) < L_GAMMA)
            {
                mask(i, j, k) = 2;
                index.emplace_back(AMREX_D_DECL(i, j, k));
            }            
        });
        TTube_ += index.size();

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::IntVect iv(AMREX_D_DECL(i,j,k));
            amrex::Box onelayer(iv, iv);
            onelayer.grow(1);
            if (mask(i,j,k) == 0)
            {
                if (nbrTTube(Psi,onelayer))
                {
                    mask(i, j, k) = 1;
                    index.emplace_back(iv);
                }
            }
        });
        NTube_ += index.size();

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (mask(i,j,k) == 0)
            {
                if (Psi(i,j,k) > 0.0)
                {
                    Psi(i, j, k) = prob_hi[0] - prob_lo[0];
                }
                else
                {
                    Psi(i, j, k) = -(prob_hi[0] - prob_lo[0]);
                }
            }
        });
    }
    PhaseFieldBC();

    Mask_.FillBoundary();

    NTubes_global = NTube_;
    amrex::ParallelDescriptor::ReduceIntSum(NTubes_global);

    // amrex::Print() << Name << " NTubes = " << NTubes_global << "\n";
}

void InterfaceAdvect::TubeAdvectLevelSet(const amrex::MultiFab &xvel, const amrex::MultiFab &yvel, const amrex::Real &dt)
{
    amrex::Real u_Mid, v_Mid, Psix, Psiy;
    amrex::Real Psix_R, Psix_L, Psiy_L, Psiy_R;

    const amrex::Real *dx = mesh_.geometry().CellSize();
    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        amrex::Array4<int> const &mask = Mask_.array(mfi);
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real const> const &u = xvel.const_array(mfi);
        amrex::Array4<amrex::Real const> const &v = yvel.const_array(mfi);
        auto &index = Index_[mfi];
        auto &source = Source[mfi];
        source.resize(index.size());

        int n = 0;
        for (auto &&cell : index)
        {
            int i = cell[0], j = cell[1], k = 0;
            if (mask(cell) != 3 && mask(cell) != 2 && mask(cell) != 1)
            {
                amrex::PrintToFile("log") << "Trouble in defining mask for tube T \n";
                amrex::PrintToFile("log") << mask(cell) << "\t" << cell << "\n";
                exit(1);
            }
            u_Mid = 0.5 * (u(i, j, k) + u(i + 1, j, k));
            v_Mid = 0.5 * (v(i, j, k) + v(i, j + 1, k));

            if (mask(cell) == 3)
            {
                WENO5_LS(Psix_L, Psix_R, Psiy_L, Psiy_R, i, j, dx, Psi);
            }
            else if (mask(cell) == 2 || mask(cell) == 1)
            {
                Psix_R = (Psi(i + 1, j, k) - Psi(i, j, k)) / dx[0];
                Psix_L = (Psi(i, j, k) - Psi(i - 1, j, k)) / dx[0];
                Psiy_R = (Psi(i, j + 1, k) - Psi(i, j, k)) / dx[1];
                Psiy_L = (Psi(i, j, k) - Psi(i, j - 1, k)) / dx[1];
            }
            if (u_Mid > 0.0)
                Psix = Psix_L;
            else
                Psix = Psix_R;
            if (v_Mid > 0.0)
                Psiy = Psiy_L;
            else
                Psiy = Psiy_R;

            source[n] = -dt * (u_Mid * Psix + v_Mid * Psiy);
            n++;
        }
    }

    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        auto &index = Index_[mfi];
        auto &source = Source[mfi];

        int n = 0;
        for (auto &&cell : index)
        {
            Psi(cell) += source[n];
            n++;
        }
    }

    PhaseFieldBC();
}

void InterfaceAdvect::TubeAdvectLevelSet_RK2(const amrex::MultiFab &xvel, const amrex::MultiFab &yvel, const amrex::Real &dt)
{
    amrex::Real u_Mid, v_Mid, Psix, Psiy;
    amrex::Real Psix_R, Psix_L, Psiy_L, Psiy_R;

    const amrex::Real *dx = mesh_.geometry().CellSize();
    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        amrex::Array4<int> const &mask = Mask_.array(mfi);
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real const> const &u = xvel.const_array(mfi);
        amrex::Array4<amrex::Real const> const &v = yvel.const_array(mfi);
        auto &index = Index_[mfi];
        auto &source = Source[mfi];
        source.resize(index.size());

        int n = 0;
        for (auto &&cell : index)
        {
            int i = cell[0], j = cell[1], k = 0;
            if (mask(cell) != 3 && mask(cell) != 2 && mask(cell) != 1)
            {
                amrex::PrintToFile("log") << "Trouble in defining mask for tube T \n";
                amrex::PrintToFile("log") << mask(cell) << "\t" << cell << "\n";
                exit(1);
            }
            u_Mid = 0.5 * (u(i, j, k) + u(i + 1, j, k));
            v_Mid = 0.5 * (v(i, j, k) + v(i, j + 1, k));

            if (mask(cell) == 3)
            {
                WENO5_LS(Psix_L, Psix_R, Psiy_L, Psiy_R, i, j, dx, Psi);
            }
            else if (mask(cell) == 2 || mask(cell) == 1)
            {
                Psix_R = (Psi(i + 1, j, k) - Psi(i, j, k)) / dx[0];
                Psix_L = (Psi(i, j, k) - Psi(i - 1, j, k)) / dx[0];
                Psiy_R = (Psi(i, j + 1, k) - Psi(i, j, k)) / dx[1];
                Psiy_L = (Psi(i, j, k) - Psi(i, j - 1, k)) / dx[1];
            }
            if (u_Mid > 0.0)
                Psix = Psix_L;
            else
                Psix = Psix_R;
            if (v_Mid > 0.0)
                Psiy = Psiy_L;
            else
                Psiy = Psiy_R;

            source[n] = -dt * (u_Mid * Psix + v_Mid * Psiy);
            n++;
        }
    }

    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real> const &frk_Psi = FRK_psi.array(mfi);
        auto &index = Index_[mfi];
        auto &source = Source[mfi];

        int n = 0;
        for (auto &&cell : index)
        {
            Psi(cell) = frk_Psi(cell) + source[n];
            n++;
        }
    }

    PhaseFieldBC();
}
/*
void InterfaceAdvect::Regularization()
{
    const amrex::Real *dx = mesh_.geometry().CellSize();
    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        auto &index = Index_[mfi];
        for (auto &&cell : index)
        {
            int i = cell[0], j = cell[1], k = 0;
            if(cell[1] == 0)
            {
                //amrex::Print()<<"cell = "<<cell<<'\n';
                //amrex::Print()<<"Psi(cell) = "<<Psi(cell)<<'\n';
                 
                if(Psi(cell) < 0.0 && std::abs(Psi(cell)) < dx[1])
                {
		    ////Regularization 1 Allow interface to break at the axis    
                    //if(Psi(i, 1 , 0) > 0.0 )
                    //{
                    //    amrex::Print()<<"Regularization applied to cell = "<<cell<<" , psi = "<<Psi(cell)<<'\n';
                    //    Psi(cell) = 0.1*dx[1];
                    //}
		    
		    else if(Psi(i, 3 , 0) > 0.0)
	            {
                        amrex::Print()<<"Regularization applied to cell = "<<cell<<" , psi = "<<Psi(cell)<<'\n';
                        Psi(cell) = 0.1*dx[1];		    
			Psi(i , j + 1 , k) = 1.1*dx[1];
			Psi(i , j + 2 , k) = 2.1*dx[1];
		    }
		    
		    else if(Psi(i, 7 , 0) > 0.0)
                    {
                        amrex::Print()<<"Regularization applied to cell = "<<cell<<" , psi = "<<Psi(cell)<<'\n';
                        Psi(cell) = 0.1*dx[1];              
                        Psi(i , j + 1 , k) = 1.1*dx[1];
			Psi(i , j + 2 , k) = 2.1*dx[1];
			Psi(i , j + 3 , k) = 3.1*dx[1];
			Psi(i , j + 4 , k) = 4.1*dx[1];
			Psi(i , j + 5 , k) = 5.1*dx[1];
			Psi(i , j + 6 , k) = 6.1*dx[1];
                    }
		    
                }
            }
        }
	const amrex::Box& bx = mfi.validbox();
	int detect_float = 0;
	int bx_x0 = bx.smallEnd(0);
	int bx_x1 = bx.bigEnd(0);
	
        for (int i = bx_x0; i<= bx_x1; i++)
        {
	    int j = 0;
	    int k = 0; 
            {
		//Detect thin strip
                if(Psi(i, j, k) < 0.0 && Psi(i - 1, k, j) > 0)
                {
	            int thin = 0;
		    for( int jj = 1; jj <= 4; jj++ )
	            {
		        if(Psi(i , jj , k) > 0)
			{
			   thin == jj;
			   amrex::Print()<<"Thin detected : "<<i<<" , "<<jj - 1<<'\n';
                           for(int jjj = 0; jjj < jj ; jjj++)
                           {
                               //Psi(i , jjj, k) = 10.0;
                               amrex::Print(-1)<<"Regularization applied to cell = "<<i<<""<<" , psi = "<<Psi(i , jjj, k)<<'\n';
                           }
			   break;
			}
		    }
                }
            }
        }
        
    }

}
*/
void InterfaceAdvect::TubeAdvectLevelSet_RK3(int RKStage, const amrex::MultiFab &xvel, const amrex::MultiFab &yvel, const amrex::Real &dt)
{
    amrex::Real u_Mid, v_Mid, Psix, Psiy;
    amrex::Real Psix_R, Psix_L, Psiy_L, Psiy_R;

    const amrex::Real *dx = mesh_.geometry().CellSize();
    
    amrex::MultiFab::Copy(psi, FRK_psi, 0, 0, 1, nghost);//copy Psi from previous time step to psi
    
    for (int step = 1; step <= RKStage; step++)
    {
        for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
        {
            amrex::Array4<int> const &mask = Mask_.array(mfi);
            //amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
            amrex::Array4<amrex::Real> Psi;
            if(step == 1)
                Psi = FRK_psi.array(mfi);
            else if(step == 2)
                Psi = RK1_psi.array(mfi);
            else if(step == 3)
                Psi = RK2_psi.array(mfi);
            else if(step == 4)
                Psi = RK3_psi.array(mfi);
            
            
            amrex::Array4<amrex::Real const> const &u = xvel.const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = yvel.const_array(mfi);
            auto &index = Index_[mfi];
            auto &source = Source[mfi];
            source.resize(index.size());
    
            int n = 0;
            for (auto &&cell : index)
            {
                int i = cell[0], j = cell[1], k = 0;
                if (mask(cell) != 3 && mask(cell) != 2 && mask(cell) != 1)
                {
                    amrex::PrintToFile("log") << "Trouble in defining mask for tube T \n";
                    amrex::PrintToFile("log") << mask(cell) << "\t" << cell << "\n";
                    exit(1);
                }
                u_Mid = 0.5 * (u(i, j, k) + u(i + 1, j, k));
                v_Mid = 0.5 * (v(i, j, k) + v(i, j + 1, k));
    
                if (mask(cell) == 3)
                {
                    WENO5_LS(Psix_L, Psix_R, Psiy_L, Psiy_R, i, j, dx, Psi);
                }
                else if (mask(cell) == 2 || mask(cell) == 1)
                {
                    Psix_R = (Psi(i + 1, j, k) - Psi(i, j, k)) / dx[0];
                    Psix_L = (Psi(i, j, k) - Psi(i - 1, j, k)) / dx[0];
                    Psiy_R = (Psi(i, j + 1, k) - Psi(i, j, k)) / dx[1];
                    Psiy_L = (Psi(i, j, k) - Psi(i, j - 1, k)) / dx[1];
                }
                if (u_Mid > 0.0)
                    Psix = Psix_L;
                else
                    Psix = Psix_R;
                if (v_Mid > 0.0)
                    Psiy = Psiy_L;
                else
                    Psiy = Psiy_R;
    
                source[n] = -dt * (u_Mid * Psix + v_Mid * Psiy);
                n++;
            }
        }
    
        for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
        {
            amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
            //amrex::Array4<amrex::Real> const &frk_Psi = FRK_psi.array(mfi);
            auto &index = Index_[mfi];
            auto &source = Source[mfi];
    
            int n = 0;
            for (auto &&cell : index)
            {
                Psi(cell) += rk_a[RKStage - 1][step - 1] * source[n];
                n++;
            }
        }

        PhaseFieldBC();
    }
}

bool InterfaceAdvect::nbrTTube(amrex::Array4<amrex::Real> const &Psi, const amrex::Box &bx)
{
    const amrex::Box &domain = mesh_.geometry().Domain();
    bool nbr_ttube = false;
    for (amrex::BoxIterator bit(bx); bit.ok(); ++bit)
    {
        amrex::IntVect iv(bit());
        if (domain.contains(iv) && fabs(Psi(iv)) < L_GAMMA)
        {
            nbr_ttube = true;
            break;
        }
    }
    
    return nbr_ttube;
}

} /*End namespace mycode */

