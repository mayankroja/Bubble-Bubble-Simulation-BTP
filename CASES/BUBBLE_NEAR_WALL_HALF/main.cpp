#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <incFSI.H>

int main(int argc, char** argv) 
{
    amrex::Initialize(argc,argv);
    {
        amrex::Real npi = 2.0*M_PI;
        amrex::Real npi2 = npi * npi;

        amrex::Real start_time = amrex::second();

        std::vector< std::unique_ptr<mycode::interface_props> > bubble(1);

        bubble[0] = std::make_unique<mycode::interface_props>();

        mycode::incFSI NS(bubble);
   
        

        /// set boundary functions 
        NS.SetXVelBoundaryFunction
        (
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real 
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);
                amrex::Real ur = bubble[0]->Int_R*bubble[0]->Int_R*bubble[0]->Int_R_dot/(r*r);
                //return ur*x/r;
                return 0.0;

            }
        );
        NS.SetXVelBoundaryFunction
        (
            1,
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);
                amrex::Real ur = bubble[0]->Int_R*bubble[0]->Int_R*bubble[0]->Int_R_dot/(r*r);
                //return ur*x/r;
                return 0.0;

            }
        );

        NS.SetYVelBoundaryFunction
        (
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real 
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);
                amrex::Real ur = bubble[0]->Int_R*bubble[0]->Int_R*bubble[0]->Int_R_dot/(r*r);
                //return ur*y/r;
                return 0.0;

            }
        );

        NS.SetYVelBoundaryFunction
        (
            1,
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);
                amrex::Real ur = bubble[0]->Int_R*bubble[0]->Int_R*bubble[0]->Int_R_dot/(r*r);
                //return ur*y/r;
                return 0.0;

            }
        );


        NS.SetPressureBoundaryFunction
        (
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);

                //return (1.0 + (bubble[0]->Int_p - 1.0)*bubble[0]->Int_R/r +
                //       0.5*(bubble[0]->Int_R*bubble[0]->Int_R_dot*bubble[0]->Int_R_dot)/r -
                //       0.5*1.0*pow(bubble[0]->Int_R*bubble[0]->Int_R*bubble[0]->Int_R_dot/(r*r),2.0));
                return 1.0;
            }
        );

        NS.SetPressureBoundaryFunction
        (0,
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);

                return 0.0;
            }
        );


        NS.SetPressureBoundaryFunction
        (1,
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);

                return 0.0;
            }
        );


        //! on boundary condition for temperature field
        NS.SetTemperatureBoundaryFunction
        (
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

        //! on boundary condition for phase field
        NS.SetPhiBoundaryFunction
        (
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

/*	NS.SetScalarsBoundaryFunction
        (
            0, 
	    0,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            0,
            1,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            0,
            2,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 1.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            0,
            3,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            1,
            0,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            1,
            1,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            1,
            2,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            1,
            3,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 1.0;
            }
        );
	

	NS.SetScalarsBoundaryFunction
        (
            2,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            2,
	    0,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 1.0;
            }
        );

	NS.SetScalarsBoundaryFunction
        (
            3,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            4,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            5,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            5,
	    0,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 1.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            6,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        NS.SetScalarsBoundaryFunction
        (
            6,
	    0,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 1.0;
            }
        );

        NS.SetScalarsBoundaryFunction
        (
            7,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
	NS.SetScalarsBoundaryFunction
        (
            8,
            [](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                return 0.0;
            }
        );
        */
        NS.setInitPressure
        (
            [&](amrex::Real x, amrex::Real y, amrex::Real t) -> amrex::Real
            {
                amrex::Real r = std::hypot(x-bubble[0]->Int_xcp, y-bubble[0]->Int_ycp);

                //return (1.0 + (bubble[0]->Int_p - 1.0)*bubble[0]->Int_R/r);

                return 1.0;
            }
        );

        //! Initialise Phase field
	/*
        {
            NS.setInitPhi
            (
                [=](amrex::Real x, amrex::Real y, amrex::Real ds) -> amrex::Real
                {
                    amrex::Real xp = -0.5;
                    amrex::Real x_phi = x - xp;
                    //return 0.5 - 0.5 * tanh(x_phi / (2.0 * ds));;
                    return 1.0;
                }
            );
        }
        */


        NS.Initialize();
        NS.Evolve();

        amrex::Real stop_time = amrex::second() - start_time;
        const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

        amrex::Print() << "Run time = " << stop_time << std::endl;
    }
    amrex::Finalize();
    return 0;
}
