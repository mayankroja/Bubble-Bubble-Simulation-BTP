#include "HypreIJ.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

HypreIJ::HypreIJ
(
    MPI_Comm comm,
    amrex::Vector<amrex::Geometry> *geom,
    amrex::Vector<amrex::BoxArray> *grid,
    amrex::Vector<amrex::DistributionMapping> *dmap,
    amrex::Vector<CFMask> *cfmask,
    amrex::Vector<amrex::iMultiFab *> &mask,
    MLTraverseIndex *tri
)
{
    LinearSolver::define(comm, geom, grid, dmap, cfmask, mask, tri);
    init();
}

HypreIJ::HypreIJ
(
    MPI_Comm comm,
    amrex::Vector<amrex::Geometry> *geom,
    amrex::Vector<amrex::BoxArray> *grid,
    amrex::Vector<amrex::DistributionMapping> *dmap,
    amrex::Vector<CFMask> *cfmask,
    MLTraverseIndex *tri
)
{
    LinearSolver::define(comm, geom, grid, dmap, cfmask, tri);
    init();
}

void HypreIJ::init()
{
    ls_lib_ = LinearSolver::HYPRE;

    amrex::ParmParse pp("hypre");
    std::string solver;
    pp.get("solver", solver);
    if (solver == "GMRES")
    {
        ls_solver_ = GMRES;
    }
    else if (solver == "AMG")
    {
        ls_solver_ = AMG;
    }
    else if (solver == "PCG")
    {
        ls_solver_ = PCG;
    }
    else if (solver == "PCG_AMG")
    {
        ls_solver_ = PCG_AMG;
    }
    else
    {
        amrex::Print() << "LinearSolver: invalid solver, options are\n";
        amrex::Print() << "\tGMRES ( Flexible GMRES with  AMG Preconditioner )\n";
        amrex::Print() << "\tAMG\n";
        amrex::Print() << "\tPCG\n";
        amrex::Print() << "\tPCG_AMG( PCG with  AMG Preconditioner )\n";
        exit(1);
    }
    preAssembleSystem();
}

HypreIJ::~HypreIJ()
{
    HYPRE_IJMatrixDestroy(HYPREMat_A);
    HYPREMat_A = NULL;
    HYPRE_IJVectorDestroy(HYPREVec_b);
    HYPREVec_b = NULL;
    HYPRE_IJVectorDestroy(HYPREVec_x);
    HYPREVec_x = NULL;

    HYPREprecond = NULL;
    HYPREsolver = NULL;
}

void HypreIJ::set_A_Val(int *ncols, int *row_id, int *col_id, amrex::Real *values)
{
    HYPRE_IJMatrixSetValues(HYPREMat_A, 1, ncols, row_id, col_id, values);
}

void HypreIJ::set_x_Val(int *row_id, amrex::Real *values)
{
    HYPRE_IJVectorSetValues(HYPREVec_x, 1, row_id, values);
}

void HypreIJ::set_b_Val(int *row_id, amrex::Real *values)
{
    HYPRE_IJVectorSetValues(HYPREVec_b, 1, row_id, values);
}

void HypreIJ::getVecVal(void *x, int nrows, int *row_id, amrex::Real *values)
{
    HYPRE_IJVectorGetValues(static_cast<HYPRE_IJVector>(x), nrows, row_id, values);
}

void HypreIJ::preAssembleSystem()
{
    int num_procs, myid;
    MPI_Comm_size(comm_, &num_procs);
    MPI_Comm_rank(comm_, &myid);

    int ilower = tridx->getProcBegin();
    int iupper = tridx->getProcEnd() - 1;

    // Ready Mat_A, Vec_x, Vec_b
    HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &HYPREMat_A);
    HYPRE_IJMatrixSetObjectType(HYPREMat_A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(HYPREMat_A);
    //
    HYPRE_IJVectorCreate(comm_, ilower, iupper, &HYPREVec_x);
    HYPRE_IJVectorSetObjectType(HYPREVec_x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(HYPREVec_x);
    //
    HYPRE_IJVectorCreate(comm_, ilower, iupper, &HYPREVec_b);
    HYPRE_IJVectorSetObjectType(HYPREVec_b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(HYPREVec_b);

    Mat_A = HYPREMat_A;
    Vec_x = HYPREVec_x;
    Vec_b = HYPREVec_b;
}

void HypreIJ::assembleSystem()
{
    HYPRE_IJMatrixAssemble(HYPREMat_A);
    HYPRE_IJVectorAssemble(HYPREVec_x);
    HYPRE_IJVectorAssemble(HYPREVec_b);
}

void HypreIJ::solverSetupAndSolve()
{
    HYPRE_IJMatrixGetObject(HYPREMat_A, (void **)&par_A);
    HYPRE_IJVectorGetObject(HYPREVec_b, (void **)&par_b);
    HYPRE_IJVectorGetObject(HYPREVec_x, (void **)&par_x);

    if (print_system)
    {
        HYPRE_IJMatrixPrint(HYPREMat_A, "IJ.out.A");
        HYPRE_IJVectorPrint(HYPREVec_b, "IJ.out.b");
        HYPRE_IJVectorPrint(HYPREVec_x, "IJ.out.x");
    }

    switch (ls_solver_)
    {
        case AMG:
            AMGSolver();
            break;
        case GMRES:
            GMRESSolver();
            break;
        case PCG:
            PCGSolver(false);
            break;
        case PCG_AMG:
            PCGSolver(true);
            break;
        default:
            amrex::PrintToFile("log") << "HYPREIJ: invalid solver option specified, using default GMRES\n";
            GMRESSolver();
            break;
    }
}

void HypreIJ::AMGSolver()
{
    //! the settings in incFSI::AMGPoissonSolve
    //! these are better than what I had previously
    HYPRE_BoomerAMGCreate(&HYPREsolver);
    HYPRE_BoomerAMGSetOldDefault(HYPREsolver);    // Falgout coarsening with modified classical interpolaiton
    HYPRE_BoomerAMGSetRelaxType(HYPREsolver, 3);  // G-S/Jacobi hybrid relaxation
    HYPRE_BoomerAMGSetRelaxOrder(HYPREsolver, 1); // uses C/F relaxation
    HYPRE_BoomerAMGSetNumSweeps(HYPREsolver, 1);  // Sweeeps on each level
    //! removing this is little bit faster
    // HYPRE_BoomerAMGSetMaxLevels(HYPREsolver, 20); // maximum number of levels
    HYPRE_BoomerAMGSetMaxIter(HYPREsolver, max_iter);
    HYPRE_BoomerAMGSetTol(HYPREsolver, rel_tolerance); // conv. tolerance
    HYPRE_BoomerAMGSetup(HYPREsolver, par_A, par_b, par_x);
    HYPRE_BoomerAMGSolve(HYPREsolver, par_A, par_b, par_x);

    // if (verbose >= 2)
    {
        HYPRE_BoomerAMGGetNumIterations(HYPREsolver, &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(HYPREsolver, &res);

        amrex::PrintToFile("log") << "\n"
                                  << "Hypre IJ BoomerAMG Iterations, Relative Residual "
                                  << num_iterations << "\t"
                                  << res << "\n";
    }

    HYPRE_BoomerAMGDestroy(HYPREsolver);
}

void HypreIJ::GMRESSolver()
{
    int restart = 30;

    /* Create solver */
    HYPRE_ParCSRFlexGMRESCreate(comm_, &HYPREsolver);

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_FlexGMRESSetKDim(HYPREsolver, restart);
    HYPRE_FlexGMRESSetMaxIter(HYPREsolver, max_iter); /* max iterations */
                                                      //      HYPRE_FlexGMRESSetPrintLevel(HYPREsolver, 2); /* print solve info */
    HYPRE_FlexGMRESSetLogging(HYPREsolver, 1);        /* needed to get run info later */

    HYPRE_FlexGMRESSetTol(HYPREsolver, rel_tolerance); /* conv. tolerance */

    /* Now set up the AMG preconditioner and specify any parameters */
    HYPRE_BoomerAMGCreate(&HYPREprecond);

    //      HYPRE_BoomerAMGSetPrintLevel(HYPREprecond, 1); /* print amg solution info */
    HYPRE_BoomerAMGSetCoarsenType(HYPREprecond, 6);
    HYPRE_BoomerAMGSetOldDefault(HYPREprecond);
    HYPRE_BoomerAMGSetRelaxType(HYPREprecond, 6); /* Sym G.S./Jacobi hybrid */
    HYPRE_BoomerAMGSetNumSweeps(HYPREprecond, 1);
    HYPRE_BoomerAMGSetTol(HYPREprecond, 0.0);   /* conv. tolerance zero */
    HYPRE_BoomerAMGSetMaxIter(HYPREprecond, 1); /* do only one iteration! */

    /* Set the FlexGMRES preconditioner */
    HYPRE_FlexGMRESSetPrecond(HYPREsolver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, HYPREprecond);

    HYPRE_ParCSRFlexGMRESSetup(HYPREsolver, par_A, par_b, par_x);
    HYPRE_ParCSRFlexGMRESSolve(HYPREsolver, par_A, par_b, par_x);

    // if (verbose >= 2)
    {
        HYPRE_FlexGMRESGetNumIterations(HYPREsolver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(HYPREsolver, &res);

        amrex::PrintToFile("log") << "\n"
                                  << "HypreIJ GMRES Iterations, Relative Residual "
                                  << num_iterations << "\t"
                                  << res << "\n";
    }

    HYPRE_ParCSRFlexGMRESDestroy(HYPREsolver);
    HYPRE_BoomerAMGDestroy(HYPREprecond);
}

void HypreIJ::PCGSolver(bool preconditioner)
{
    /* Create solver */
    HYPRE_ParCSRPCGCreate(comm_, &HYPREsolver);

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_PCGSetMaxIter(HYPREsolver, max_iter);  /* max iterations */
    HYPRE_PCGSetTol(HYPREsolver, rel_tolerance); /* conv. tolerance */
    HYPRE_PCGSetTwoNorm(HYPREsolver, 1);         /* use the two norm as the stopping criteria */
    // HYPRE_PCGSetPrintLevel(HYPREsolver, 2); /* print solve info */
    HYPRE_PCGSetLogging(HYPREsolver, 1); /* needed to get run info later */

    if (preconditioner)
    {
        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&HYPREprecond);
        // HYPRE_BoomerAMGSetPrintLevel(HYPREprecond, 1); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(HYPREprecond, 6);
        HYPRE_BoomerAMGSetOldDefault(HYPREprecond);
        HYPRE_BoomerAMGSetRelaxType(HYPREprecond, 6); /* Sym G.S./Jacobi hybrid */
        HYPRE_BoomerAMGSetNumSweeps(HYPREprecond, 1);
        HYPRE_BoomerAMGSetTol(HYPREprecond, 0.0);   /* conv. tolerance zero */
        HYPRE_BoomerAMGSetMaxIter(HYPREprecond, 1); /* do only one iteration! */

        /* Set the PCG preconditioner */
        HYPRE_PCGSetPrecond(HYPREsolver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, HYPREprecond);
    }  

    /* Now setup and solve! */
    HYPRE_ParCSRPCGSetup(HYPREsolver, par_A, par_b, par_x);
    HYPRE_ParCSRPCGSolve(HYPREsolver, par_A, par_b, par_x);

    {
        HYPRE_PCGGetNumIterations(HYPREsolver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(HYPREsolver, &res);

        amrex::PrintToFile("log") << "\n"
                                  << "HypreIJ PCG_AMG Iterations, Relative Residual "
                                  << num_iterations << "\t"
                                  << res << "\n";
    }

    HYPRE_ParCSRPCGDestroy(HYPREsolver);

    if (preconditioner)
        HYPRE_BoomerAMGDestroy(HYPREprecond);
}

} /*End namespace mycode */

