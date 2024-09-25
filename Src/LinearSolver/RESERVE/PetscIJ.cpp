#include "PetscIJ.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

PetscIJ::PetscIJ
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
    PETSC_COMM_WORLD = comm;
    init();
}

PetscIJ::PetscIJ
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
    PETSC_COMM_WORLD = comm;
    init();
}

void PetscIJ::init()
{
    ls_lib_ = LinearSolver::PETSC;
    PetscInitialize(0, 0, 0, 0);

    amrex::ParmParse pp("petsc");
    std::string solver;
    pp.get("solver", solver);
    if (solver == "GMRES")
    {
        ls_solver_ = GMRES;
    }
    else if (solver == "CG")
    {
        ls_solver_ = CG;
    }
    else if (solver == "BICGS")
    {
        ls_solver_ = BICGS;
    }
    else
    {
        amrex::Print() << "LinearSolver: invalid solver, options are\n";
        amrex::Print() << "\tGMRES\n";
        // amrex::Print() << "\tAMG\n";
        amrex::Print() << "\tCG\n";
        amrex::Print() << "\tBICGS\n";
        exit(1);
    }
    preAssembleSystem();
}

PetscIJ::~PetscIJ()
{
    MatDestroy(&PETSCMat_A);
    PETSCMat_A = nullptr;

    VecDestroy(&PETSCVec_b);
    PETSCVec_b = nullptr;

    VecDestroy(&PETSCVec_x);
    PETSCVec_x = nullptr;

    KSPDestroy(&PETSCsolver);
    PETSCsolver = nullptr;

    PetscFinalize();
}

void PetscIJ::set_A_Val(int *ncols, int *row_id, int *col_id, amrex::Real *values)
{
    MatSetValues(PETSCMat_A, 1, row_id, *ncols, col_id, values, INSERT_VALUES);
}

void PetscIJ::set_x_Val(int *row_id, amrex::Real *values)
{
    VecSetValues(PETSCVec_x, 1, row_id, values, INSERT_VALUES);
}

void PetscIJ::set_b_Val(int *row_id, amrex::Real *values)
{
    VecSetValues(PETSCVec_b, 1, row_id, values, INSERT_VALUES);
}

void PetscIJ::getVecVal(void *x, int nrows, int *row_id, amrex::Real *values)
{
    VecGetValues(static_cast<Vec>(x), nrows, row_id, values);
}

void PetscIJ::preAssembleSystem()
{
    int num_procs, myid;
    MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &myid);

    int ilower = tridx->getProcBegin();
    int iupper = tridx->getProcEnd() - 1;
    int nNodes_proc = iupper - ilower + 1;
    int nNodes_wolrd = tridx->getWorldSize();

    // create Mat_A, Vec_x, Vec_b
    {
        // estimated amount of block diag elements
        int d_nz = regular_stencil_size / 2;
        // estimated amount of block off diag elements
        int o_nz = d_nz / 2;

        MatCreate(PETSC_COMM_WORLD, &PETSCMat_A);
        MatSetType(PETSCMat_A, MATMPIAIJ);

        MatSetSizes(PETSCMat_A, nNodes_proc, nNodes_proc, nNodes_wolrd, nNodes_wolrd);
        MatMPIAIJSetPreallocation(PETSCMat_A, d_nz, NULL, o_nz, NULL);
        //Maybe an over estimate of the diag/off diag #of non-zero entries, so we turn off malloc warnings
        MatSetUp(PETSCMat_A);
        MatSetOption(PETSCMat_A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);

        VecCreateMPI(PETSC_COMM_WORLD, nNodes_proc, nNodes_wolrd, &PETSCVec_x);
        VecDuplicate(PETSCVec_x, &PETSCVec_b);
    }

    Mat_A = PETSCMat_A;
    Vec_x = PETSCVec_x;
    Vec_b = PETSCVec_b;
}

void PetscIJ::assembleSystem()
{
    MatAssemblyBegin(PETSCMat_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(PETSCMat_A, MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(PETSCVec_x);
    VecAssemblyEnd(PETSCVec_x);
    //
    VecAssemblyBegin(PETSCVec_b);
    VecAssemblyEnd(PETSCVec_b);
}

void PetscIJ::solverSetupAndSolve()
{
    KSPCreate(PETSC_COMM_WORLD, &PETSCsolver);
    KSPSetOperators(PETSCsolver, PETSCMat_A, PETSCMat_A);

    // set ksp solver type
    // if (ls_solver_ == GMRES)
    // {
    //     KSPSetType(PETSCsolver, KSPGMRES);
    //     KSPGMRESSetOrthogonalization(PETSCsolver, KSPGMRESModifiedGramSchmidtOrthogonalization);
    //     //KSPGMRESModifiedGramSchmidtOrthogonalization();
    // }
    // else if (ls_solver_ == CG)
    // {
    //     KSPSetType(PETSCsolver, KSPCG);
    // }
    // else if (ls_solver_ == BICGS)
    // {
    //     KSPSetType(PETSCsolver, KSPBCGS);
    // }

    // Set up preconditioner
    PC pc;
    KSPGetPC(PETSCsolver, &pc);

    // Classic AMG
    PCSetType(pc, PCGAMG);
    PCGAMGSetType(pc, PCGAMGAGG);
    PCGAMGSetNSmooths(pc, 0);

    // PCSetType(pc, PCJACOBI);

    if (is_singular_)
    {
        MatNullSpace nullspace;
        MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
        MatNullSpaceRemove(nullspace, PETSCVec_b);
        MatSetNullSpace(PETSCMat_A, nullspace);
        MatNullSpaceDestroy(&nullspace);
    }

    KSPSetTolerances(PETSCsolver, rel_tolerance, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    KSPSolve(PETSCsolver, PETSCVec_b, PETSCVec_x);
    // if (verbose >= 2)
    {
        int niters;
        amrex::Real res;
        KSPGetIterationNumber(PETSCsolver, &niters);
        KSPGetResidualNorm(PETSCsolver, &res);
        amrex::PrintToFile("log") << "\n"
                                    << "PETSc Iterations, Residual Norm\t" << niters << "\t" << res << "\n";
    }

    if (print_system)
    {
        /// write mat A
        {
            PetscViewer viewer;
            PetscViewerASCIIOpen(comm_, "MatA.txt", &viewer);
            MatView(PETSCMat_A, viewer);
            PetscViewerDestroy(&viewer);
        }

        /// write vec b
        {
            PetscViewer viewer;
            PetscViewerASCIIOpen(comm_, "Vecb.txt", &viewer);
            VecView(PETSCVec_b, viewer);
            PetscViewerDestroy(&viewer);
        }

        /// write vec x
        {
            PetscViewer viewer;
            PetscViewerASCIIOpen(comm_, "Vecx.txt", &viewer);
            VecView(PETSCVec_x, viewer);
            PetscViewerDestroy(&viewer);
        }
    }
}

} /*End namespace mycode */
