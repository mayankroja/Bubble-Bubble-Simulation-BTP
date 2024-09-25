#include "incFSI.H"
#include <CFMask.H>
#include <AMReX_PlotFileUtil.H>

namespace mycode
{

    //!Restart
    void incFSI::WriteCheckpointFile(const std::string &checkpointname)
    {
        // chk00010            write a checkpoint file with this root directory
        // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
        //                     the BoxArrays at each level
        // chk00010/Level_0/
        // chk00010/Level_1/
        // etc.                these subdirectories will hold the MultiFab data at each level of refinement
        // checkpoint file name, e.g., chk00010
        
        amrex::Print() << "writing checkpoint file " << checkpointname << "\n";
        amrex::PrintToFile("log") << "writing checkpoint file " << checkpointname << "\n";
        const int nlevels = finest_level + 1;
        // ---- prebuild a hierarchy of directories
        // ---- dirName is built first.  if dirName exists, it is renamed.  then build
        // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
        // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
        // ---- after all directories are built
        // ---- ParallelDescriptor::IOProcessor() creates the directories
        amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);
        // write Header file
        if (amrex::ParallelDescriptor::IOProcessor())
        {
    
            std::string HeaderFileName(checkpointname + "/Header");
            amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out |
                                                        std::ofstream::trunc |
                                                        std::ofstream::binary);
            if (!HeaderFile.good())
            {
                amrex::FileOpenFailed(HeaderFileName);
            }
    
            HeaderFile.precision(17);
    
            // write out title line
            HeaderFile << "Checkpoint file for incFSI\n";
    
            // write out finest_level
            HeaderFile << finest_level << "\n";
    
            // write out Iter
            HeaderFile << Iter;
            HeaderFile << "\n";
    
            // write out dt
            //! note that dt is corresponding to next step
            HeaderFile << dt;
            HeaderFile << "\n";
    
            // write out time
            //! note that time is corresponding to Iter
            HeaderFile << Time;
            HeaderFile << "\n";
    
            // write the BoxArray
            for (int lev = 0; lev <= finest_level; ++lev) 
            {
                //boxArray(lev).writeOn(HeaderFile);
                grids[lev].writeOn(HeaderFile);
                HeaderFile << '\n';
            }

    
            //! for coupled problem forces on interfaces will be required
            // no. of interfaces
            if(N_IF > 0)
            {
                for (int lev = 0; lev <= finest_level; ++lev)
                {
                    HeaderFile << interfaces[lev].size();
                    HeaderFile << '\n';
                    //amrex::Print(-1)<<interfaces[lev].size()<<'\n';
                    //! write vol and P_interface for advecting interfaces
                    {
                        int nsolid = 0;
                        for (auto &&solid : interfaces[lev])
                        {
                            if (solid->isAdvectLS())
                            {
                                HeaderFile << solid->Volume_0() << " " << solid->Volume() << " " << solid->getP_interface() << " " << solid->getP_interface0();
                                nsolid++;
                            }
                            HeaderFile << '\n';
                        }
                    }
                }
            }
        }

        //! write data
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            amrex::VisMF::Write(xvel[lev], amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "u"));
            amrex::VisMF::Write(yvel[lev], amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "v"));
            amrex::VisMF::Write(Pressure[lev], amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "P"));
            if(TempField) amrex::VisMF::Write(Theta[lev], amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "T"));
            if(PhaseField) amrex::VisMF::Write(Phi[lev], amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Phi"));
	    if(DamageModel)
            {
	        for(int iscalar = 0; iscalar < nscalar; iscalar++)
		{
	            const std::string& scalarfile = amrex::Concatenate("Scalar",iscalar,2);
		    amrex::VisMF::Write(Scalars[lev][iscalar], amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", scalarfile));
		}
	    }
        }
        if(N_IF > 0)
        {
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                for (auto &&solid : interfaces[lev])
                {
                    if (solid->isAdvectLS())
                    {
                        std::string var_name = solid->name() + "_Psi";
                        amrex::VisMF::Write(solid->Psi(), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", var_name));
                    }
                }
            }
        }
    }
    
    namespace
    {
        // utility to skip to next line in Header
        void GotoNextLine(std::istream &is)
        {
            constexpr std::streamsize bl_ignore_max{100000};
            is.ignore(bl_ignore_max, '\n');
        }
    
        void GetCheckFileName(std::string &check_file)
        {
            //! check status in InputR/status.dat and Input/status.dat
            int stat = 0, statR = 0;
            amrex::Real T = 0.0, TR = 0.0;
            // Input/status
            {
                std::string File("Input/status.dat");
                std::string word;
                std::ifstream dataFile(File);
    
                std::string lineContents;
                while (getline(dataFile, lineContents))
                {
                    std::stringstream ss(lineContents);
                    ss >> word; // time
                    T = std::stod(word);
                    ss >> word; // status
                    stat = std::stoi(word);
                }
                dataFile.close();
            }
            // InputR/status
            {
                std::string File("InputR/status.dat");
                std::string word;
                std::ifstream dataFile(File);
    
                std::string lineContents;
                while (getline(dataFile, lineContents))
                {
                    std::stringstream ss(lineContents);
                    ss >> word; // time
                    TR = std::stod(word);
                    ss >> word; // status
                    statR = std::stoi(word);
                }
                dataFile.close();
            }
    
            if (stat == 1 && statR == 1)
            {
                if (T > TR)
                {
                    check_file = "Input/chk";
                }
                else
                {
                    check_file = "InputR/chk";
                }
            }
            else if (stat == 1)
            {
                check_file = "Input/chk";
            }
            else if (statR == 1)
            {
                check_file = "InputR/chk";
            }
            else
            {
                amrex::Abort("chk files are currupted.. can not restart");
            }
        }
    
    } // namespace
    
    void incFSI::ReadCheckpointFile()
    {
        //GetCheckFileName(chk_file);
        chk_file = amrex::Concatenate("Restart/chk",chk_int_read,5);
        amrex::PrintToFile("log") << "Restart from checkpoint " << chk_file << "\n";
    
        // Header
        std::string File(chk_file + "/Header");
    
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());
    
        amrex::Vector<char> fileCharPtr;
        amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);
    
        std::string line, word;
    
        // read in title line
        std::getline(is, line);
    
        // read in finest_level
        //int finest_level;
        is >> finest_level;
        GotoNextLine(is);
    
        // read in istep
        is >> Iter;
    
        // read in dt
        is >> dt;
    
        // read in time
        is >> Time;
        
        for (int lev = 0; lev <= finest_level; ++lev) {
        
            // read in level 'lev' BoxArray from Header
            amrex::BoxArray ba;
            ba.readFrom(is);
            GotoNextLine(is);
        
            // create a distribution mapping
            amrex::DistributionMapping dm { ba, amrex::ParallelDescriptor::NProcs() };
        
            // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
            SetBoxArray(lev, ba);
            SetDistributionMap(lev, dm);
        
            AllocateMemory(lev, ba, dm);
            if(N_IF > 0) MakeInterface(lev, ba, dm);
        }
        
        //! note the box array will be read from mesh in the main.cpp
        //! here this ba is unused
        if(N_IF > 0)
        { 
            //! Read vol and P_interface for advecting interfaces
            for (int lev = 0; lev <= finest_level; ++lev) 
            {
                int n;
                is >> n;
                GotoNextLine(is);

                int nsolid = 0;
                for (auto &&solid : interfaces[lev])
                {
                    if (solid->isAdvectLS())
                    {
                        std::getline(is, line);
                        {
                            std::istringstream lis(line);
                            lis >> word;
                            solid->setVolume0(std::stod(word));
                            lis >> word;
                            solid->setVolume(std::stod(word));
                            lis >> word;
                            //P_interface0[nsolid] = std::stod(word);//Read from input in MakeInterface
                            lis >> word;
                            solid->setP_interface(std::stod(word));
                        }
                        nsolid++;
                    }
                }
            }
        }  

        //! read u, v, and Pressure
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            xvel[lev].setVal(0.0);
            yvel[lev].setVal(0.0);
            xvel_old[lev].setVal(0.0);
            yvel_old[lev].setVal(0.0);
            U[lev].setVal(0.0);
            Pressure[lev].setVal(0.0);
            Pstar[lev].setVal(0.0);

            amrex::VisMF::Read(xvel[lev], amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", "u"));
            amrex::VisMF::Read(yvel[lev], amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", "v"));
            amrex::VisMF::Read(Pressure[lev], amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", "P"));
            if(TempField) amrex::VisMF::Read(Theta[lev], amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", "T"));
            if(PhaseField) amrex::VisMF::Read(Phi[lev], amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", "Phi"));
	    if(DamageModel)
	    {
		for(int iscalar = 0; iscalar < nscalar ; iscalar++)
		{
		    const std::string& scalarfile = amrex::Concatenate("Scalar",iscalar,2);
	            amrex::VisMF::Read(Scalars[lev][iscalar], amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", scalarfile));
	        }
	    }
        }

        if(N_IF > 0)
        {
            for (int lev = 0; lev <= finest_level; ++lev)  
            {
        //! if fixed interface then in scheme these are not called so needs to be called here
                //! read Psi
                for (auto &&solid : interfaces[lev])
                {
                    if (solid->isAdvectLS())
                    {
                        std::string var_name = solid->name() + "_Psi";
                        amrex::VisMF::Read(solid->Psi(), amrex::MultiFabFileFullPrefix(lev, chk_file, "Level_", var_name));
                    }
                    mask_[lev]->GhostCellIdentfication();
		    if(lev == finest_level)
	            {
                        ComputeCutFaceVel(lev);
                        ComputeCutFacePressure(lev);
		    }
                }
            }
        }
    }
    
    void incFSI::ChkFile()
    {
        const std::string& dir = amrex::Concatenate("Restart/chk",Iter,5);
        //std::string status_file = dir + "/status.dat";
        std::string status_file = "status.dat";
        //amrex::Print()<<"status_file = "<<status_file<<"\n"; 
        //! status of chk file before writing is 0
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            //! write out the status of chk file
            std::ofstream ofs(status_file, std::ofstream::out);
            ofs << Time << "\t0\n";
            ofs.close();
        }
    
        WriteCheckpointFile(dir);
        //! old chk files will be renamed to new directories named like chk00350.old.46576787980.
        //! need to remove these old files
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::string cmd = "rm -r  Restart/*.old.* 2> /dev/null";
            system(cmd.c_str());
    
            //! status of chk file after writing is 1
            std::ofstream ofs(status_file, std::ofstream::out);
            ofs << Time << "\t1\n";
            ofs.close();
        }
    }

} // namespace mycode
