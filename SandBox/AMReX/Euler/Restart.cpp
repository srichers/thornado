/*
Restart.cpp

Originally created by nataraj for one MultiFab

Modified by sjdunham for multiple MultiFabs (not working yet)

*/

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <Restart.H>

using namespace amrex;

extern "C"
{
  void writecheckpointfile
         ( int StepNo[], int FinestLevel, Real dt[], Real time[],
           BoxArray** BA, MultiFab** MF_uGF, MultiFab** MF_uCF,
           MultiFab** MF_uPF, MultiFab** MF_uAF )
  {

    // chk00010            Write a checkpoint file with this root directory
    // chk00010/Header     This contains information you need to save
    //                       (e.g., FinestLevel, time, etc.), and also
    //                       the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                These subdirectories will hold the MultiFab
    //                     data at each level of refinement

    // Checkpoint file name, e.g., chk00010

    const std::string& checkpointname
                         = amrex::Concatenate( chk_file, StepNo[0] );

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nLevels = FinestLevel+1;

    // ---- Prebuild a hierarchy of directories
    // ---- dirName is built first. If dirName exists, it is renamed. Then build
    // ----   dirName/subDirPrefix_0 .. dirName/subDirPrefix_nLevels-1
    // ---- If callBarrier is true, call ParallelDescriptor::Barrier()
    // ----   after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy( checkpointname, "Level_", nLevels, true );

    // Write Header file
    if ( ParallelDescriptor::IOProcessor() )
    {

      std::string HeaderFileName( checkpointname + "/Header" );
      std::ofstream HeaderFile( HeaderFileName.c_str(), std::ofstream::out   |
				                        std::ofstream::trunc |
				                        std::ofstream::binary );
      if( ! HeaderFile.good() )
      {
        amrex::FileOpenFailed( HeaderFileName );
      }

      HeaderFile.precision(17);

      VisMF::IO_Buffer io_buffer( VisMF::IO_Buffer_Size );
      HeaderFile.rdbuf() -> pubsetbuf( io_buffer.dataPtr(), io_buffer.size() );

      // Write out title line
      HeaderFile << "Checkpoint file for IsentropicVortex\n";

      // Write out FinestLevel
      HeaderFile << FinestLevel << "\n";

      // Write out array of StepNo
      for( int i = 0; i < nLevels; ++i )
      {
        HeaderFile << StepNo[i] << " ";
      }
      HeaderFile << "\n";

      // Write out array of dt
      for(int i = 0; i < nLevels; ++i)
      {
        HeaderFile << dt[i] << " ";
      }
      HeaderFile << "\n";

      // Write out array of time
      for( int i = 0; i < nLevels; ++i )
      {
        HeaderFile << time[i] << " ";
      }
      HeaderFile << "\n";

      // Write the BoxArray at each level
      for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
      {
        BoxArray& BA1 = *BA[iLevel];
        BA1.writeOn(HeaderFile);
        HeaderFile << '\n';
      }

    } // End of Header file writing

    // Write the MultiFab data to, e.g., chk00010/Level_0/
    for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
    {
      MultiFab& MF_uGF1 = *MF_uGF[iLevel];
      MultiFab& MF_uCF1 = *MF_uCF[iLevel];
      MultiFab& MF_uPF1 = *MF_uPF[iLevel];
      MultiFab& MF_uAF1 = *MF_uAF[iLevel];
      VisMF::Write( MF_uGF1, amrex::MultiFabFileFullPrefix
                                      ( iLevel, checkpointname,
                                        "Level_", "Geometry" ) );
      VisMF::Write( MF_uCF1, amrex::MultiFabFileFullPrefix
                                      ( iLevel, checkpointname,
                                        "Level_", "Conserved" ) );
      VisMF::Write( MF_uPF1, amrex::MultiFabFileFullPrefix
                                      ( iLevel, checkpointname,
                                        "Level_", "Primitive" ) );
      VisMF::Write( MF_uAF1, amrex::MultiFabFileFullPrefix
                                      ( iLevel, checkpointname,
                                        "Level_", "Auxiliary" ) );
    }

  } // End of WriteCheckpointFile function

  void readheaderandboxarraydata
         ( int finest_level[], int stepno[], Real dt[], Real time[],
           BoxArray** ba, DistributionMapping** dm )
  {

    ParmParse pp( "amr" );
    pp.query( "restart", restart_chkfile );

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File( restart_chkfile + "/Header" );

    VisMF::IO_Buffer io_buffer( VisMF::GetIOBufferSize() );

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile( File, fileCharPtr );
    std::string fileCharPtrString( fileCharPtr.dataPtr() );
    std::istringstream is( fileCharPtrString, std::istringstream::in );

    std::string line, word;

    // Read in title line
    std::getline( is, line );

    // Read in finest_level
    is >> finest_level[0];
    GotoNextLine( is );

    // Read in array of istep
    std::getline( is, line );
    {
      std::istringstream lis( line );
      int i = 0;
      while( lis >> word )
      {
        stepno[i++] = std::stoi( word );
      }
    }

    // Read in array of dt
    std::getline( is, line );
    {
      std::istringstream lis( line );
      int i = 0;
      while( lis >> word )
      {
        dt[i++] = std::stod( word );
      }
    }

    // Read in array of time
    std::getline( is, line );
    {
      std::istringstream lis( line );
      int i = 0;
      while( lis >> word )
      {
        time[i++] = std::stod( word );
      }
    }

    // Read in level 'lev' BoxArray from Header
    for( int lev = 0; lev <= finest_level[0]; ++lev )
    {
      BoxArray& ba1=*ba[lev];
      ba1.readFrom( is );
      ba[lev] = &ba1;
      GotoNextLine( is );

      // Create a distribution mapping
      DistributionMapping dm1{ ba1, ParallelDescriptor::NProcs() };
      *dm[lev]=dm1;
    }

  } // End of readheaderandboxarraydata function

  void readmultifabdata( int FinestLevel, MultiFab** MF, int iMF) 
  {
    ParmParse pp( "amr" );
    pp.query( "restart", restart_chkfile );

    // Header
    std::string File( restart_chkfile + "/Header" );

    VisMF::IO_Buffer io_buffer( VisMF::GetIOBufferSize() );

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile( File, fileCharPtr );
    std::string fileCharPtrString( fileCharPtr.dataPtr() );
    std::istringstream is( fileCharPtrString, std::istringstream::in );

    // Read in the MultiFab data
    std::string MF_Name;
    switch( iMF )
    {
      case 0:
        MF_Name = "Geometry";
	break;
      case 1:
	MF_Name = "Conserved";
	break;
      case 2:
	MF_Name = "Primitive";
	break;
      case 3:
	MF_Name = "Auxiliary";
	break;
      default:
        std::cout << "Invalid." << std::endl;
    }

    for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
    {
      VisMF::Read( *MF[iLevel], 
        amrex::MultiFabFileFullPrefix
                 ( iLevel, restart_chkfile, "Level_", MF_Name ) );
    }

  } // End of readmultifabdata function

  void GotoNextLine ( std::istream& is )
  {
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore( bl_ignore_max, '\n' );
  }

} // End of extern "C" block
