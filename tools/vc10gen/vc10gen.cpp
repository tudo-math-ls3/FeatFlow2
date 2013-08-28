#include <Windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <direct.h>

using namespace std;

// path check file
#define KERNEL_PROJECT_PATH "./kernel/kernel.vfproj"

class VC10Gen
{
public:
  // build ids
  vector<string> build_ids;
  // kernel GUID string
  string kernel_guid;
  // project GUID string
  string project_guid;
  // project name
  string project_name;
  // project path relative to root dir
  string project_path;
  // root path relative to project dir
  string root_path;

  // data file list
  set<string> dat_list;
  // source file list
  set<string> src_list;

public:
  explicit VC10Gen(string _project_path = "", string _project_name = "") :
    kernel_guid(read_kernel_guid()),
    project_guid(gen_random_guid()),
    project_path(_project_path),
    project_name(_project_name)
  {
    // generate build ids
    build_ids.push_back("dbg");
    build_ids.push_back("dbg-omp");
    build_ids.push_back("opt");
    build_ids.push_back("opt-omp");
  }

  static string read_kernel_guid()
  {
    cout << "Checking Kernel Project File '" << KERNEL_PROJECT_PATH << "'..." << endl;
    // try to open the stream
    ifstream ifs(KERNEL_PROJECT_PATH, std::ios_base::in);
    if(ifs.is_open() && ifs.good())
    {
      string line;
      line.reserve(512);
      size_t pos0;
      while(!ifs.eof() && ifs.good())
      {
        getline(ifs, line);
        if((pos0 = line.find("ProjectIdGuid=\"{")) == line.npos)
          continue;
        return line.substr(pos0 + 15, 38);
      }
    }
    throw string("Failed to parse '" + string(KERNEL_PROJECT_PATH) + "'");
  }

  static int gen_random_int()
  {
    static int s((int)time(NULL));
    static int x(362436069);
    static int y(521288629);
    static int z(88675123);
    int t = s ^ (s << 11);
    s = x;
    x = y;
    y = z;
    return z = z ^ (z >> 19) ^ t ^ (t >> 8);
  }

  // generates a random GUID
  static string gen_random_guid()
  {
    // generate 128 random bits
    int b0(gen_random_int());
    int b1(gen_random_int());
    int b2(gen_random_int());
    int b3(gen_random_int());

    // modify b1 and b2
    b1 = (b1 & ~0xF000) | 0x4000;
    b2 = (b2 & ~0x4000) | 0x8000;

    // print to stream
    ostringstream oss;
    oss << setfill('0') << hex << uppercase;
    oss << "{";
    oss << setw(8) << b0 << "-";
    oss << setw(4) << ((b1 >> 16) & 0xFFFF) << "-";
    oss << setw(4) << (b1 & 0xFFFF) << "-";
    oss << setw(4) << (b2 & 0xFFFF) << "-";
    oss << setw(4) << ((b2 >> 16) & 0xFFFF);
    oss << setw(8) << b3;
    oss << "}";

    return oss.str();
  }

  bool check_project_dir(int& depth)
  {
    vector<string> dirs;
    depth = 0u;

    // separate paths
    size_t n0(0);
    while(n0 != project_path.npos)
    {
      // find separators
      size_t n1 = project_path.find_first_of('\\', n0);
      size_t n2 = project_path.find_first_of('/', n0);
      size_t n3(0);
      if((n1 == project_path.npos) && (n2 == project_path.npos))
      {
        dirs.push_back(project_path.substr(n0));
        break;
      }
      else if(n1 == project_path.npos)
        n3 = n2;
      else if(n2 == project_path.npos)
        n3 = n1;
      else
        n3 = min(n1, n2);

      // separate strings
      dirs.push_back(project_path.substr(n0, n3-n0));
      n0 = n3 + 1u;
    }

    // eliminate empty paths
    for(size_t i(0); i < dirs.size(); )
    {
      if(dirs[i].empty())
        dirs.erase(dirs.begin() + i);
      else
        ++i;
    }

    depth = int(dirs.size());

    // create directories
    string path(".");
    for(size_t i(0); i < dirs.size(); ++i)
    {
      path += "/" + dirs[i];

      // check whether the directory exists
      DWORD attribs = GetFileAttributes(path.c_str());
      if((attribs == INVALID_FILE_ATTRIBUTES) || ((attribs & FILE_ATTRIBUTE_DIRECTORY) == 0))
      {
        cout << "ERROR: Application path '" << path << "' does not exist!" << endl;
        return false;
      }
    }

    // okay
    return true;
  }

  bool find_sources(bool add_files = false)
  {
    WIN32_FIND_DATA find_data;
    HANDLE find_handle;

    // find cpp files
    memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
    find_handle = FindFirstFile(string("./" + project_path + "/src/*.f90").c_str(), &find_data);
    if(find_handle !=  INVALID_HANDLE_VALUE)
    {
      bool found(true);
      while(found)
      {
        src_list.insert(find_data.cFileName);
        found = FindNextFile(find_handle, &find_data) != FALSE;
      }
      FindClose(find_handle);
    }

    if(src_list.empty())
      return false;

    cout << endl << "The following source files have been found:" << endl;
    set<string>::iterator it(src_list.begin()), jt(src_list.end());
    for(; it != jt; ++it)
      cout << "- src/" << *it << endl;
    cout << endl;

    if(add_files)
      return true;

    string cmd;
    cout << "Do you want these files to be included in the project file?" << endl;
    cout << "Type 'y' for yes or 'n' for no" << endl;
    cout << ">";
    cin >> cmd;
    cout << endl;

    return (cmd == "y") || (cmd == "yes");
  }

  bool find_data_files(bool add_files = false)
  {
    WIN32_FIND_DATA find_data;
    HANDLE find_handle;

    // find cpp files
    memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
    find_handle = FindFirstFile(string("./" + project_path + "/data/*.dat").c_str(), &find_data);
    if(find_handle !=  INVALID_HANDLE_VALUE)
    {
      bool found(true);
      while(found)
      {
        dat_list.insert(find_data.cFileName);
        found = FindNextFile(find_handle, &find_data) != FALSE;
      }
      FindClose(find_handle);
    }

    if(dat_list.empty())
      return false;

    cout << endl << "The following data files have been found:" << endl;
    set<string>::iterator it(dat_list.begin()), jt(dat_list.end());
    for(; it != jt; ++it)
      cout << "- data/" << *it << endl;
    cout << endl;

    if(add_files)
      return true;

    string cmd;
    cout << "Do you want these files to be included in the project file?" << endl;
    cout << "Type 'y' for yes or 'n' for no" << endl;
    cout << ">";
    cin >> cmd;
    cout << endl;

    return (cmd == "y") || (cmd == "yes");
  }

  bool write_solution()
  {
    // build solution filename
    string sln_path(project_path + "\\" + project_name + ".sln");
    cout << "Writing solution file '" << sln_path << "'..." << endl;
    ofstream ofs(sln_path, ios_base::out|ios_base::trunc);
    if(!ofs.is_open() || !ofs.good())
    {
      cout << "ERROR: Failed to create '" + sln_path + "'" << endl;
      return false;
    }

    // write utf-8 bom
    ofs << char(-17) << char(-69) << char(-65) << endl;

    // write header
    ofs << "Microsoft Visual Studio Solution File, Format Version 11.00" << endl;
    ofs << "# Visual Studio 2010" << endl;

    // include app project
    ofs << "Project(\"{6989167D-11E4-40FE-8C1A-2192A86A7E90}\") = \"" << project_name << "\", \"";
    ofs << project_name << ".vfproj\", \"" << project_guid << "\"" << endl;
    ofs << "\tProjectSection(ProjectDependencies) = postProject" << endl;
    ofs << "\t\t" << kernel_guid << " = " << kernel_guid << endl;
    ofs << "\tEndProjectSection" << endl;
    ofs << "EndProject" << endl;

    // include kernel project
    ofs << "Project(\"{6989167D-11E4-40FE-8C1A-2192A86A7E90}\") = \"kernel\", \"";
    ofs << root_path << "\\kernel\\kernel.vfproj\", \"" << kernel_guid << "\"" << endl;
    ofs << "EndProject" << endl;

    // write solution configurations
    ofs << "Global" << endl;
    ofs << "\tGlobalSection(SolutionConfigurationPlatforms) = preSolution" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << *it << "|Win32 = " << *it << "|Win32" << endl;
      ofs << "\t\t" << *it << "|x64 = " << *it << "|x64" << endl;
    }
    ofs << "\tEndGlobalSection" << endl;

    // write project configurations
    ofs << "\tGlobalSection(ProjectConfigurationPlatforms) = postSolution" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << project_guid << "." << *it << "|Win32.ActiveCfg = " << *it << "|Win32" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|Win32.Build.0 = " << *it << "|Win32" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x64.ActiveCfg = " << *it << "|x64" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x64.Build.0 = " << *it << "|x64" << endl;
    }
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << kernel_guid << "." << *it << "|Win32.ActiveCfg = " << *it << "|Win32" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|Win32.Build.0 = " << *it << "|Win32" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x64.ActiveCfg = " << *it << "|x64" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x64.Build.0 = " << *it << "|x64" << endl;
    }
    ofs << "\tEndGlobalSection" << endl;

    // write solution properties
    ofs << "\tGlobalSection(SolutionProperties) = preSolution" << endl;
    ofs << "\t\tHideSolutionNode = FALSE" << endl;
    ofs << "\tEndGlobalSection" << endl;

    // end-of-file
    ofs << "EndGlobal" << endl;
    ofs.close();

    return true;
  }

  void write_project_config(std::ofstream& ofs, const string bid, const string platform)
  {
    bool dbg = (bid.find("dbg") != bid.npos);
    bool omp = (bid.find("omp") != bid.npos);
    bool x64 = (platform.compare("x64") == 0);

    const string bid2 = string("if12-") + (x64 ? "x64" :  "x86") + "-" + bid;

    ofs << "\t\t<Configuration Name=\"" << bid << "|" << platform << "\" OutputDirectory=\"$(ProjectDir)\" ";
    ofs << "IntermediateDirectory=\"$(ProjectDir)\\obj\\$(ProjectName)." << bid2 << "\">" << endl;

    // write compiler options
    ofs << "\t\t\t\t<Tool Name=\"VFFortranCompilerTool\" SuppressStartupBanner=\"true\"";
    // debug info: line info only for 'opt' mode, otherwise full debug info
    ofs << " DebugInformationFormat=\"" << (dbg ? "debugEnabled" : "debugLineInfoOnly") << "\"";
    // optimisation: full for 'opt' mode, otherwise disabled
    ofs << " Optimization=\"optimize" << (dbg ? "Disabled" : "Full") << "\"";
    // auto-parallelisation: only for 'opt-omp' mode
    ofs << ((!dbg && omp) ? " Parallelization=\"true\"" : "");
    // i/o buffering: enable for 'opt' mode, disable for 'dbg' mode
    ofs << (!dbg ? " BufferedIO=\"true\"" : "");
    // prepreocessor stuff
    ofs << " Preprocess=\"preprocessYes\"";
    ofs << " AdditionalIncludeDirectories=\"$(ProjectDir)" << root_path << "\"";
    // openmp: enable for 'omp' modes, disable otherwise
    ofs << (omp ? " OpenMP=\"OpenMPParallelCode\"" : " OpenMPConditionalCompilation=\"false\"");
    // check routine interfaces in 'dbg' mode
    ofs << (dbg ? " WarnInterfaces=\"true\"" : "");
    // fp exception: throw exceptions in 'dbg' mode, produce NaNs and INFs in 'opt' mode
    ofs << " FloatingPointExceptionHandling=\"" << (dbg ? "fpe0" : "fpe1") << "\"";
    // fp model: fast for 'opt' mode; 'source' (i.e. precise) for 'dbg' mode
    ofs << (dbg ? " FloatingPointModel=\"source\"" : "");
    // always flush denormals to zero
    ofs << " FlushDenormalResultsToZero=\"true\"";
    // always use C calling convention
    ofs << " CallingConvention=\"callConventionCRef\"";
    // path for debug database file
    ofs << " PdbFile=\"$(OutDir)\\$(TargetName).pdb\"";
    if(dbg)
    {
      // various other checks only for debug mode
      ofs << " Traceback=\"true\" NullPointerCheck=\"true\" BoundsCheck=\"true\"";
      ofs << " UninitializedVariablesCheck=\"true\" ArgTempCreatedCheck=\"true\"";
      ofs << " RuntimeLibrary=\"rtMultiThreadedDebug\" Interfaces=\"true\"";
    }
    ofs << "/>" << endl;

    // write linker options
    ofs << "\t\t\t\t<Tool Name=\"VFLinkerTool\" OutputFile=\"$(OutDir)\\$(ProjectName)." << bid2 << ".exe\"";
    // no incremental linking
    ofs << " LinkIncremental=\"linkIncrementalNo\" SuppressStartupBanner=\"true\"";
    // no manifest
    ofs << " GenerateManifest=\"false\"";
    // always generate debug info
    ofs << " GenerateDebugInformation=\"true\"";
    // we're always linking a console app
    ofs << " SubSystem=\"subSystemConsole\"";
    // we never use ipo...
    ofs <<" InterproceduralOptimizations=\"false\"/>" << endl;

    // write the rest
    ofs << "\t\t\t\t<Tool Name=\"VFResourceCompilerTool\"/>" << endl;
    ofs << "\t\t\t\t<Tool Name=\"VFMidlTool\" SuppressStartupBanner=\"true\""
      << (x64 ? " TargetEnvironment=\"midlTargetAMD64\"" : "") << "/>" << endl;
    ofs << "\t\t\t\t<Tool Name=\"VFCustomBuildTool\"/>" << endl;
    ofs << "\t\t\t\t<Tool Name=\"VFPreLinkEventTool\"/>" << endl;
    ofs << "\t\t\t\t<Tool Name=\"VFPreBuildEventTool\"/>" << endl;
    ofs << "\t\t\t\t<Tool Name=\"VFPostBuildEventTool\"/>" << endl;
    ofs << "\t\t\t\t<Tool Name=\"VFManifestTool\" SuppressStartupBanner=\"true\"/></Configuration>";
    // Note: The missing '<< endl' in the last line is intended.
  }

  bool write_project()
  {
    // build project filename
    string prj_path(project_path + "\\" + project_name + ".vfproj");
    cout << "Writing project file  '" << prj_path << "'..." << endl;
    ofstream ofs(prj_path, ios_base::out|ios_base::trunc);
    if(!ofs.is_open() || !ofs.good())
    {
      cout << "ERROR: Failed to create '" + prj_path + "'" << endl;
      return false;
    }

    // write file header
    ofs << "<\?xml version=\"1.0\" encoding=\"UTF-8\"\?>" << endl;
    ofs << "<VisualStudioProject ProjectCreator=\"Intel Fortran\" "
      << "Keyword=\"Console Application\" Version=\"11.0\" "
      << "ProjectIdGuid=\"" << project_guid << "\">" << endl;

    // write platforms
    ofs << "\t<Platforms>" << endl;
    ofs << "\t\t<Platform Name=\"Win32\"/>" << endl;
    ofs << "\t\t<Platform Name=\"x64\"/></Platforms>" << endl;

    // write configurations
    ofs << "\t<Configurations>";
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it)
    {
      ofs << endl;
      write_project_config(ofs, *it, "Win32");
      ofs << endl;
      write_project_config(ofs, *it, "x64");
    }
    ofs << "</Configurations>" << endl;

    // start file list
    ofs << "\t<Files>";

    // write data file list
    if(!dat_list.empty())
    {
      ofs << endl << "\t\t<Filter Name=\"Data Files\" Filter=\"dat\">";
      for(set<string>::iterator it(dat_list.begin()); it != dat_list.end(); ++it)
        ofs << endl << "\t\t<File RelativePath=\".\\data\\" << *it << "\"/>";
      ofs << "</Filter>";
    }

    // write source file list
    if(!src_list.empty())
    {
      ofs << endl << "\t\t<Filter Name=\"Source Files\" Filter=\"f90;for;f;fpp;ftn;def;odl;idl\">";
      for(set<string>::iterator it(src_list.begin()); it != src_list.end(); ++it)
        ofs << endl << "\t\t<File RelativePath=\".\\src\\" << *it << "\"/>";
      ofs << "</Filter>";
    }

    // end file list
    ofs << "</Files>" << endl;

    // end-of-file
    ofs << "\t<Globals/></VisualStudioProject>" << endl;
    ofs.close();

    return true;
  }

  void main(bool add_files)
  {
    // write kernel and project GUID
    cout << "Kernel  GUID: " << kernel_guid << endl;
    cout << "Project GUID: " << project_guid << endl;

    // read project path
    if(project_path.empty())
    {
      cout << endl << "Please enter the path to the project directory:" << endl;
      cout << ">";
      cin >> project_path;
      if(!cin.good())
        return;
      cout << endl;
    }

    // read project name
    if(project_name.empty())
    {
      cout << "Please enter the project name:" << endl << ">";
      cin >> project_name;
      if(!cin.good())
        return;
      cout << endl;
    }

    // generate project directories
    int path_depth(0);
    if(!check_project_dir(path_depth))
      return;

    // build relative kernel project path
    root_path = "..";
    for(int i(1); i < path_depth; ++i)
      root_path += "\\..";

    // build header and source lists
    if(!find_sources(add_files))
      src_list.clear();
    if(!find_data_files(add_files))
      dat_list.clear();

    // write solution file
    if(!write_solution())
      return;

    // write project file
    if(!write_project())
      return;

    cout << endl << "Project files for '" << project_name << "' have been written successfully" << endl;
  }
};


// application main entrypoint
int main(int argc, char** argv)
{
  // print header
  cout << endl << "FeatFlow2 Visual Studio 2010 Project File Generator" << endl << endl;

  // print warning
  cout << "WARNING: This tool is still under construction, so use with care!" << endl << endl;

  string prj_name, prj_path;
  if(argc > 1)
  {
    prj_path = argv[1];
    if(argc > 2)
      prj_name = argv[2];
    else
    {
      size_t pos1 = prj_path.find_last_of('\\');
      size_t pos2 = prj_path.find_last_of('/');
      size_t pos = prj_path.npos;
      if((pos1 != pos) && (pos2 != pos))
        pos = (pos1 < pos ? pos1 : pos2);
      else if(pos1 != pos)
        pos = pos1;
      else if(pos2 != pos)
        pos = pos2;
      else
      {
        cout << "ERROR: Could not extract project name from '" << prj_path << "' !" << endl << endl;
        return 1;
      }
    }
  }

  try
  {
    VC10Gen generator(prj_path, prj_name);
    generator.main(argc > 1);
  }
  catch(string& s)
  {
    cout << "ERROR: " << s << endl;
  }
  return 0;
}
