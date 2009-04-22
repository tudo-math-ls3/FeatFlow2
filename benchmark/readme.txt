The basic regression test suite -- short introduction for dummys
----------------------------------------------------------------
The FEAT regression test suite allows to do a couple of benchmarks
to the code every night to ensure everything compiles and runs
correctly. The basic concept is that a couple of test applications
run on a couple of machines producing some output which is compared
to reference values. EMails are generated when some tests fail so that
the user responsible for an improper checkin can be blamed.

General structure of the benchmarks - Quick Start Guide
-----------------------------------------------------------
a) Benchmark applications are compiled in the directory of the
   benchmark application, e.g.
     "Featflow2/benchmark/apps_cc2d". 
    
b) The actual applications on the other hand are always started
   in the main benchmark directory, i.e.
     "Featflow2/benchmark". 
     
c) The subdirectory
     "Featflow2/benchmark/data"
   contains a couple of subdirectories for every application containing
   the data files of the corresponding test(s). The data files for the
   standard cc2d test can be found e.g. in the directory
     "Featflow2/benchmark/data/apps_cc2d"
   This is in slight contrast to the standard directory structure
   of a usual application where the "data/" subdirectory contains
   the actual data.

d) For every benchmark application, the directory
     "Featflow2/benchmark/tests"
   contains a file that configures the test definitions that should
   be performed with an application. The test definitions for the
   standard cc2d application can be found e.g. in the file
      "Featflow2/benchmark/tests/cc2d.fbdef"
   Each test starts with a parameter "testid=..." that defines the
   ID of the test. The parameter "appl =..." defines the name of the
   application. It should coincide with the name of the .fbconf file.

e) When an application is started, the following data is passed to it:

   - The name/path of the data file specified by the "datfile=..." 
     parameter in the sections of the corresponding .fbdef file
     is passed as 1st command line parameter.

     For the cc2d application with test-id CC2D_001 e.g., the "datfile"
     parameter is defined as
       "datfile  = ./data/apps_cc2d/master.dat"
     which points to the master file that controls the application
     relative to the benchmark directory.
     The application is then executed by calling
       "feat2Benchmark-cc2d ./data/apps_cc2d/master.dat"
       
   - The parameter "mglevels = x,x,x,..." defines a list of one or
     multiple levels. For every value in the list, a separate benchmark
     is started. The current level is then passed as environment
     variable
       $MGLEVEL = ...
     to the application.
     
   - Below each test-id, the user can specify an arbitrary number
     of variables which are passed as environment variables to the
     application.
     
     Example for cc2d.conf:
       testid   = CC2D_001
       datfile  = ./data/apps_cc2d/master.dat
       mglevels = 4,5
       element  = 11
       eps      = 1E-10
       
     Defines test case "CC2D_001". The benchmark defines the environment
     variables
       $ELEMENT = 11
       $EPS     = 1E-10
     and calls the application two times with the environment variable
     $MGLEVEL set to $MGLEVEL=4 and $MGLEVEL=5.
     
   - Parameters defined in one test id are inherited by the next one
     until redefinition.
     
     Example for cc2d.conf:
     
       appl = cc2d
     
       testid   = CC2D_001
       datfile  = ./data/apps_cc2d/master.dat
       mglevels = 4,5
       element  = Q1T
       eps      = 1E-10

       testid   = CC2D_002
       element  = Q2

       testid   = CC2D_003
       element  = Q1T
       eps      = 1E-5

       testid   = CC2D_004
       element  = Q2
       
     This defines the following test cases:
     * CC2D_001, cc2d, $MGLEVEL=4, $ELEMENT=Q1T, $EPS=1E-10
     * CC2D_001, cc2d, $MGLEVEL=5, $ELEMENT=Q1T, $EPS=1E-10
     * CC2D_002, cc2d, $MGLEVEL=4, $ELEMENT=Q2, $EPS=1E-10
     * CC2D_002, cc2d, $MGLEVEL=5, $ELEMENT=Q2, $EPS=1E-10
     * CC2D_003, cc2d, $MGLEVEL=4, $ELEMENT=Q1T, $EPS=1E-5
     * CC2D_003, cc2d, $MGLEVEL=5, $ELEMENT=Q1T, $EPS=1E-5
     * CC2D_004, cc2d, $MGLEVEL=4, $ELEMENT=Q2, $EPS=1E-5
     * CC2D_004, cc2d, $MGLEVEL=5, $ELEMENT=Q2, $EPS=1E-5

f) During the execution, the following variables in the environment
   are automatically set:
     $LOGDIR = Directory, where the application must write log files to.
   The application must write a file named
     "$LOGDIR/benchmarkresultfile"
   to this directory that contains deterministic benchmark results
   that should be compared to reference results. Timing results
   and similar nondeterministic results should not be written to this
   file as the file is later on diff'ed to a file with reference
   results.
   
g) SUMMARY: The data flow for executing a benchmark application is
   as follows:
   
   Commands:
   
     cd Featflow2/benchmark
     ./configure
     
     -> The configure script list all directories
        "kernel_*", "apps_*" and "area51_*" and extract the application
        name from the "*".
     -> Let's assume, we have the following subdirectories:
          "apps_cc2d"
          "apps_pp2d"
          "kernel_triatest"
        Then this results in three applications:
          "cc2d"
          "pp2d"
          "triatest"
        Application names must be unique without the
        kernel_ / apps_ / area51_ -qualifier in front;
        therefore, "apps_cc2d" + "area51_cc2d" are not allowed!
     -> In every application subdirectory, a GNUmakefile is created
        for that particular application by calling the application
        specific "configure" script in that directory.
        
     make feat2Benchmark-[application-name]
     
     -> This command starts compiling the application 
        [application-name]. It generates an executable
           "feat2benchmark-[application-name]"
        in the benchmark directory.
     -> Example:
          make feat2benchmark-cc2d
        Compiles the application "cc2d" in the subdirectory
        apps_cc2d. The result is the executable
          "Featflow2/bechmark/feat2benchmark-cc2d"
        
      make [testcase-file] run
      
      -> Starts the tests listed in the file 
           "Featflow2/benchmark/[testcase-file].fbconf"
           
         More precisely:
         
      -> At first, all files
           "Featflow2/benchmark/tests/*.fbconf"
         are parsed to find all testcases.
         Each testcase has a unique id "test-id" and is associated to a
         specific application defined by the "appl =" parameter
         in the .fbconf file.
         
         Example: The test id "CC2D_001" in "tests/cc2d.fbconf"
           defines a test with the application "appl = cc2d".
           
      -> The file [testcase-file].fbconf contains a list of all
         test-id's that should be executed. The "run" parameter 
         instructs the make command to start all these test cases.
         
         Example: To start the testcase with id "CC2D_001", one can
           use the following two commands:
           
             echo CC2D_001 > cc2d_problem.fbconf
             make cc2d_problem run
           
           The test id "CC2D_001" defines that the underlying
           application is "cc2d" because of "appl = cc2d". Its 
           executable is therefore named "feat2benchmark-cc2d". The 
           application is started in the benchmark folder by the
           following command:
           
             "./feat2benchmark-cc2d ./data/apps_cc2d/master.dat"
             
           where "./data/apps_cc2d/master.dat" comes from the 
           "datfile  =" parameter associated to "CC2D_001".
           
       -> During the execution, the environment variable $LOGFILE
          points to a unique directory that must receive logfiles.
          A file named 
            "$LOGDIR/benchmarkresultfile"
          must be manually created by the application containing 
          deterministic values associated to this test case (on the 
          current architecture with the current build target). The file 
          must contain exactly those results that should be compared to
          reference results.
          
       -> The result file
            "$LOGDIR/benchmarkresultfile"
          is copied to the result directory and named
            "results/[build-target]/test[test-id].results".

          For example, the "CC2D_001" testcase on architecture
          "cc2d-opteron-linux" generates the file
            "$LOGDIR/benchmarkresultfile"
          which is copied to to
            "results/cc2d-opteron-linux/testCC2D_001.results".

       -> In a last step, the result files of all test cases in
            "results/[build-target]/*".
          are diff'ed to the reference results in
            "refsol/[build-target]/*"
          to find differences.
          
          For example, in the "CC2D_001" testcase the files
            "results/cc2d-opteron-linux/testCC2D_001.results".
          and
            "refsol/cc2d-opteron-linux/testCC2D_001.results".
          are compared.
          
       -> To update reference results, the user can copy the
          new result files from "results/..." to "refsol/...".
          

Automatically checking the repository
-------------------------------------
In the subdirectory ./bin there is a script called 
'runregressiontest_feat2'. This script is designed to be called in a
cron job every night. It checks out the current repository, compiles
the benchmark applications, compares the results and notifies the
users whether the tests run successfully. There are a couple of
options available in this script to configure its behaviour.

When the scrip is started by

  "runregressiontest_feat2 -t [testcase-file]",
  
it does basically the following:

a) It retrieves a clean checkout of the repository
b) It switches to the benchmark directory
     "Featflow2/benchmark"
c) It lists all possible build-id's on that machine
d) For every build-id, the corresponding commands are executed to start
   the benchmark:
   
     ./configure --id=[currentid]
       -> Set up the makefiles
       
     make benchmark  
       -> Builds all applications
       
     make [testcase-file] run
       -> Runs all testcases; these are defined as a list of id's
          in the file "[testcase-file].fbconf"
          
       -> Example: 
            "runregressiontest_feat2 -t alltests",
          leads to
            "make alltests run"
          searching for the tests in "alltests.fbconf".

     Instead of "alltests", one could also specify another set
     of tests like "dailytests". For every test set, a corresponding
     .fbconf file must exist (e.g. "dailytests.fbconf").
