The basic regression test suite -- short introduction for dummys
----------------------------------------------------------------
The FEAT regression test suite allows to run a couple of benchmarks,
interactively and non-interactively via a cronjob at night to validate
the code, i.e. to ensure that everything compiles and runs
correctly. The basic concept is that a couple of test applications run
on a couple of machines producing some output which is compared
against reference solutions. E-mails are generated in case
applications fail to compile or tests fail to reproduce the reference
solution. Given that the number of commits per day is reasonably small
it is easy to pinpoint the user whose improper checkin is to blame.

General structure of the benchmarks - Quick Start Guide
-----------------------------------------------------------
a) Every benchmark application has its own subdirectory, e.g.
     "Featflow2/benchmark/apps_cc2d"
   that may contain source files, but needs at least a configure
   script in order to be able to generate a Makefile.
    
b) The benchmark applications are, in contrast to the standard
   behaviour, all stored directly in the main benchmark directory,
   i.e.  
     "Featflow2/benchmark".
     
c) The subdirectory
     "Featflow2/benchmark/data"
   contains a couple of subdirectories for every application containing
   the data files of the corresponding test(s). The data files for the
   standard cc2d test can be found e.g. in the directory
     "Featflow2/benchmark/data/apps_cc2d"
   This slightly differs from the standard directory structure of a
   FEAT application that stores its data in a "data/" subdirectory.

d) For every benchmark application, the directory
     "Featflow2/benchmark/tests"
   contains one or more files that contains test definitions for a
   given application. The test definitions for the standard cc2d
   application can be found, e.g., in the file
      "Featflow2/benchmark/tests/cc2d.fbdef"
   Each test starts with a parameter "testid=..." that defines the
   test ID. The parameter "appl =..." defines the name of the
   application. The name of the .fbconf file should coincide with it
   or at least contain it.

e) When an application is started, the following data is passed to it:

   - The name/path of the data file specified by the "datfile=..." 
     parameter in the sections of the corresponding .fbdef file
     is passed as 1st command line parameter.

     For the cc2d application with test ID CC2D_001, e.g., the "datfile"
     parameter is defined as
       "datfile  = ./data/apps_cc2d/master.dat"
     which points to the master file that controls the application
     relative to the benchmark directory.
     The application is then executed by calling
       "feat2benchmark-cc2d ./data/apps_cc2d/master.dat"
       
   - The parameter "mglevels = x,x,x,..." defines a list of one or
     several refinement levels. For every value in the list, a
     separate benchmark is started. The current level is then passed
     as environment variable
       $MGLEVEL = ...
     to the application.
     
   - Below each test ID, the user can specify an arbitrary number
     of variables of arbitrary name which are passed as environment
     variables to the application.
     
     Example for cc2d.conf:
       testid   = CC2D_001
       datfile  = ./data/apps_cc2d/master.dat
       mglevels = 4,5
       element  = 11
       eps      = 1E-10
       
     Defines test case "CC2D_001". The three variables 'testid',
     'datfile', 'mglevels' have a special meaning that have already
     been explained. The remaining ones are defined as environment
     variables in the benchmark control script 'runtests'
       $ELEMENT = 11
       $EPS     = 1E-10
     This script, when run, calls the application two times with the
     environment variable $MGLEVEL set to $MGLEVEL=4 and $MGLEVEL=5.
     
   - For convenience, parameters defined in one test id are inherited
     by the next one until redefinition. This to avoid having to
     retype the same settings over and over again.
     
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
     $LOGDIR = directory to application is supposed to write its log
               files to
     $RESULTFILE = Name of the result file for benchmark data.
   $LOGDIR is set to 'logs/' + 'testid'. The directory is automatically
   created along with the 'runtests' script. Having a unique name
   (because obviously all test IDs are unique), several tests can be
   run simultaneously, e.g. in a queueing environment, with
   interference, i.e. it won't happen that two instances of the same
   application try to write to the same files. 
   Every benchmark application is supposed to create a file named
     $LOGDIR/$RESULTFILE = "$LOGDIR/benchmarkresultfile"
   that contains deterministic benchmark results that allow comparison
   of this program run to previous runs by means of reference
   results. Obviously, timing results and similar nondeterministic
   results should not be written to this file as the file is later on
   diff'ed to a file with reference results.
   
g) SUMMARY: The data flow for executing a benchmark application is
   as follows:
   
   Commands:
   
     cd Featflow2/benchmark
     ./configure
     
     -> The configure script walks through all directories matching
        "kernel_*", "apps_*" and "area51_*" and creates a Makefile in
        every one of them (invoking the "configure" script in that
        subdirectory) in order to be able to build the application.
     -> Benchmark application names are built according to the
        following rules:
        - They have a common prefix 'feat2benchmark-'
        - The suffix is the part matching '*' of the respective
          subdirectory names:
          Example: Let's assume, we have the following subdirectories:
             "apps_cc2d"
             "apps_pp2d"
             "kernel_triatest"
          Then this results in three applications:
             "cc2d"
             "pp2d"
             "triatest"
          which results in the following binaries being created:
             "feat2benchmark-cc2d"
             "feat2benchmark-pp2d"
             "feat2benchmark-triatest"

        - Obviously, it is mandatory that application names are
          unique, i.e.  having two subdirectories named "apps_cc2d"
          and "area51_cc2d" is not allowed as you'd get a conflict in
          the binary names.

        
     make feat2benchmark-[application-name]
     
     -> This command starts compiling the application 
        [application-name]. It generates an executable
           "feat2benchmark-[application-name]"
        in the benchmark directory.
     -> Example:
          make feat2benchmark-cc2d
        Compiles the application "cc2d" in the subdirectory
        apps_cc2d. The result is the executable
          "Featflow2/bechmark/feat2benchmark-cc2d"

 
      make benchmark

      -> Is a generic command that compiles all benchmark applications
         "feat2benchmark-*".        


      make [testcase-file]         or        make [testID]
      
      -> Create benchmark run control script 'runtests' for all tests
         listed in the file
           "Featflow2/benchmark/[testcase-file].fbconf"
           
         More precisely:
         
      -> At first, all files matching
           "Featflow2/benchmark/tests/*.fbconf"
         are parsed to find all testcases.
         Each testcase has a unique test ID named "testid" and is
         associated to a specific application defined by the "appl ="
         parameter in the .fbconf file, and to an application class
         defined by the "class =" parameter in the .fbconf file.
         
         Example: The test ID "CC2D_001" in "tests/cc2d.fbconf"
                  defines a test to be run with help of the
                  application "appl = cc2d". The associated test class
                  is "DEFAULT".
         
         The subroutine "fb_runAndEvalBenchmarkTest()" in the file
         "runtests.template" defines the different test classes.
         The test class "DEFAULT" is the most popular.
         
         Other test classes may be added here to support an extended
         verification of solutions, e.g. the calculation of error
         reduction rates or similar. The DEFAULT implementation
         just executes the application for a couple of levels
         and a specific configuration file.
           
      -> The file [testcase-file].fbconf contains a list of all
         test IDs that should be executed. Comment lines start
         with the hash character.
         Alternatively, it is possible to directly specify the
         test ID on the command line.


      ./runtests                   or        make run

      -> The make target "run" is provided for convenience only and
         does nothing more than started the benchmark execution by
         invoking the command './runtests'.


      In queueing environment like LiDO-1, LiDO-2, NEC SX-8/9 or on the 
      D-GRID cluster (DGRZR) one most likely does *not* want to run the
      simulations on the gateway node one is logged into. Instead, a
      jobfile needs to be created per test ID and submitted to the queue.
      Currently, Feat2 only has one of these automatic submission scripts:
      
      bin/lido_schedule_tests [testcase-file]         or
      bin/lido_schedule_tests [testID]

      -> The script does basically the same as "make [testcase-file|testID]",
         but beyond that also automagically determines the most
         appropriate queue based on 
         * configured settings for walltime (hardcoded in the header
           of bin/lido_schedule_tests or specified by the parameter "-t"),
         * the necessary number of processes (which is trivially 1 given
           that Feat2 supports serial applications only) 
         * the interconnect (trivial for Feat2 as well: ethernet)
         Then, the script tries very hard to submit the job, in
         particular in case the queue is so full that no more jobs are
         accepted at first.
      -> One can specify multiple testcase-ID's or multiple testcase files
         in the command line of lido_schedule_tests, e.g.
           bin/lido_schedule_tests CC2D_001 CC2D_002 ...
         or, resp.,
           bin/lido_schedule_tests tests1.fbconf tests2.fbconf ...
         with test1.fbconf, test2.fbconf, ... files containing test cases.
         All tests are independently scheduled and executed in parallel!
      -> The function fb_setMailSender() in the file 
         include/lib_for_xxx_schedule_tests maps the username to a valid
         email address. If lido_schedule_tests stops with an error that
         the user is unknown, one has either to add a mapping here or
         tell lido_schedule_tests a valid email address using the
         "-m email address" parameter on the command line.


   Summary to set up and run benchmark in non-queueing environments:
     
     cd Featflow2/benchmark
     ./configure
     make benchmark alltests run


   Summary to set up and run on LiDO-1:

     cd Featflow2/benchmark
     ./configure
     make benchmark
     bin/lido_schedule_tests alltests.fbconf


   How is the benchmark run?
          
           The test ID "CC2D_001" defines that the underlying
           application is "cc2d" because of "appl = cc2d". So, the
           binary "feat2benchmark-cc2d" is used. The application is
           started in the benchmark folder by the following command:
           
             "./feat2benchmark-cc2d ./data/apps_cc2d/master.dat"
             
           where "./data/apps_cc2d/master.dat" comes from the 
           "datfile  =" parameter associated to "CC2D_001".
           
       -> During the execution, the environment variable $LOGDIR
          points to a unique directory that must receive logfiles.
          The environment variable $RESULTFILE defines the name of the
          benchmark log file where the application must write
          log data to, which is later compared to reference results.
          So, the application has to write all crucial, deterministic 
          values associated to this test case (on the current 
          architecture with the current build ID) to a file named 
            "$LOGDIR/$RESULTFILE" ( = "$LOGDIR/benchmarkresultfile" )
          
       -> The content of the result file
            "$LOGDIR/$RESULTFILE"
          is appended to the file
            "$LOGDIR/$RESULTFILE.all-level". 
          After looping over all refinement levels given in $MGLEVELS
          the resulting file is stored to the result directory and
          named:
            "results/[buildID]/test[testID].results".

          For example, the "CC2D_001" test ID on architecture
          "cc2d-opteron-linux" run exactly one time on level 4 and
          generates the file
            "$LOGDIR/$RESULTFILE"
          which is copied to
            "$LOGDIR/$RESULTFILE.all-level"
          and subsequently stored under
            "results/cc2d-opteron-linux/testCC2D_001.results".

       -> In a last step, the result files of all test cases in
            "results/[buildID]/*".
          are diff'ed to the reference results in
            "refsol/[buildID]/*"
          to find differences.
          
          For example, in the "CC2D_001" testcase the files
            "results/cc2d-opteron-linux/testCC2D_001.results".
          and
            "refsol/cc2d-opteron-linux/testCC2D_001.results".
          are compared.
          
       -> To update reference results, the user can copy the
          new result files from "results/..." to "refsol/...".
          

Automatic regression test of the repository
-------------------------------------------
In the subdirectory ./bin there is a script called
'runregressiontest_feat2'. This script is designed to be called in a
cron job every night. It checks out the current repository, detect
whether (relevant) changes have been made since the last time a
working copy has been made, compiles the benchmark applications,
compares the results and notifies the users whether the tests run
successfully. There are a couple of options available in this script
to configure its behaviour.

When the script is started by

  "runregressiontest_feat2 --checkout --compare-feat2-checkouts -t [testcase-file]",
  
it does basically the following:

a) It retrieves a clean checkout of the repository
b) It switches to the benchmark directory
     "Featflow2/benchmark"
c) It determines the default build ID for the current host and runs
   the benchmark for it:
   
     ./configure
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

     On queueing systems like LiDO a separate script is invoked.
     'bin/lido_schedule_tests' generates a jobfile for every test ID
     and submits it to the appropriate queue.

     The script then waits till all jobs have finished (which is
     trivial in case of sequential execution, a bit less in case a
     queueing environment is being used), processes the output, adds
     hyperlinks between deviated tests and sends this result via
     e-mail to the configured recipients.

d) The benchmark can be run at night for non-default build IDs by
   passing a distinct build ID. Example

    "runregressiontest_feat2 --checkout --compare-feat2-checkouts -t [testcase-file] --id=pc64-opteron-linux-gcc-goto"


Preparing an application for the use with the regression system
---------------------------------------------------------------
An application has to be prepared to be able to run with the regression
test system. The rule is as follows:

 * All data which should be automatically checked shall be written
   into the file 

         "$LOGDIR/$RESULTFILE" ( = "$LOGDIR/benchmarkresultfile" )

 * Timings and date/time data shall not be written to the result-file.
 
As indicated above, if called for multiple levels, the benchmark system
automatically concatenates the results of all levels into a file

         "$LOGDIR/$RESULTFILE.all-level"

and compare this file to the reference.

Setting up a Featflow2 application for the use with the benchmark system
involves a modification in the output rules of the genoutput.f90 module:


* Possibility 1: Log the complete application output

 This should be used with care since often, timings and/or date/time information
 also belongs to the output. However, it allows to log also kernel output.
 Two steps have to be done in the application:
 
 1.) Switch on the benchmark output file with a proper call to "output_init".
 2.) Redirect all output to the benchmark log file by changing cbenchLogPolicy.

 In practice, this looks as follows. Assume that in a data file, there 
 is a variable
   
   ...
   sbenchlog = "$LOGDIR/$RESULTFILE"
   ...
   
 The routine "output_init" has to be called with this log file specified
 as benchmark log file, e.g.
   
   call parlst_getvalue_string_direct (...,"sbenchlog",sbenchlog,bdequote=.true.)
   call output_init ("","",sbenchlog)
   cbenchLogPolicy = 2
   

* Possibility 2: Selective output to the benchmark log file

 This is similar to the above way, but changes the output policy for every
 "output_line" command. In practice, this looks as follows .Assume again that 
 in a data file, there is a variable
   
   ...
   sbenchlog = "$LOGDIR/$RESULTFILE"
   ...
   
 The routine "output_init" has to be called with this log file specified
 as benchmark log file, e.g.
   
   call parlst_getvalue_string_direct (...,"sbenchlog",sbenchlog,bdequote=.true.)
   call output_init ("","",sbenchlog)

 In the call to "output_line", the benchmark log policy has to be specified
 in the output-mode parameter.
 
   call output_line("A message only to the benchmark log.", &
       OU_CLASS_MSG, OU_MODE_STD + OU_MODE_BENCHLOG)

 Only those messages with OU_MODE_BENCHLOG set are written to the benchmark
 log file. Therefore, the application programmer has just to add this
 constant to the output mode in all output commands which should be written
 to the benchmark log.
 
 Remark: Spefifying OU_MODE_BENCHLOG without a benchmark log file set by
 "output_init" does not harm. There is just no file written to the hard disc.
