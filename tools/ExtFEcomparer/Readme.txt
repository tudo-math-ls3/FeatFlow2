Version 2.1.2   19.11.2014

Added: Now it is also possible to specify an element-list. This means:
i.e. if your element-pair is not provided you specify a list with element identifiers.
Identifier i is used for component i.

Added: it is possible to overwrite the parameters with the command line
and specifiying a master-dat with the command-line.
For a master.dat type ./ExtFEcomparer --smaster=./path/to/my/master.dat
To overwrite a parameter use --D<section>.variable:<entry>=value

Added: A Benchmark-Modus which does not need any read-in-data - good for
testing new features/checking they still give right results.

Modified: Use another projection for the UCD-Output which is more general and also
supports 1D and (for future purpose) 3D

Fixed: Output of L2-Norms to a file had a wrong comment (but did the right work)

Version 2.1.1   7.11.2014

To sum up the features and present the changelog, we need a notation:
f_k = Component k of function f, for example if we solve 2D-Navier-Stokes,
f_1 = Velocity in x direction, f_2 = Velocity in y direction, f_3 = pressure
f is created by a representing vector of coefficents (i.e. the solution vector
of the Navier-Stokes-Equations from cc2d), the domain (represented by 
*.tri and *.prm files) and the finite elements that were used, represented
by a name or (in case of Navier-Stokes) an identifier for the element pair.
All this has to be present to actually use the tool (for 2 functions).
Hint: If you want to analyze only 1 function, enter it 2 times as input:
One time as f, one time as g

Summary of features:

1) Calculations for 1D and 2D 
1.1) Calculate ||f_k||_L2(\Omega)
1.2) Calculate ||g_k||_L2(\Omega)
1.3) Calculate ||f_k - g_l||_L2(\Omega)
1.4) Calculate (d^i/dx^i) f_k(x) |i|=0,1 i: multiindex
1.5) Calculate (d^j/dx^j) g_k(x) |j|=0,1 j: multiindex
1.6) Calculate (d^i/dx^i)f_k(x) - (d^j/dx^j)g_l(x), |i|=0,1; |j| = 0,1 i,j: multiindex

2) Output
2.1) Export the vector of coefficents in an other format
2.2) UCD-Output of the mesh
2.3) UCD-Output of the input functions

Changelog:
-Bugfixes in the UCD-Output + fine-tuning the UCD-Output
-Added a warning if the vector does not match the discretisation (only for the format
 int N, char sVarname, ...)
Fixed typos in the logfile_settings.dat - the logfiles were deactivated because of the typo


On request I will add more features, if someone is willing to test I will add
1) More formats for UCD-Output
2) Discretisation with 1 element in 2D (results from poisson app for example)
3) H1-Norm-Calculations
4) Calculating the L2-Norm of the divergence of the input functions

Malte Schuh

Version 2.1   4.11.2014

To sum up the features and present the changelog, we need a notation:
f_k = Component k of function f, for example if we solve 2D-Navier-Stokes,
f_1 = Velocity in x direction, f_2 = Velocity in y direction, f_3 = pressure
f is created by a representing vector of coefficents (i.e. the solution vector
of the Navier-Stokes-Equations from cc2d), the domain (represented by 
*.tri and *.prm files) and the finite elements that were used, represented
by a name or (in case of Navier-Stokes) an identifier for the element pair.
All this has to be present to actually use the tool (for 2 functions).
Hint: If you want to analyze only 1 function, enter it 2 times as input:
One time as f, one time as g

Summary of features:

1) Calculations for 1D and 2D 
1.1) Calculate ||f_k||_L2(\Omega)
1.2) Calculate ||g_k||_L2(\Omega)
1.3) Calculate ||f_k - g_l||_L2(\Omega)
1.4) Calculate (d^i/dx^i) f_k(x) |i|=0,1 i: multiindex
1.5) Calculate (d^j/dx^j) g_k(x) |j|=0,1 j: multiindex
1.6) Calculate (d^i/dx^i)f_k(x) - (d^j/dx^j)g_l(x), |i|=0,1; |j| = 0,1 i,j: multiindex

2) Output
2.1) Export the vector of coefficents in an other format
2.2) UCD-Output of the mesh
2.3) UCD-Output of the input functions

Changelog:
- Permanently activated the posibility to write out the vector again.
  This might seem weird, but has 2 reasons:
  1) If you get numbers that you don't trust, you can check if
     the vector is read in correct
  2) If you have a vector that is unformatted/binary, this format is - from what
     I learned - processor dependent. If you run your code on computer A and want/must
     use the ExtFEcomparer on another computer (ie if you get results from computer A and B)
     then you can run the ExtFEcomparer on both computers to convert the vectors
     to formatted form, afterwards you can run the ExtFEcomparer on one computer
     and compare the two results.

- Added the function to write out the mesh.
  This is mostly for debugging, if you get numbers you don't trust you can write
  out the mesh that was used and check if it is reconstructed in a correct way
  or if it results in trouble when refining it
  Note that this will only work if you enter 2 vectors with the right size as well

- Added the function to do a UCD-Output (paraview) of the input functions
  This has two reasons:
  1) If you get numbers you don't understand, you can have a look on your
     input solutions
  2) Maybe you just saved the vector and later decided you need a graphical output
     This output cannot be as good as a perfectly set up output for you problem
     because the code does not know i.e. that in component 1 we have velocity_x and
     in component 2 we have velocity_y, so there is no chance to bundle that in one
     output vector - but if you are good in manipulating files, you can later manipulate
     this to one vector.
     Another point is that the code does not know if component 3 is pressure or not, so
     it will print out this one with the vertex_based variant and not by an element_based
     like it is usually done. So it does look different - but it is better than nothing.

Changes under the hood:
Restructured the code in a way that every analysis tool is some sort of it's own module
and clearly seperated from the others. This means in particular: We have one routine
"init_postprocessing" that is called by the main routine. In this routine we have calls to
"init_postprocessing_ucd",... which do the real work. Idea behind this: keep the code readable
and make it easier to add new features/find errors.

General structure of the code:
There is one routine called ExtFEcomparer_calculate in the core which* is the one that
"does something" with the FE-functions and produces output in the terminal (only there!)
All information what shall be done and the results are stored in the postprocessing structure
already which is the only thing that is seen by the postprocessing routines.
There the file-output is generated.
Why is this way chosen? Simply because the calculation-routines are not easy, we do not
want to have too much code there so that it is easier to debug in any sense.
Also we want to minimize the write-commands in the calculation routines as it makes it
easier to parallize them if there every data is already present/can be catched easily.
Another benifit is that this - in my opinion - is the best way to easily work on the output,
i.e. if you just want to work on the output, you look in the output-routine and can assume
that any information you need is in there, you don't have to search at which place in the
calculation routines the variables you need have the right values.


Future plans: tuning the ucd-output + adding ucd-output of f-g on a region of interest
+ adding the option to calculate flux-values

*to be more precise: There are again just routines which call (in the same file) the real
working routines.

Malte Schuh

Version 2.0 - 21.10.2014

This is complete rewrite of the Version 1.0. Many things changed only under the hood, 
some change for the user

Changes for the user:
1) Added new vector types that can be read in:  Unformatted (i.e. non-human-readable), with the file-format:  
int N, char svarname, followed by N values - and this one for every variable
Previously only every vector that was written out by the subroutine
write_BlockVectorHR(...) could be read in
2) Added Support for discretisations with only 1 element type
3) Added 1D
4) Added more calculation options. The cost of this is that
 now you have to know which variable in which vector is what - which should not be too hard,
 i.e. for cc2d variable 1 is the velocity in x direction, variable 2 the velocity in y direction and 
  variable 3 the pressure.
5) Added the possibility to set up every calculation by itself: Own mask function and own cubature rule

Changes under the hood:
- Added functions that allow a 1D,2D and 3D support to be realised very easy
- deleted all not-neccesary data-types
- restructured the read-in of the vectors
- restructured the creation of the discretisation 
-restructured the creation of the parametrisation

Still necessary: If you want to analyze only 1 function, you still
have to enter 2 vectors - but you can set the calculations in a way
that the second one is totally ignored. 

Hint: If you get results you do not trust, switch of OpenMP



Version 1.0 - 16.6.2014

This is an Extended Finite Elemente Comparer and some sort of postprocessing tool.
When you did your calculations and saved the solution vector of the resulting equation
you can do some postprocessing with this tool.
This tool covers the most gernal 2D-case that is possible:
You have a function f1 and a function f2 which were calculated using different finite elements
and even different geometries (which is interesting for penalty methods,
testing which influence a different grid has,...).
In fact, all you need so you can work with this tool is some geometrical
common area.


It is coded by Malte Schuh (malte.schuh@math.tu-dortmund.de)

To describe what it can do, we need some definitions first
f* always refers to a FEM-function
u always refers to the entire velocity
p always refers to the pressure.

This tool can calculate the following norms
- ||u1 - u2||_L2(\Omega1)
- ||p1 - p2||_L2(\Omega1)
- ||u1||_L2(\Omega2)
- ||p1||_L2(\Omega2)
- ||u2||_L2(\Omega3)
- ||p2||_L2(\Omega3)

However, you need to be able to specify \Omega1, \Omega2 and \Omega3. This is usually done with
an IF-statement in the accoding configuration file, i.e.:
If (x < 0.1,1,0) specify \Omega as the subset of the mesh which have an x-component lower than 0.1
To find out all supported expressions read the parser documentation.
You also can to specify the cubature rule that is used for f1 and f2, i.e. G4x4 or AUTO_G2 - see the Featflow-documentation
to find out which cubature rules are suppored.

The tool can also calculate point values and the difference, i.e.
f1(p)
f2(p)
|f1(p) - f2(p)|.

At the moment point values and values of the first derivatives are supported.

If you want it, the results are written out into a text-file.


At the moment this tool supports only 2D and only vectors which are written out using:
- formatted output (human readable, not processor dependend)
- no sorting strategy (which is standard in i.e. cc2d)

Important notes: 
1) This tool really needs 2 input vectors, if you put in only 1 then it will crash.
However I left the option for a "hack": If you want to analyze only 1 function, you can
enter it 2 times as input file but deactivate the calculation of the things you do not
want to calculate.
2) This tool does not check if you input is correct, to be precise it will not check if you
entered a solution on level 7 but set NLMAX for this function only to 3 or if you calculated
the solution i.e. with Q2/P1 but entered that it was calculated with Q1~(EM30) / Q1~(EM30) / Q0.
When you do this and analyze this function/compare it to another function, this code has
unspecified behavior - this totally depends on the type of input-error.
However, this is not a bug - it is the result of note 1): If you want to analyze 1 function (and
not compare it to another function) then you can feed the program with any input for
the second function and deactivate that you calculate anything with the second function.


If you need more features or have any questions, do not hesitate to contact me.

Malte Schuh
malte.schuh@math.tu-dortmund.de
