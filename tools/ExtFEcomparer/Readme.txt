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