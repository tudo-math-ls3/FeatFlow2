Version 1.0 - 16.6.2014

This is an Extended Finite Elemente Comparer and some sort of postprocessing tool.
When you did your calculations and saved the solution vector of the resulting equation
you can do some postprocessing with this tool.
This tool covers the most gernal 2D-case that is possible:
You have a function f1 and a function f2 which were calculated using different finite elements
and even different geometries (which is interesting for penalty methods,
testing which influence a different grid has,...)


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

Important: This tool really needs 2 input vectors, if you put in only 1 then it will crash.
However I left the option for a "hack": If you want to analyze only 1 function, you can
enter it 2 times as input file but deactivate the calculation of the things you do not
want to calculate.


If you need more features or have any questions, do not hesitate to contact me.

Malte Schuh
malte.schuh@math.tu-dortmund.de