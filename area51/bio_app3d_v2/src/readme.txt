Testing for the analytical solution (code validation)

1) Initial (=1) and boundary conditions are done.
15:28   28_09_2009

2) Seems like ready to go framework for the analytical solution
17:11	28_09_2009

3) Something with f_c is not right.
14:54   29_09_2009

4) Testing the 3D case
17:27	19.10.2009

5) simlpe chemo=x+y+z, u=0 works.
16:20	21.10.2009

6) Test c_t=\Delta c + f_c   
Somewhere is a mistake!!!!
17:23	21.10.2009 

7) Laplacian for c seems to work. But some problems in convergence.
10:11 	22.10.2009

8) Laplacians for chemo and cells work. Initial chemo=cell=x*(x-16)*...*(z-16)
16:21 	28.10.2009

9) Solve 
(M + D_2 \delta t L + \alpha \delta t M) c^{n+1} = M c^n + \beta \delta t M u^n
(M + D_1 \delta t L) u^{n+1} = M u^{n}
with initial chemo=cell=x*(x-16)*...*(z-16) 
(so everything but convective-chemotactical term)
14:51	02.11.2009

10) Solve (correct)
(M + D_2 \delta t L + \alpha \delta t M) c^{n+1} = M c^n 
(M + D_1 \delta t L) u^{n+1} = M u^{n}
with initial chemo=cell=x*(x-16)*...*(z-16) 
14:51	02.11.2009

11) Ready to go code, but without convective-chemoattractive term.
(M + D_2 \delta t L + \alpha \delta t M) c^{n+1} = M c^n + \beta \delta t M u^n
(M + D_1 \delta t L) u^{n+1} = M u^{n}
with initial chemo=cell=x*(x-16)*...*(z-16) 
12:10	03.11.2009

12) $grad c$ is not evaluated, but written as it is (1,1,1).
u=x*(16-x); c=x+y+z
16:34	10.11.2009

13) The same as 13, but everywhere are lumped matrices used.
16;07   12.11.2009

14) Cleaned code. Gradient c stil doesn't work.
15:23   13.11.2009

15) Before testing for Prof. Turek.
16:27   24.11.2009

16) Working code. Before remaking it into coupled version.
16:03   07.12.2009

17) Relaxation is done.
10:58   08.12.2009

18) Norm |analyt - numer| can be meassured
12:52   08.12.2009

19) Before working with Mu^n in chemotactic equation.
10:50   09.12.2009

20) good code to test difference between (x*(16-x),\phi) and (u,\phi), where u=x*(16-x).
11:49   09.12.2009

21) Just tedious validation of THIS code. No mistakes occured so far.
11:18   12.01.2010

22) Before switching to Robert's modified code.
11:02 	03.02.2010

23) Analytical-numerical difference with the help of L2 norm is realized.
12:49	04.02.2010

24) The H1 norm for ||c-c_h|| and ||u-u_h|| is realised.
13:53	04.02.2010

25) Before realising Q_2 elements.
10:24	11.02.2010 

26) Frist good results for the validation of the 3D bio code are obtained.
16:02	17.02.2010
