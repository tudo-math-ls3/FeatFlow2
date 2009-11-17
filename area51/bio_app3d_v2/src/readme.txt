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

