read2dprm prm bench1.prm
read2dtri tri bench1.tri --boundary prm
meshhierarchy mhier --mesh tri --boundary prm --levels 5
fehierarchy fehier --meshhierarchy mhier --boundary prm --components 3 --quadelements EL_EM30 EL_EM30 EL_Q0 --quadcub CUB_G4x4 CUB_G4x4 CUB_G4x4
mlevelprjhierarchy mlprj --fehierarchy fehier

createblockvector vec1 --fehierarchy fehier --level 5
createblockvector vec2 --fehierarchy fehier --level 4

interpolatevector vec1 5 vec2 4 --fehierarchy fehier --mlprjhierarchy mlprj
