The SYSUTILS library here is used for capsuling machine-dependent
routines. In particular this capsules the timing-routine ETIME into
another routine ZTIME which is used in FeatFlow. ETIME is standard
on most machine architectures, but there are some architectures not
supporting this command. In this case the routine ZTIME has to be
modified to the appropriate machine dependent GetTime-routine.