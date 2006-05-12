#!/usr/local/bin/gawk
#
# Extracts from the Featflow benchmark log file the timing 
# information in the form suitable for submission 
# to the Featflow benchmark webpage.
#
# To call this script manually, first start the Featflow benchmark,
# then switch to the Featflow installation directory and type the
# following command:
#
#    awk -f bin/extract.awk bench-xxxxxxxxxxxxxxx.log
#
# with "bench-xxxxxxxxxxxxxxx.log" being the name of the log file
# that was created by the benchmark. The script will print the
# timing results on screen. Make sure that your console shows
# 132 columns, otherwise you can't see all information on the screen :-)

BEGIN {i1=0; i2="x"; dim=3;}

/Machine/ {host=$2; id=$4}
/convective part/ {if($4==1) type="upw"; else type="sd";}
/ILEV,NVT,NMT,/ {dim=2;}
/ILEV,NVT,NAT,/ {dim=3;}
/total time :/ {cpu=$4;}
/grid  time :/ {grid=$4;}
/post  time :/ {post=$4;}
/lin.  time :/ {lin=$4;}
/mavec time :/ {mavec=$5;}
/konv. time :/ {konv=$5;}
/bdry  time :/ {bdry=$5;}
/LC    time :/ {lc=$5;}
/ILU   time :/ {ilu=$5;}
/mg    time :/ {mg=$4; i2="c";}
/U-mg  time :/ {umg=$4; i2="p";}
/P-mg  time :/ {pmg=$4; i2="p";}

/MULTIGRID COMPONENTS/ {i1=1;}

{
  if(i1==1) {
    if(i2=="c") {
      printf("cc%1dd-%s \t%d\t%d\t%d\t%d\t%d\n",
             dim,type,int(cpu+0.5),int(grid+post+0.5),
             int(mavec+lc+0.5),int(konv+0.5),int(mg+0.5));
    }
    if(i2=="p") {
      printf("pp%1dd-%s \t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
             dim,type,int(cpu+0.5),int(grid+post+0.5),
             int(mavec+lc+0.5),int(ilu+0.5),int(konv+0.5),
             int(umg+0.5),int(pmg+0.5));
    }
    i1=0; i2="x";
  }
}

END   {
  printf("#host: %s   id: %s\n",host,id);
}

