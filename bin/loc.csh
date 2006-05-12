#!/bin/tcsh

cd ..
set lol = 0
set loc = 0

foreach lauf ( `find . -name "*.f"` )

  set tmp = `wc $lauf`
  set l = `echo $tmp | awk '{print $1}'`
  set c = `echo $tmp | awk '{print $3}'`

@ lol = $lol + $l
@ loc = $loc + $c

end

echo Zeilen: $lol  Zeichen: $loc
