#!/bin/bash

if [  -f ./dcball_par_all.dat ]; then
    rm ./dcball_par_all.dat
fi

parameter_map=(215 315 415 515 615 715 815 915 1015 1215 1415 1615 1815 2015 2215 2415 2615 2815 3015 3215 3415 3615 3815 4015 4215 4415 4615 4815 5015 5215 5415 5615 5815 6015 6215 6415 6615 6815 7015 7215 7415 7615 7815 8015 8215 8415 8615 8815 9015 9215 9415 9615) # 9815 10015)

for map in "${parameter_map[@]}";do
    root -q -b -l -s 'dcball_newfit.C("~tczank/Documentos/smearisr/Newisr'${map}'_smeared.root")' > dcball_par.${map}.dat
    sed -i '1,2d' dcball_par.${map}.dat
done

for map in "${parameter_map[@]}";do
cat dcball_par.${map}.dat >> dcball_par_all.dat
rm dcball_par.${map}.dat
done

