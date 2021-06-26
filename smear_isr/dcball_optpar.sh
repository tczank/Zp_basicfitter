#!/bin/bash

if [  -f ./opt_dcball_par_all.dat ]; then
    rm ./opt_dcball_par_all.dat
fi

parameter_map=(215 315 415 515 615 715 815 915 1015 1215 1415 1615 1815 2015 2215 2415 2615 2815 3015 3215 3415 3615 3815 4015 4215 4415 4615 4815 5015 5215 5415 5615 5815 6015 6215 6415 6615 6815 7015 7215 7415 7615 7815 8015 8215 8415 8615 8815 9015 9215 9415 9615) # 9815 10015)

for map in "${parameter_map[@]}";do
# at home
    root -q -b -l -s 'dcball_paropt.C("~tczank/Documentos/New_isr_smeared/Newisr'${map}'_smeared.root")' > opt_dcball_par.${map}.dat
# not home
#    root -q -b -l -s 'dcball_paropt.C("~tczank/Documentos/smearisr/Newisr'${map}'_smeared.root")' > opt_dcball_par.${map}.dat
    sed -i '1,2d' opt_dcball_par.${map}.dat
done

for map in "${parameter_map[@]}";do
cat opt_dcball_par.${map}.dat >> opt_dcball_par_all.dat
rm opt_dcball_par.${map}.dat
done

