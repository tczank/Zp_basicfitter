#!/bin/bash

if [  -f ./whizard0eff_all.dat ]; then
    rm ./whizard0eff_all.dat
fi

parameter_map=(712 1212 1712 2212 2712 3212 3712 4212 4712 5212 5712 6212 6712 7212 7712 8212 8712 9212 )

for map in "${parameter_map[@]}";do
    root -q -b -l -s 'fitreborn_v4.C("../whizard0/darkz_whizard0_'${map}'_all.root")' > whizard0eff.${map}.dat
    sed -i '1,2d' whizard0eff.${map}.dat
done

for map in "${parameter_map[@]}";do
cat whizard0eff.${map}.dat >> whizard0eff_all.dat
rm whizard0eff.${map}.dat
done

