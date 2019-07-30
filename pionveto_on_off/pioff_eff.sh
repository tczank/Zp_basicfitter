#!/bin/bash

if [  -f ./pioff_eff_all.dat ]; then
    rm ./pioff_eff_all.dat
fi

parameter_map=(212 712 1212 1712 2212 2712 3212 3712 4212 4712 5212 5712 6212 6712 7212 7712 8212 8712 9212 9712 10000 )

for map in "${parameter_map[@]}";do
    root -q -b -l -s 'fitreborn_pivoff.C("../../pion_veto_off/darkz_isr_'${map}'_all_pvoff.root")' > pioff_eff.${map}.dat
    sed -i '1,2d' pioff_eff.${map}.dat
done

for map in "${parameter_map[@]}";do
cat pioff_eff.${map}.dat >> pioff_eff_all.dat
rm pioff_eff.${map}.dat
done

