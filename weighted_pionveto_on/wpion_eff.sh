#!/bin/bash

if [  -f ./wpion_eff_all.dat ]; then
    rm ./wpion_eff_all.dat
fi

parameter_map=(212 712 1212 1712 2212 2712 3212 3712 4212 4712 5212 5712 6212 6712 7212 7712 8212 8712 9212 )

for map in "${parameter_map[@]}";do
    root -q -b -l -s 'fitreborn_pivon.C("../../weighted_pion_veto/darkz_isr_'${map}'_all_wpvon.root")' > pion_eff.${map}.dat
    sed -i '1,2d' pion_eff.${map}.dat
done

for map in "${parameter_map[@]}";do
cat wpion_eff.${map}.dat >> wpion_eff_all.dat
rm wpion_eff.${map}.dat
done

