#!/bin/bash

if [  -f ./bias_mc_all.dat ]; then
    rm ./bias_mc_all.dat
fi

parameter_map=(212 712 1212 1712 2212 2712 3212 3712 4212 4712 5212 5712 6212 6712 7212 7712 8212 8712 9212 )

for map in "${parameter_map[@]}";do
    root -q -b -l -s 'bias_study_mc.C("../../weighted_pion_veto/darkz_isr_'${map}'_all_wpvon.root")' > bias_mc.${map}.dat
    sed -i '1,2d' bias_mc.${map}.dat
done

for map in "${parameter_map[@]}";do
cat bias_mc.${map}.dat >> bias_mc_all.dat
rm bias_mc.${map}.dat
done

