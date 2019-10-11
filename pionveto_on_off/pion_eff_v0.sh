#!/bin/bash

if [  -f ./pion_eff_all_newbg.dat ]; then
    rm ./pion_eff_all_newbg.dat
fi

parameter_map=(212 712 1212 1712 2212 2712 3212 3712 4212 4712 5212 5712 6212 6712 7212 7712 8212 8712 9212 )

for map in "${parameter_map[@]}";do
    root -q -b -l -s 'fitreborn_pivon_v0.C("../../weighted_pion_veto/darkz_isr_'${map}'_all_wpvon.root")' > pion_eff.${map}_newbg.dat
    sed -i '1,2d' pion_eff.${map}_newbg.dat
done

for map in "${parameter_map[@]}";do
cat pion_eff.${map}_newbg.dat >> pion_eff_all_newbg.dat
rm pion_eff.${map}_newbg.dat
done

gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=../../weighted_pion_veto_on/fitsum.pdf ../../weighted_pion_veto_on/*root.eps
