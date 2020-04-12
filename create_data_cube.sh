#!/bin/bash
shopt -s expand_aliases #allow the alias to work
# setup idl environment
# /export/data/apps/idl85/envi53/bin/envi_setup.bash

fileseq="2xSc_07_EC_BH_vw1e4_ReposNoRFBNoRP"
snapID=151
# create RSS data for all 3 dithering set
for ds in $(seq 0 2)
do
    echo "creating RSS for dithering set: $ds"
    ./create_RSS.sh $fileseq $snapID fib_radius=0.79 spec_style="star_age" rtfaceon=0 tauv=1.0 one_comp_dust=0 with_metal=1 cir_fib=1 with_PSF=1 real_fibpos=1 dither_set=$ds
done



# ------------------ This bit does not work very well --------------
# ------------------ manually run the python script in the last line if necessary
# wait the last batches of cell spectra calculation
# If there is no any idl calculation is ongoing, jump out the loop
# An calculation should not last for more than 5 mins.
# If we detect idl process for more than 5mins, it's likely that there are other ongoing calculation.
# Just ignore it and continue
for i in $(seq 0 60)
do
    idl_count="$(ps aux|grep idl|wc -l)"
    echo $idl_count
    # two idl lmgrd process and grep idl will persist when find idl_count
    # so idl_count >=3 at all times.
    if [ $idl_count -ge 4 ]
    then
        echo "IDL ongoing"
        sleep 5
    else
        echo "IDL done! Jump out!"
        break
    fi
done
python make_datacube_from_cell_spectra.py $fileseq $snapID style="$spec_style"  tauv=$tauv faceon=$rtfaceon cell_size=$cell_size with_metal=$with_metal
