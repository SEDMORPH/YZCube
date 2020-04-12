#!/bin/bash
shopt -s expand_aliases #allow the alias to work
# setup idl environment
# /export/data/apps/idl85/envi53/bin/envi_setup.bash

fileseq="2xSc_07_EC_BH_vw1e4_ReposNoRFBNoRP"
snapID=151
redshift=0.04
spec_style="star_age"
rtfaceon=0
tauv=1.0
one_comp_dust=0
with_metal=1
cir_fib=1
with_PSF=1
real_fibpos=1
# step 1. get fiber location
python some_dependencies/fiber_locations.py $redshift
# step 2. get PSF weight
python some_dependencies/PSF_mass_weight.py $redshift
# step 3. create RSS(cell spectra) with circluar mock fibers or sqaure one. An example could be
# create RSS data for all 3 dithering set
for ds in $(seq 0 2)
do
    echo "creating RSS for dithering set: $ds"
    # ./create_RSS.sh $fileseq $snapID fib_radius=0.79 spec_style=$spec_style rtfaceon=$rtfaceon tauv=$tauv one_comp_dust=$one_comp_dust with_metal=$with_metal cir_fib=$cir_fib with_PSF=$with_PSF real_fibpos=$real_fibpos dither_set=$ds
done



#step 4. make datacube from RSS
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
python make_datacube_from_RSS_fits.py $fileseq $snapID style="$spec_style"  tauv=$tauv faceon=$rtfaceon with_metal=$with_metal
