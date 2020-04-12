#!/bin/bash
. ~/.bashrc
# setup idl environment
# . /export/data/apps/idl85/envi53/bin/envi_setup.bash
# export PATH=$PATH:/export/data/apps/idl85/idl85/bin
# shopt -s expand_aliases # allow the alias to work
# declare -a xlist=(-5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0)
# declare -a ylist=(-5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0)

fileseq=$1
snapID=$2
for ARGUMENT in "$@"
do

    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=)

    case "$KEY" in
            cell_size)              cell_size=${VALUE} ;;
	    fib_radius)              fib_radius=${VALUE} ;;
	    redshift)              redshift=${VALUE} ;;
            spec_style)              spec_style=${VALUE} ;;
            rtfaceon)    rtfaceon=${VALUE} ;;
            one_comp_dust)    one_comp_dust=${VALUE} ;;
            with_metal)    with_metal=${VALUE} ;;
            with_PSF)    with_PSF=${VALUE} ;;
            tauv)    tauv=${VALUE} ;;
            cir_fib)    cir_fib=${VALUE} ;; #using Circular fibers
            real_fibpos)    real_fibpos=${VALUE} ;; #using real fibers position
            dither_set)    dither_set=${VALUE} ;; # which dither_set is used (0 or 1 or 2), set a negative value if square fibers are used
            *)
    esac


done

threadnum=12
sleep_time=5

option_str="/cell_spectra"
if [ ! -z ${cell_size+"nothing"} ]
then
    option_str="$option_str, cell_size=$cell_size"
fi
if [ ! -z ${fib_radius+"nothing"} ]
then
    option_str="$option_str, fib_radius=$fib_radius"
fi
if [ ! -z ${redshift+"nothing"} ]
then
    option_str="$option_str, redshift=$redshift"
fi
if [ ! -z ${tauv+"nothing"} ]
then
    option_str="$option_str, tauv=$tauv"
else
    echo "default tauv"
fi

if [ $spec_style = "SEDmorph" ]
then
    option_str="$option_str"
else
    option_str="$option_str, spec_style='$spec_style'"
fi
if [ $rtfaceon -gt 0 ]
then
    option_str="$option_str,/rtfaceon"
fi
if [ $one_comp_dust -gt 0 ]
then
    option_str="$option_str,/one_comp_dust"
fi
if [ $with_metal -gt 0 ]
then
    option_str="$option_str,/with_metal"
fi
if [ $with_PSF -gt 0 ]
then
    option_str="$option_str,/with_PSF"
fi
if [ $cir_fib -gt 0 ]
then
    option_str="$option_str,/cir_fib"
fi
echo $option_str


if [ $real_fibpos -gt 0 ]
then
  xlist=`python Get_xylist.py get_x=True real_fibpos=True dither_set=$dither_set`
  ylist=`python Get_xylist.py get_y=True real_fibpos=True dither_set=$dither_set`
  declare -a xlist=($xlist)
  declare -a ylist=($ylist)
else
  xlist=`python Get_xylist.py cell_size=$cell_size xymin=-5 xymax=5`
  ylist=$xlist
fi

# echo ${xlist[0]}
if [ $real_fibpos -gt 0 ]
then
  echo "working"
  n_fib=${#xlist[@]}
  for i in $(seq 0 $n_fib)
  do
    x=${xlist[$i]}
    y=${ylist[$i]}
    echo "cell_pos(x,y):" $x, $y
    command_str="sedm2_run, '$fileseq',snap=$snapID, cell_x_offset=$x,cell_y_offset=$y,$option_str"
    echo $command_str
    nohup idl -e "$command_str" &
    # idl -e "$command_str"
    #keep only 10 threads for idl, do not use too much threads
    idl_count="$(ps aux|grep idl|wc -l)"
    while [ $idl_count -ge $threadnum ]
    do
	echo $i "done   " $idl_count "ongoing"
        sleep $sleep_time
        idl_count="$(ps aux|grep idl|wc -l)"
    done
  done

else
  for x in $xlist
  # for x in "${xlist[@]}"
  do
      for y in $ylist
      # for y in "${ylist[@]}"
      do
      	echo "cell_pos(x,y):" $x, $y
      	command_str="sedm2_run, '$fileseq',snap=$snapID, cell_x_offset=$x,cell_y_offset=$y, $option_str"
      	echo $command_str
      	nohup idl -e "$command_str" &
      	# idl -e "$command_str"


      	#keep only 10 threads for idl, do not use too much threads
      	idl_count="$(ps aux|grep idl|wc -l)"
      	while [ $idl_count -ge $threadnum ]
      	do
      	    sleep $sleep_time
      	    idl_count="$(ps aux|grep idl|wc -l)"
      	done # while
      done # for y
  done # for x
fi # if real_fibpos


