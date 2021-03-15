# CV=0.47 --> QRS = 132 ± 19
# CV=0.4 --> QRS = 203 ± 19 
CV_value=0.36
for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`
do
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --myoCV $CV_value --RV_electrode --folder_name eikonal_default_CV_$CV_value --np 20 --overwrite-behaviour overwrite 
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished CV 0.36. Now crop base!"
