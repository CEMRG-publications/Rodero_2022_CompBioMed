## ---- Default simulations

# Healthy

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 20`;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_default --np 20 --overwrite-behaviour overwrite;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --RV_electrode --folder_name eikonal_default --np 20 --overwrite-behaviour overwrite;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished default healthy"

# HF

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_default --overwrite-behaviour overwrite --np 20;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --RV_electrode --folder_name eikonal_default --overwrite-behaviour overwrite --np 20;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished default HF"

# ### Changing FEC layer height

# healthy

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 20`;
do for FEC_values in 33 100;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --FEC $FEC_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_FEC_$FEC_values --np 20 --overwrite-behaviour overwrite;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --FEC $FEC_values --RV_electrode --folder_name eikonal_FEC_$FEC_values --overwrite-behaviour overwrite --np 20;
done;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished FEC layer healthy"


# --- HF

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`;
do for FEC_values in 33 100;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --FEC $FEC_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_FEC_$FEC_values --np 20 --overwrite-behaviour overwrite;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --FEC $FEC_values --RV_electrode --folder_name eikonal_FEC_$FEC_values --overwrite-behaviour overwrite --np 20;
done;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished FEC layer HF"

### Changing CV - fibre direction

# --- Healthy

# for heart in 01 02 03 04 05 06 07 08 09 `seq 10 20`;
# do for CV_values in 0.07 0.8;
# do for vein in AN AL LA PL PO;
# do for elec_choice in `seq 1 8`;
# do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --myoCV $CV_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_CV_$CV_values --overwrite-behaviour overwrite --np 20;
# done;
# done;
# /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --myoCV $CV_values --RV_electrode --folder_name eikonal_CV_$CV_values --overwrite-behaviour overwrite --np 20;
# done;
# done

# --- HF

# for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`;
# do for CV_values in 0.07 0.8;
# do for vein in AN AL LA PL PO;
# do for elec_choice in `seq 1 8`;
# do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --myoCV $CV_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_CV_$CV_values --overwrite-behaviour overwrite --np 20;
# done;
# done;
# /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --myoCV $CV_values --RV_electrode --folder_name eikonal_CV_$CV_values --overwrite-behaviour overwrite --np 20;
# done;
# done


### Changing cross-fibre anisotropy

# --- Healthy

# for heart in 01 02 03 04 05 06 07 08 09 `seq 10 20`;
# do for kxf_values in 0.6 1;
# do for vein in AN AL LA PL PO;
# do for elec_choice in `seq 1 8`;
# do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --xfibre_aniso $kxf_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_kxf_$kxf_values --overwrite-behaviour overwrite --np 20;
# done;
# done;
# /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --xfibre_aniso $kxf_values --RV_electrode --folder_name eikonal_kxf_$kxf_values --overwrite-behaviour overwrite --np 20;
# done;
# done

# # --- HF

# for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`;
# do for kxf_values in 0.6 1;
# do for vein in AN AL LA PL PO;
# do for elec_choice in `seq 1 8`;
# do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --xfibre_aniso $kxf_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_CV_$CV_values --overwrite-behaviour overwrite --np 20;
# done;
# done;
# /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --xfibre_aniso $kxf_values --RV_electrode --folder_name eikonal_CV_$CV_values --overwrite-behaviour overwrite --np 20;
# done;
# done


# ### Changing endocardial anisotropy

# # --- Healthy

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 20`;
do for kFEC_values in 3 5 10;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --kFEC $kFEC_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_kFEC_$kFEC_values --overwrite-behaviour overwrite --np 20;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --kFEC $kFEC_values --RV_electrode --folder_name eikonal_kFEC_$kFEC_values --overwrite-behaviour overwrite --np 20;
done;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished kFEC healthy"

# --- HF

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`;
do for kFEC_values in 3 5 10;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --kFEC $kFEC_values --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_kFEC_$kFEC_values --overwrite-behaviour overwrite --np 20;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --kFEC $kFEC_values --RV_electrode --folder_name eikonal_kFEC_$kFEC_values --overwrite-behaviour overwrite --np 20;
done;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished kFEC HF"


# # --- With scar

# # --- Healthy

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 20`;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --with_scar --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_scar --overwrite-behaviour overwrite --np 20;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --with_scar --RV_electrode --folder_name eikonal_scat --overwrite-behaviour overwrite --np 20;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished scar healthy"

# # --- HF

for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`;
do for vein in AN AL LA PL PO;
do for elec_choice in `seq 1 8`;
do /home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --with_scar --PHI_lead $vein --ba2ap_lead $elec_choice --folder_name eikonal_scar --overwrite-behaviour overwrite --np 20;
done;
done;
/home/crg17/Desktop/scripts/multipole/python/run_EP_multipole.py --current_case $heart --HF --with_scar --RV_electrode --folder_name eikonal_scar --overwrite-behaviour overwrite --np 20;
done

/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Finished scar HF"
