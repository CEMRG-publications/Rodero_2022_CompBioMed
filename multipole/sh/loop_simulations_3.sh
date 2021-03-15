SA_folder="default"

for heart in 10 11;
do /home/crg17/Desktop/scripts/multipole/bin/build_HAC_table.o /data/SA_multipole/$SA_folder/h/$heart /media/crg17/"Seagate Backup Plus Drive"/CT_cases/h_case$heart/meshing/1000um/BiV;
done

# change_param=CV;
# for change_val in 0.07 0.45;
# do for heart in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
# do /home/crg17/Desktop/scripts/multipole/build_HAC_table.o /data/SA_multipole/$change_param"_"$change_val/HF/$heart /media/crg17/"Seagate Backup Plus Drive"/CT_cases/HF_case$heart/meshing/1000um/BiV;
# done;
# done;

# change_param=FEC;
# for change_val in 33 100;
# do for heart in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
# do /home/crg17/Desktop/scripts/multipole/build_HAC_table.o /data/SA_multipole/$change_param"_"$change_val/HF/$heart /media/crg17/"Seagate Backup Plus Drive"/CT_cases/HF_case$heart/meshing/1000um/BiV;
# done;
# done;


# change_param=kFEC;
# for change_val in 5 10;
# do for heart in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
# do /home/crg17/Desktop/scripts/multipole/build_HAC_table.o /data/SA_multipole/$change_param"_"$change_val/HF/$heart /media/crg17/"Seagate Backup Plus Drive"/CT_cases/HF_case$heart/meshing/1000um/BiV;
# done;
# done;


# change_param=kxf;
# for change_val in 0.6 1;
# do for heart in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
# do /home/crg17/Desktop/scripts/multipole/build_HAC_table.o /data/SA_multipole/$change_param"_"$change_val/HF/$heart /media/crg17/"Seagate Backup Plus Drive"/CT_cases/HF_case$heart/meshing/1000um/BiV;
# done;
# done;

