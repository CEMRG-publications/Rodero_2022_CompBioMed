for SA_folder in default
do

# Healthy

for heart in 10 11
do sed -i "1s/.*/lead AN AL LA PL PO/" /data/SA_multipole/$SA_folder/h/$heart/multipole_AT1090.dat
sed -i "1s/.*/lead AN AL LA PL PO/" /data/SA_multipole/$SA_folder/h/$heart/multipole_QRS.dat
sed -i "1s/.*/lead AN AL LA PL PO/" /data/SA_multipole/$SA_folder/h/$heart/multipole_TATLV.dat
done

# HF

#for heart in 01 02 03 04 05 06 07 08 09 `seq 10 24`
#do sed -i "1s/.*/lead AN AL LA PL PO/" /data/SA_multipole/$SA_folder/HF/$heart/multipole_AT1090.dat
#sed -i "1s/.*/lead AN AL LA PL PO/" /data/SA_multipole/$SA_folder/HF/$heart/multipole_QRS.dat
#sed -i "1s/.*/lead AN AL LA PL PO/" /data/SA_multipole/$SA_folder/HF/$heart/multipole_TATLV.dat
#done

done
