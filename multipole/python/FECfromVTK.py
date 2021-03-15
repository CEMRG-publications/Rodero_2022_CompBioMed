#!/usr/bin/env python2

import os
import shutil

def mkdir(path2create, option=""):

    if option == "-p":
        if not os.path.exists(path2create):
            os.mkdir(path2create)
    else:
        shutil.rmtree(path2create)
        os.mkdir(path2create)
    

def correct_fibres_header(vtkpath, vtkname, outpath="", outname=""):

    if outpath ==  "":
        outpath = vtkpath
    if outname == "":
        outname = vtkname + "_new"

    infile = open(vtkpath + '/' + vtkname + ".vtk", 'r')
    outfile = open(outpath + '/' + outname + ".vtk", "w+")

    for line in infile.readlines():
        if line == "VECTORS fibres float\n":
            outfile.write("VECTORS fiber float\n")
        elif line == "VECTORS sheets float\n":
            outfile.write("VECTORS sheet float\n")
        elif line == "SCALARS ID int 1\n":
            outfile.write("SCALARS elemTag int 1\n")
        else:
            outfile.write(line)


    infile.close()
    outfile.close()

def extract_UVC_files(vtkpath,vtkname,nPts,nElem):

    mkdir(vtkpath + "/BiV/UVC", option="-p")

    # Header + pts + header + cells + space + header + elem + header + elem + header + fibres + header + sheets + header
    
    start_Z = 5 + nPts + 1 + nElem + 1 + 1 + nElem + 3 + nElem + 1 + nElem + 1 + nElem + 3
    start_RHO = start_Z + nPts + 2
    start_PHI = start_RHO + nPts + 2
    start_V = start_PHI + nPts + 2

    # Longitudinal coordinate

    UVC_Z=[None] * nPts
    UVC_RHO = UVC_Z
    UVC_PHI = UVC_Z
    UVC_V = UVC_Z

    UVC = UVC_Z

    count = 0
    print("starts in ")
    print(start_Z)
    with open(vtkpath + '/' + vtkname + ".vtk") as fp:
        for i, line in enumerate(fp):
            #print("starts in ")
            #print(start_Z)
            #print("\n")
            if i >= start_Z and i < start_Z + nPts:
                if '-100.' not in line:
                    print("Wrote something\n")
                    UVC[count] = line
                    count = count + 1

    with open(vtkpath + "/BiV/UVC/COORDS_Z.dat", 'w+') as f:
        for item in UVC:
            f.write("%s" % item)

    UVC_Z = UVC

    # Transmural coordinate

    count = 0

    with open(vtkpath + '/' + vtkname + ".vtk") as fp:
        for i, line in enumerate(fp):
            if i >= start_RHO and i < start_RHO + nPts:
                if '-100.' not in line:
                    UVC[count] = line
                    count = count + 1

    with open(vtkpath + "/BiV/UVC/COORDS_RHO.dat", 'w+') as f:
        for item in UVC:
            f.write("%s" % item)
    
    UVC_RHO = UVC
    # Rotational coordinate

    count = 0

    with open(vtkpath + '/' + vtkname + ".vtk") as fp:
        for i, line in enumerate(fp):
            if i >= start_PHI and i < start_PHI + nPts:
                if '-100.' not in line:
                    UVC[count] = line
                    count = count + 1

    with open(vtkpath + "/BiV/UVC/COORDS_PHI.dat", 'w+') as f:
        for item in UVC:
            f.write("%s" % item)
    
    UVC_PHI = UVC
    # Ventricular coordinate

    count = 0

    with open(vtkpath + '/' + vtkname + ".vtk") as fp:
        
            if i >= start_V and i < start_V + nPts:
                if '-100.' not in line:
                    UVC[count] = line
                    count = count + 1

    with open(vtkpath + "/BiV/UVC/COORDS_V.dat", 'w+') as f:
        for item in UVC:
            f.write("%s" % item)
    
    UVC_V = UVC
    #COMBINED 

    with open(vtkpath + "/BiV/UVC/COMBINED_COORDS_Z_RHO_PHI_V.dat", 'w+') as f:
        for i in range(len(UVC_Z)):
            f.write("%s %s %s %s\n" % (UVC_Z[i], UVC_RHO[i], UVC_PHI[i], UVC_V[i]))

if __name__ == "__main__":
    
    heart_vec = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]

    npts_vec =  [481066,691916,527572,563509,637354,539616,529708,502241,451528,445112,440291,536011,562910,674011,573966,490365,572656,438898,536860,556012,511855,537441,407578,462282]
    nelem_vec = [2349414,3490090,2636456,2873599,3266530,2656723,2598475,2489060,2227614,2176295,2169085,2635144,2788537,3423028,2866740,2360470,2807390,2192445,2663864,2712110,2544717,2679368,1991085,2257529]


    for i in [6]:
        i = 6
        heart_name = heart_vec[i]
        heart = heart_name
        heart_path = "/data/harddrive/CT_cases/HF_case" + heart_name + "/meshing/1000um"
        npts = npts_vec[i]
        nelem = nelem_vec[i]
    
        # Change the name of the headers of the fibrs and elems

        #correct_fibres_header(heart_path,heart_name)
        #"""
        # Convert from vtk to elem file

        #os.system("meshtool convert -imsh=" + heart_path + "/" + heart_name
        #        +"_new.vtk -omsh=" + heart_path + "/HF_case" + heart_name + " -ofmt=carp_txt" )
        #"""
        # Extract UVC files in the UVC folder in the BiV folder

        extract_UVC_files(heart_path,heart_name + "_new",npts,nelem)

        # Extract the BiV mesh

       # os.system("meshtool extract mesh -msh=" + heart_path + "/HF_case" + heart_name
       #+ " -tags=1,2,25,26,27,28 -submsh=" + heart_path + "/BiV/BiV -ifmt=carp_txt -ofmt=carp_txt")

        # Extract the vtx file

        #os.system("/home/crg17/Desktop/scripts/multipole/HF_elem2vtx.o " + heart)



