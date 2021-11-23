import glob
import numpy as np
import os
import pathlib
import shutil

import  files_manipulations

def copy_and_scale_4ch_lv(cohort="h", heart_case="01", factor=1.0):

    old_mesh_dir = os.path.join("/media","crg17","Seagate Backup Plus Drive","CT_cases",cohort + "_case" + heart_case,"meshing","1000um")
    old_mesh_dir_slashes = os.path.join("/media","crg17","Seagate\ Backup\ Plus\ Drive","CT_cases",cohort + "_case" + heart_case,"meshing","1000um")

    new_mesh_dir = os.path.join("/media","crg17","Seagate Backup Plus Drive","CT_cases",cohort + "_case" + heart_case + "_" + str(factor),"meshing","1000um","cavities")
    new_mesh_dir_slashes = os.path.join("/media","crg17","Seagate\ Backup\ Plus\ Drive","CT_cases",cohort + "_case" + heart_case + "_" + str(factor),"meshing","1000um")

    if cohort == "h":
        old_mesh_dir += "/cavities"

    pathlib.Path(new_mesh_dir).mkdir(parents=True, exist_ok=True)

    os.system(os.path.join("/home", "common", "cm2carp", "bin", "return_carp2original_coord.pl ") +
              os.path.join(old_mesh_dir_slashes, cohort + "_case" + heart_case) + ".pts " + str(factor) + " 0 0 0 > " +
              os.path.join(new_mesh_dir_slashes, cohort + "_case" + heart_case) + "_" + str(factor) + ".pts")

    shutil.copy(os.path.join(old_mesh_dir, "LV_endo_closed") + ".surf",
                os.path.join(new_mesh_dir, "LV_endo_closed") + ".surf")

def copy_and_scale_biv(cohort="h", heart_case="all", factor=1.25, old_suffix="noPVTV", new_suffix="big"):

    if heart_case == "all":
        heart_cases = ["0" + str(i) for i in range(1, 10)] + [str(i) for i in range(10, 21)]
        if cohort == "HF":
            heart_cases += ["21", "22", "23", "24"]
        for heart_case_i in heart_cases:
            copy_and_scale_biv(cohort=cohort,heart_case=heart_case_i,factor=factor,old_suffix=old_suffix, new_suffix=new_suffix)
    else:
        print("Scaling " + cohort + " " + heart_case)
        mesh_dir = os.path.join("/media","crg17","Seagate Backup Plus Drive","CT_cases",cohort + "_case" + heart_case,
                                "meshing","1000um","BiV","meshes")
        mesh_dir_slashes = os.path.join("/media","crg17","Seagate\ Backup\ Plus\ Drive","CT_cases",
                                        cohort + "_case" + heart_case,"meshing","1000um","BiV","meshes")

        os.system(os.path.join("/home", "common", "cm2carp", "bin", "return_carp2original_coord.pl ") +
                  os.path.join(mesh_dir_slashes, "BiV_FEC_w5_h33_retagged_" + old_suffix + ".pts") + " " +
                  str(factor) + " 0 0 0 > " +
                  os.path.join(mesh_dir_slashes, "BiV_FEC_w5_h33_retagged_" + new_suffix + ".pts"))

        print("Writing " + os.path.join(mesh_dir, "BiV_FEC_w5_h33_retagged_" + new_suffix + ".elem"))

        shutil.copy(os.path.join(mesh_dir, "BiV_FEC_w5_h33_retagged_" + old_suffix + ".elem"),
                    os.path.join(mesh_dir, "BiV_FEC_w5_h33_retagged_" + new_suffix + ".elem"))

        shutil.copy(os.path.join(mesh_dir, "BiV_FEC_w5_h33_retagged_" + old_suffix + ".lon"),
                    os.path.join(mesh_dir, "BiV_FEC_w5_h33_retagged_" + new_suffix + ".lon"))

def vol_3dstokes(pts_file, surf_file):
    """Function to calculate the area of a surface or the volume of a closed
    surface.

    Args:
        pts_file (pts): Pts object corresponding to the coordinates of the
        surface.
        surf_file (surf): surface object to calculate the area or volume from.
    Returns:
        float : vol.
    """

    d12 = np.array([np.nan, np.nan, np.nan])
    d13 = np.array([np.nan, np.nan, np.nan])
    vol = np.array([])

    for i in range(surf_file.size):
        d13[0] = pts_file.p1[surf_file.i1[i]] - \
                 pts_file.p1[surf_file.i3[i]]
        d13[1] = pts_file.p2[surf_file.i1[i]] - \
                 pts_file.p2[surf_file.i3[i]]
        d13[2] = pts_file.p3[surf_file.i1[i]] - \
                 pts_file.p3[surf_file.i3[i]]

        d12[0] = pts_file.p1[surf_file.i1[i]] - \
                 pts_file.p1[surf_file.i2[i]]
        d12[1] = pts_file.p2[surf_file.i1[i]] - \
                 pts_file.p2[surf_file.i2[i]]
        d12[2] = pts_file.p3[surf_file.i1[i]] - \
                 pts_file.p3[surf_file.i2[i]]

        cr = np.cross(d13, d12)
        crNorm = np.linalg.norm(cr)

        area_tr = 0.5 * crNorm

        zMean = (pts_file.p3[surf_file.i1[i]] + \
                 pts_file.p3[surf_file.i2[i]] + \
                 pts_file.p3[surf_file.i3[i]]) / 3.

        nz = cr[2] / crNorm

        vol = np.append(vol, area_tr * zMean * nz)

    return vol

def compute_volume(cohort="h", heart_case="01", factor=1.0):

    heart_name = cohort + "_case" + heart_case + "_" + str(factor)
    path2fourch = os.path.join("/media","crg17","Seagate Backup Plus Drive","CT_cases",heart_name,"meshing","1000um")

    copy_and_scale_biv(cohort=cohort, heart_case=heart_case, factor=factor)

    lvendo_closed_pts = files_manipulations.pts.read(os.path.join(path2fourch, heart_name + ".pts"))
    lvendo_closed_surf = files_manipulations.surf.read(os.path.join(path2fourch, "cavities", "LV_endo_closed.surf"))

    temp_vol = vol_3dstokes(pts_file=lvendo_closed_pts, surf_file=lvendo_closed_surf)

    LV_chamber_vol = np.abs(sum(temp_vol) * 1e-12)  # In mL

    return round(LV_chamber_vol, 2)

def show_mean_volume(cohort="h", factor=1.0):

    heart_cases = ["0" + str(i) for i in range(1,10)] + [str(i) for i in range(10,21)]
    if cohort == "HF":
        heart_cases += ["21", "22", "23", "24"]
    volume_list = []

    for heart_case in heart_cases:
        heart_vol = compute_volume(cohort=cohort, heart_case=heart_case, factor=factor)
        volume_list.append(heart_vol)

    print(np.mean(volume_list))
    print(np.std(volume_list))