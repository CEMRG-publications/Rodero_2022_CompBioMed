#!/usr/bin/env python2

import os
import sys
import math
import subprocess
import copy
import numpy as np
import sets
import scipy.linalg as spla
import scipy.signal as spsig
import argparse

def read_vtx(filename):
    return np.loadtxt(filename, dtype=int, skiprows=2, ndmin=1)    

def CARPtoInit(args):

    # for a simulation with ekbatch you need an init file containing the stimuli 
    # and the CVs of the regions
    # from a CARP simulation, you will need:
    # - vtx file for a stimulus
    # - a set of tags with the conduction velocities

    # supply vtx file for stimulus
    # vtxName = '/data/Dropbox/HFCohort/13_case20/Mesh/800um/electrode_endo_rv.vtx'

    tags_myo = [1, 2]
    tags_FEC = [3, 4, 25, 26]
    tags_scar = [100]

    vtx = []
    nVtx = 0
    if args.stimulus is not None:
        for vtxFile in args.stimulus:
            temp = read_vtx(vtxFile)
            vtx.append(temp)
            nVtx += temp.shape[0]

    CV_FEC = float(args.kFEC)*float(args.CV_l)
    # ----------------------------------------------------------------------------------- #
    # write .init file
    f = open(args.output + '.init','w')

    # header
    f.write('vf:0 vs:0 vn:0 vPS:0\n') # Default properties
    f.write('retro_delay:0 antero_delay:0\n') # If there's no PS, it's ignored.
    # number of stimuli and regions
    f.write('%d %d\n' % (int(nVtx), int(len(tags_myo)) + len(tags_FEC) + len(tags_scar)))
    # stimulus
    for i in range(len(vtx)):
        if len(vtx[i]) == 1:
            f.write('%d %f\n' % (vtx[i],0))
        else:
            for n in vtx[i]:
                f.write('%d %f\n' % (int(n),0))
    # ek regions
    for i,t in enumerate(tags_myo):
        f.write('%d %f %f %f\n' % (int(t),args.CV_l,args.CV_t,args.CV_t))
    for i,t in enumerate(tags_FEC):
        f.write('%d %f %f %f\n' % (int(t),CV_FEC,CV_FEC,CV_FEC))
    for i,t in enumerate(tags_scar):
        f.write('%d %f %f %f \n' % (int(t), args.CV_scar, args.CV_scar, args.CV_scar))
    

    f.close()

def main(args):

    CARPtoInit(args)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--stimulus', 
                        action='append', 
                        help='Path and name stimulus vtx file', 
                        default=None)
    parser.add_argument('--CV_l', 
                        type=float,
                        help='give bulk conduction velocity', 
                        default=0.6)
    parser.add_argument('--CV_t', 
                        type=float,
                        help='give anisotropy ratio', 
                        default=0.35)
    parser.add_argument('--kFEC', 
                        type=float,
                        help='give Purkinje conduction velocity', 
                        default=6)
    parser.add_argument('--CV_scar', 
                        type=float,
                        help='give scar conduction velocity', 
                        default=0)
    parser.add_argument('--output',type=str, default='output',
                        help='output init file')

    args = parser.parse_args()

    main(args)

