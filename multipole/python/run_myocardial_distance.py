#!/usr/bin/python

EXAMPLE_DESCRIPTIVE_NAME = 'Multipole simulation' 
EXAMPLE_AUTHOR = ('Cristobal Rodero')

import os
from datetime import date
from glob import glob

from carputils import settings
from carputils import tools
from carputils import ep
from carputils import model
#from carputils import testing
from carputils.resources import petsc_block_options

import sys
#~ MPI_EXEC       ='mpiexec'


def scripthelp():
    print('Script to run an EP simulation.')
    print('Mandatory arguments: --current_case, --PHI_lead, --ba2ap_lead')
    print('Optional arguments: --mesh_size, --duration, --HF, --kFEC, --myoCV, --FEC, --LV_FEC_tag, --RV_FEC_tag')
    print('Use the help of each argument for more information.')

def parser():
    parser = tools.standard_parser()
    parser.add_argument('--path',
                        default='',
                        help='Path to the mesh.')
    parser.add_argument('--simulation_name',
                        default='',
                        help='Name of the simulation.')
    return parser


def jobID(args):
    """
    Generate name of top level output directory.
    """
    return '{}/{}'.format(args.path,args.simulation_name)

@tools.carpexample(parser, jobID)

def run(args, job):


    meshname   = '{}/BiV'.format(args.path)

    
    cmd = ['-simID', job.ID]

    cmd += ['-meshname', meshname]
    
    stims = setupStimuli(args)
    propOpts, numStims, stims = setupPropagation('R-E+', 1, stims)   
    """
      Returns
    -------
    propOpts: list
        List of CARP command line options to set up propagation model
    numStims: int
        updated number of stimuli
    stims: list
        list of stimuli commands needed for setting up the propagation model
    """

    cmd += ['-bidomain',0]  # It's solving monodomain.

    stimOpts = [ '-num_stim', numStims ] + stims
    cmd += stimOpts 

    actTime = activationTime(args)
    cmd += actTime

    conduction = conductivity(args)
    cmd += conduction

    imp = impRegions(args)
    cmd += imp
    
    # When to save:

    cmd += ['-tend', 200, # length of simulation (ms)
            '-dt', 100., # Time step size (us)
            '-tsav', 1.] # times to save state (ms)


    cmd += ['-experiment',6] # Eikonal Solve, so only depolarization.

    # Run main CARP simulation
    # ------------------------
    job.carp(cmd)

# ============================================================================
#    EP FUNCTIONS
# ============================================================================

def setupPropagation(propagation, numStims, stims):

    propOpts, numStims, stims = ep.setupPropagation(propagation, numStims, stims, '')

    return propOpts, numStims, stims

def setupSolverEP(args):

    sopts = [ '-parab_solve',        1,
              '-pstrat',             1,
              '-pstrat_i',           1,
              '-pstrat_backpermute', 1,
              '-mapping_mode',       1,
              '-redist_extracted',   0,
              '-mass_lumping',       1 ]

    return sopts

# --- set conduction velocity and conductivities -----------------------------

def conductivity(args):  
    
    num_regs = 1 

    gregion = ['-num_gregions', num_regs]

    gregion += ['-gregion[0].name', 'Ventricular myocardia',
                '-gregion[0].num_IDs', 4,
                '-gregion[0].ID[0]', 1,
                '-gregion[0].ID[1]', 2,
                '-gregion[0].ID[2]', 25,
                '-gregion[0].ID[3]', 26]
    

    ekregion = ['-num_ekregions', num_regs]

    ekregion += ['-ekregion[0].ID', 0,
                 '-ekregion[0].vel_f', 1,
                 '-ekregion[0].vel_s',  1,
                 '-ekregion[0].vel_n',  1]
    


    return gregion + ekregion

# --- set cellular electrical dynamics ---------------------------------------

def impRegions(args):

    num_imp_regions = 1

    imp = ['-num_imp_regions', num_imp_regions] 

    # Ventricles 
    imp += ['-imp_region[0].im',      'TT2',
            '-imp_region[0].num_IDs', 4,
            '-imp_region[0].ID[0]',   1,
            '-imp_region[0].ID[1]',   2,
            '-imp_region[0].ID[2]',   25,
            '-imp_region[0].ID[3]',   26]
    return imp

# --- set activation time computation ----------------------------------------
def activationTime(args):

    activation = ['-num_LATs',            1,
                  '-lats[0].measurand',   0,      # Vm (0) or extracellular potential (1)
                  '-lats[0].all',         0,      # All activation times (1) or only the 1st (0)
                  '-lats[0].threshold', -60,
                  '-lats[0].method',      1 ]     # Threshold crossing

    return activation

# --- set stimulus -----------------------------------------------------------

def setupStimuli(args):

    # Generate electrical trigger options
    Istim    = 150  # strength
    Dstim    = 2    # duration
    num_stims = 1


    stim = ['-stimulus[0].stimtype', 0,
            '-stimulus[0].start',    0.0,
            '-stimulus[0].strength', Istim,
            '-stimulus[0].duration', Dstim]
    
    stim += ['-stimulus[0].vtx_file', '{}/BiV.epi.surf.vtx'.format(args.path)]
            
    return stim
# --- set stimulus -----------------------------------------------------------

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        scripthelp()
    else:
        run()
