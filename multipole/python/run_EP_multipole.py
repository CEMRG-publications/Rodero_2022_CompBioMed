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
    parser.add_argument('--current_case',
                        default='-1',
                        help='Case number.')
    parser.add_argument('--duration',
                        default=250.0,
                        help='Duration of experiment')
    parser.add_argument('--HF',
                        action='store_true')
    parser.add_argument('--kFEC',
                        default='7', # Adhoc to obtain reasonable AT
                        help='Conduction velocity of the FEC wrt to the myocardium.'
                        )
    parser.add_argument('--myoCV',
                        default='0.36', #Hand tuned so the HF cases match the QRS from the literature
                        help='Bulk conduction velocity of the myocardium.')
    parser.add_argument('--xfibre_aniso',
                        default='0.29')
    parser.add_argument('--FEC',
                        default='70')
    parser.add_argument('--LV_FEC_tag',
                        default='25')
    parser.add_argument('--RV_FEC_tag',
                        default='26')
    parser.add_argument('--PHI_lead',
                        choices=['AN','AL','LA','PL','PO']
                        )
    parser.add_argument('--ba2ap_lead',
                        choices=['1','2','3','4','5','6','7','8']
                        )
    parser.add_argument('--RV_electrode',
                        action='store_true')
    parser.add_argument('--RV_midseptum',
                        action='store_true')
    parser.add_argument('--folder_name',
                        default='eikonal')
    parser.add_argument('--with_scar',
                        action='store_true')
    return parser


def jobID(args):
    """
    Generate name of top level output directory.
    """
    if not args.HF:
        outPath	 = '/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case' + args.current_case + '/simulations/multipole/' + args.folder_name
    if args.HF:
        outPath = '/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case' + args.current_case +'/simulations/multipole/' + args.folder_name

    cmd = 'mkdir -p ' + outPath
    os.system(cmd)

    if args.RV_electrode:
        outName = outPath + '/BiV.RV_endo.apex'
    elif args.RV_midseptum:
        outName = outPath + '/BiV.midseptum'
    else:
        outName = outPath + '/{}'.format(args.PHI_lead + '_' + args.ba2ap_lead)

    return outName

@tools.carpexample(parser, jobID)

def run(args, job):

    if not args.HF:
        meshdir	 = '/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case' + args.current_case + '/meshing/1000um/BiV'
        meshdir4bash = '/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases/h_case' + args.current_case + '/meshing/1000um/BiV'
    if args.HF:
        meshdir  = '/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case' + args.current_case + '/meshing/1000um/BiV'
        meshdir4bash = '/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases/HF_case' + args.current_case + '/meshing/1000um/BiV'
    meshname   = '{}/BiV'.format(meshdir)

    cmd = 'cp ' + meshdir4bash + '/FEC/BiV_FEC_w5_h' + args.FEC + '.elem ' + meshdir4bash + '/BiV.elem'
    os.system(cmd)
    cmd = 'cp ' + meshdir4bash + '/FEC/BiV_FEC_w5_h' + args.FEC + '_retagged.elem ' + meshdir4bash + '/BiV.elem'
    os.system(cmd)

    if args.with_scar:
        cmd = 'cp ' + meshdir4bash + '/BiV_FEC70_scar5mm.elem ' + meshdir4bash + '/BiV.elem'
        os.system(cmd)

    cmd  = tools.carp_cmd()
    
    cmd += ['-simID', job.ID]

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

    cmd += ['-tend', args.duration, # length of simulation (ms)
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
    
    num_regs = 2 # Myocardia and FECs

    if args.with_scar:
        num_regs += 1

    gregion = ['-num_gregions', num_regs]

    gregion += ['-gregion[0].name', 'Ventricular myocardia',
                '-gregion[0].num_IDs', 2,
                '-gregion[0].ID[0]', 1,
                '-gregion[0].ID[1]', 2]
    
    gregion += ['-gregion[1].name', 'Fast endocardial conduction layer',
                '-gregion[1].num_IDs', 2,
                '-gregion[1].ID[0]', 25,
                '-gregion[1].ID[1]', 26]
    
    if args.with_scar:
        gregion += ['-gregion[2].name', 'Scar',
                    '-gregion[2].num_IDs', 1,
                    '-gregion[2].ID[0]', 100]

    ekregion = ['-num_ekregions', num_regs]

    ekregion += ['-ekregion[0].ID', 0,
                 '-ekregion[0].vel_f', args.myoCV,
                 '-ekregion[0].vel_s',  float(args.xfibre_aniso)*float(args.myoCV),#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4608725
                 '-ekregion[0].vel_n',  float(args.xfibre_aniso)*float(args.myoCV)]
    
    ekregion += ['-ekregion[1].ID', 1,
                 '-ekregion[1].vel_f', float(args.kFEC)*float(args.myoCV),
                 '-ekregion[1].vel_s', float(args.kFEC)*float(args.myoCV),
                 '-ekregion[1].vel_n', float(args.kFEC)*float(args.myoCV)]
                 
    if args.with_scar:
        ekregion += ['-ekregion[2].ID', 2,
                    '-ekregion[2].vel_f', 0.0,
                    '-ekregion[2].vel_s', 0.0,
                    '-ekregion[2].vel_n', 0.0]


    return gregion + ekregion

# --- set cellular electrical dynamics ---------------------------------------

def impRegions(args):

    num_imp_regions = 1

    imp = ['-num_imp_regions', num_imp_regions] 

    # Ventricles 
    imp += ['-imp_region[0].im',      'TT2',
            '-imp_region[0].num_IDs', 4 + float(args.with_scar),
            '-imp_region[0].ID[0]',   1,
            '-imp_region[0].ID[1]',   2,
            '-imp_region[0].ID[2]',   25,
            '-imp_region[0].ID[3]',   26]

    if args.with_scar:
        imp += ['-imp_region[0].ID[4]', 100]

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

    if args.HF:
        leadsPath = '/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case' + args.current_case +'/simulations/multipole'
    else:
        leadsPath = '/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case' + args.current_case + '/simulations/multipole'

    stim = ['-stimulus[0].stimtype', 0,
            '-stimulus[0].start',    0.0,
            '-stimulus[0].strength', Istim,
            '-stimulus[0].duration', Dstim]
    
    if args.RV_electrode:
        stim += ['-stimulus[0].vtx_file', leadsPath + "/BiV.RV_endo.apex"]
    elif args.RV_midseptum:
        stim += ['-stimulus[0].vtx_file', leadsPath + "/BiV.midseptum"]
    else:
        stim += ['-stimulus[0].vtx_file', leadsPath + "/" + args.PHI_lead + '_' + args.ba2ap_lead]
            
    return stim
# --- set stimulus -----------------------------------------------------------

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        scripthelp()
    else:
        run()
