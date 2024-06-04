# YuLiu98 in 2024-06-03
# Run AutoNEB calculation on NEB images by FCP-ABACUS-ASE
# part of ATST-Tools scripts

from ase.optimize import FIRE, BFGS
from ase.io import read, write
from ase.parallel import world, parprint, paropen
#from pathlib import Path
from abacus_autoneb import FcpAutoNEB

# neb setting
mpi = 16
omp = 4
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "AutoNEBrun"
# algorism = 'eb' # default for AutoNEB
algorism = "improvedtangent" # IT-NEB
init_chain = "init_neb_chain.traj"
climb = True
fmax = [0.20, 0.05]  # eV / Ang, 2 fmax for others and last CI-NEB in all-images
n_simul = world.size # only for autoneb, number of simultaneous calculation
n_images = 10 # only for autoneb, max number of all image, which should be reached
smooth_curve = False # True to do more neb step to smooth the curve.
k = 0.10 # eV/Ang^2, force constant of spring, 0.05 is from VTST-Tools

# abacus parameters
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}/"
basis_dir = f"{lib_dir}/"
pp = {
      'C':'C_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      'Pt':'Pt_ONCV_PBE-1.0.upf',
      }
basis = {
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb'
         ,}
kpts = [2, 1, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'dft_functional': 'pbe',
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'd3_bj',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-6,
    'scf_nmax': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'cal_force': 1,
    'cal_stress': 1,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
    'efield_pos_max': 0.0
}

# FCP parameters
fcp_parameters = { # innercalc = cal_abacus,  # innercalc must NOT be set, since it is automatically set to abacus
    'fcptxt': 'log-fcp.txt',
    'U': 1.5,
    'NELECT': 246,
    'C': 1/80,
    'FCPmethod': 'Newton-fitting',
    'FCPconv': 0.01,
    'NELECT0': 246, 
    'adaptive_lr': False,
    'work_ref': 4.6,
    'max_FCP_iter': 10000
}


if __name__ == "__main__": 
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = FcpAutoNEB(init_chain, fcp_parameters=fcp_parameters, 
                     parameters=parameters, algorism=algorism, 
                     directory=neb_directory, k=k,
                     n_simul=n_simul, n_max=n_images, 
                     abacus=abacus,  mpi=mpi, omp=omp, )
    neb.run(optimizer=neb_optimizer, climb=climb, 
                fmax=fmax, smooth_curve=smooth_curve)