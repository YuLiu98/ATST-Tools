import numpy as np
import sys
from ase.io import Trajectory, read, write
from ase.mep.neb import NEBTools

def neb2dimer(neb_traj:list, n_max:int=0, out_traj='dimer_init.traj', out_vec='displacement_vector.npy', 
                step_before_TS:int=1, step_after_TS:int=1):
    '''
    Transform neb chain to dimer init
    '''
    # you should use the latest neb chain or other single neb chain
    # but a full neb chain is permitted, in which the latest neb chain will be used
    if n_max == 0:
        print("=== n_max set to 0, get n_images by using NEBTools ===")
        n_images = NEBTools(neb_traj)._guess_nimages()
    elif (n_max >= 0) and (type(n_max) == int):
        n_images = n_max + 2
    else:
        raise ValueError("n_max must be a non-negative integer")
    if len(neb_traj) < n_images:
        raise ValueError(f"n_images ={n_images} is larger than the length of neb_traj ={len(neb_traj)}")
    elif len(neb_traj) > n_images:
        print(f"=== n_images ={n_images} is smaller than the length of neb_traj ={len(neb_traj)} ! Only the last {n_images} images are used ===")
    # used neb_chain is the final traj
    neb_chain = neb_traj[-n_images:]
    # get TS information from NEB chain
    barrier = NEBTools(neb_chain).get_barrier()[0]
    fmax = NEBTools(neb_chain).get_fmax()
    raw_barrier = max([image.get_potential_energy() for image in neb_chain])
    TS_info = [(ind, image) 
            for ind, image in enumerate(neb_chain) 
            if image.get_potential_energy() == raw_barrier][0]
    print(f"=== Locate TS in {TS_info[0]} of 0-{n_images-1} images  ===")
    print(f"=== TS Barrier: {barrier:.4f} (eV) ===")
    print(f"=== TS Fmax: {fmax:.4f} (eV/A) ===")
    print(f"=== TS images is output as {out_traj} for dimer init ===")
    # output TS of neb for dimer init
    write(out_traj, TS_info[1], format='traj')
    # output displancement vector by using the nearest two images of TS
    ind_before_TS = TS_info[0] - step_before_TS
    ind_after_TS = TS_info[0] + step_after_TS
    img_before = neb_chain[ind_before_TS]
    img_after = neb_chain[ind_after_TS]
    print(f"=== Displacement vector is generated by position minus from {ind_after_TS} image to {ind_before_TS} image ===")
    print("=== Notice: The displacement vector is normalized for good performance ===")
    image_vector = (img_after.positions - img_before.positions)
    vec_norm = np.linalg.norm(image_vector)
    displacement_vector = image_vector / vec_norm
    print(f"=== Displacement vector is output as {out_vec} for dimer init ===")
    np.save(out_vec,displacement_vector)
    return
    
if __name__ == "__main__":
    msg = '''
neb2dimer.py is to make dimer inputfile by using neb result
Usage: 
    python neb2dimer.py [traj_file] ([n_max])
        Notice that n_max = n_images - 2
These sctipt will output two files:
    dimer_init.traj: TS of neb chain for dimer init-structure
    displacement_vector.npy: displacement vector of dimer init
    '''
    if len(sys.argv) < 2:
        print(msg)
    elif len(sys.argv) == 2:
        if sys.argv == "--help" or sys.argv == "-h":
            print(msg)
        else:
            traj_file = sys.argv[1]
            neb_traj = Trajectory(traj_file)
            neb2dimer(neb_traj)
    else:
        traj_file = sys.argv[1]
        n_max = int(sys.argv[2])
        neb_traj = Trajectory(traj_file)
        neb2dimer(neb_traj, n_max)