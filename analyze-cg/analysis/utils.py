import mdtraj as md
import numpy as np

# imports all python files
from .molecules import *
from .angles import *
from .directors import *
from .s2 import *
from .height import *
from .load import *

def calc_all_directors(frame, masses, residues):
    """ Calculates directors for all residues in a frame. This is
    a wrapper for the calc_director function which only works for
    a single residue

    Parameters:
    -----------
    frame : mdtraj.Trajectory
        frame to analyze
    masses : list
        list of masses corresponding to each bead in the frame

    Returns:
    --------
    directors : list
        list of directors
    """
    r = [residue for residue in residues if residue.name in molecule]
    masses = np.array(masses)

    def tail_worker(residue, atom_idxs):
        """ worker function for calculating a director. This allows for
        list comprehension

        Parameters:
        -----------
        atom_idxs : list
            list of indices corresponding to a the tail

        Returns:
        --------
        director : list
            returns a single director for a tail
        """

        atoms = frame.top.select('resid {}'.format(residue.index))
        atoms = np.array(atoms).take(atom_idxs)
        res_coords = frame.xyz[0, atoms]
        res_mass = masses.take([atoms])
        com = calc_com(res_coords, res_mass)
        centered_coords = res_coords - com
        moi = calc_moi(centered_coords, res_mass)
        w, v = np.linalg.eig(moi)
        director = v[:,np.argmin(w)]
        return director

    directors = [tail_worker(residue, atom_indices) for residue in r for atom_indices in molecule[residue.name][1]]
    directors = np.array(directors)
    return directors

def calc_order_parameter(directors):
    Q = calc_q(directors)
    S2 = calc_s2(Q)
    return S2

def calc_tilt_angle(directors):
    vec1 = directors
    vec2 = [[0, 0, 1]] * len(vec1)
    tilt = calc_angle(vec1, vec2)
    return tilt
