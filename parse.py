import json
import os
import re
import sys

import numpy as np
import pandas as pa
import spglib
from scipy.constants import physical_constants

# bohr radius in Angstrom
br = physical_constants['Bohr radius'][0] / physical_constants['Angstrom star'][0]
pseudo_dir = '/scratch/SSSP_efficiency_pseudos'

SIRIUS_JSON = {
    "control": {
        "processing_unit": "gpu",
        "std_evp_solver_type": "lapack",
        "gen_evp_solver_type": "lapack",
        "verbosity": 1,
    },

    "parameters": {
        "electronic_structure_method": "pseudopotential",
        "xc_functionals": ["XC_GGA_X_PBE", "XC_GGA_C_PBE"],
        "smearing_width": 0.025,
        "use_symmetry": True,
        "num_mag_dims": 1,
        "gk_cutoff": 6.0,
        "pw_cutoff": 27.00,
        "energy_tol": 1e-8,
        "potential_tol": 1e-8,
        "num_dft_iter": 100,
        "ngridk": "None"
    },

    "iterative_solver": {
        "type": "davidson",
        "min_occupancy": 1e-5
    },

    "unit_cell": {
        "lattice_vectors": [],
        "atom_coordinate_units": "au",
        "atom_types": [""],
        "atom_files": {},
        "atoms": {
        }
    },

    "mixer": {
        "beta": 0.95,
        "type": "broyden1",
        "max_history": 8
    }
}


def irreducible_kpoints(fname_sirius_json='sirius.json'):
    sirius_config = json.load(open(fname_sirius_json, 'r'))
    pos = sirius_config['unit_cell']['atoms']
    C = np.array(sirius_config['unit_cell']['lattice_vectors'])
    positions = list(pos.values())
    indicator = [np.ones(len(x)) * i for i, x in enumerate(positions)]
    indicator = np.hstack(indicator)

    # assume atom pos given in a.u.
    pos= np.vstack(positions)
    rpos = np.linalg.solve(C.T, pos.T).T
    cell = (C, pos, indicator)
    mesh = sirius_config['parameters']['ngridk']

    # Gamma centre mesh
    mapping, grid=spglib.get_ir_reciprocal_mesh(
        mesh, cell, is_shift = [0, 0, 0])

    # Irreducible k-points
    print('ngridk              :', mesh)
    print('Number of ir-kpoints: %d' % len(np.unique(mapping)))
    # print(grid[np.unique(mapping)] / np.array(mesh, dtype=float))


def to_list(arr):
    return [list(x) for x in arr]


def load_kpoints(fname):
    data=np.loadtxt(fname, comments = '#')
    return [int(x) for x in data[:3]], list(data[3:])


def load_cell(fname):
    """
    assuming cell is given in Angstrom
    """
    return [list(x) for x in np.loadtxt(fname) / br]


if __name__ == '__main__':
    sirius_json=SIRIUS_JSON

    if len(sys.argv) > 1:
        dirname=sys.argv[1]
    else:
        dirname='./'

    pos=pa.read_csv(os.path.join(dirname, 'POS'),
                      delimiter = r'\s+', skiprows = 1, header = None)
    # store atom positions in atomic units
    pos_dict={k: to_list(np.array(data[[1, 2, 3]]/br))
                for k, data in pos.groupby(0)}

    sirius_json['unit_cell']['atom_types'] = list(pos_dict.keys())
    sirius_json['unit_cell']['atoms'] = pos_dict
    sirius_json['unit_cell']['lattice_vectors'] = load_cell(
        os.path.join(dirname, 'CELL'))
    sirius_json['unit_cell']['atom_files'] = [
        k + '.json' for k in pos_dict.keys()]
    print("natoms: ", sum([len(x) for x in pos_dict.values()]))

    ngridk, shiftk = load_kpoints(os.path.join(dirname, 'KPOINTS'))
    sirius_json['parameters']['ngridk'] = ngridk
    sirius_json['parameters']['shiftk'] = shiftk

    with open('sirius.json', 'w') as f:
        json.dump(sirius_json, f, indent=2)

    irreducible_kpoints()

    # get number of electrons, assuming UPF is called ELEM.json, where elem is the atom label
    nelectrons = 0
    elem_charge = {}
    for atom in sirius_json['unit_cell']['atom_types']:
        pseudo_json = json.load(open(os.path.join(pseudo_dir, atom+'.json')))
        z = int(pseudo_json['pseudo_potential']['header']['z_valence'])
        na = len(pos_dict[atom])
        nelectrons += z * na
        elem_charge[atom + '_' + str(na)] = z
    print('num electrons:', nelectrons)
    print('elements: ', elem_charge)
