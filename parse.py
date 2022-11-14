import json
import os
import re
import sys

import re
import numpy as np
import pandas as pa
import spglib
from scipy.constants import physical_constants
from readers import read_system, read_electrons, read_species, to_float

# bohr radius in Angstrom
br = physical_constants['Bohr radius'][0] * 1e10

# pseudo
pseudo_dir = '/home/simon/work/sg15_hgh'

# convert QE smearing name to sirius name
smearing_qe2sirius = {'gauss': 'gaussian',
                      'gauss': 'gaussian',
                      'f-d': 'fermi_dirac',
                      'fermi-dirac': 'fermi_dirac',
                      'fd': 'fermi_dirac',
                      'marzari-vanderbilt': 'cold',
                      'm-v': 'cold',
                      'mv': 'cold',
                      'cold': 'cold',
                      'methfessel-paxton': 'methfessel_paxton',
                      'm-p': 'methfessel_paxton',
                      'mp': 'methfessel_paxton'}

SIRIUS_JSON = {
    "control": {
        "processing_unit": "cpu",
        "std_evp_solver_type": "lapack",
        "gen_evp_solver_type": "lapack",
        "verbosity": 1,
    },
    "parameters": {
        "electronic_structure_method": "pseudopotential",
        "xc_functionals": ["XC_GGA_X_PBE", "XC_GGA_C_PBE"],
        "smearing_width": np.nan,
        "smearing": "TOOD",
        "use_symmetry": True,
        "num_mag_dims": '<MISSING>',
        "gk_cutoff": np.nan,
        "pw_cutoff": np.nan,
        "num_dft_iter": 40,
    },
    "unit_cell": {
        "lattice_vectors": [],
        "atom_coordinate_units": "au",
        "atom_types": [""],
        "atom_files": {},
        "atoms": {
        }
    },
}


def irreducible_kpoints(fname_sirius_json='sirius.json'):
    """
    TODO: this function ignores magnetic states
    """
    sirius_config = json.load(open(fname_sirius_json, 'r'))
    pos = sirius_config['unit_cell']['atoms']
    C = np.array(sirius_config['unit_cell']['lattice_vectors'])
    if 'ngridk' not in sirius_config['parameters']:
        return

    # positions is a list of lists of atom positions, the first dimension
    # is for the number of atom types
    positions = list(pos.values())
    indicator = [np.ones(len(x)) * i for i, x in enumerate(positions)]
    indicator = np.hstack(indicator)

    # assume atom pos given in a.u.
    pos = np.vstack([np.array(x)[:, :3] for x in positions])
    rpos = np.linalg.solve(C.T, pos[:,:3].T).T
    cell = (C, rpos, indicator)

    mesh = sirius_config['parameters']['ngridk']
    if (rpos > 1).any() or (rpos < 0).any():
        print('WARNING: atoms not inside unit cell')

    # Gamma centre mesh
    mapping, grid = spglib.get_ir_reciprocal_mesh(
        mesh, cell, is_shift=[0, 0, 0])

    # Irreducible k-points
    print('ngridk              :', mesh)
    print('Number of ir-kpoints: %d' % len(np.unique(mapping)))
    # print(grid[np.unique(mapping)] / np.array(mesh, dtype=float))


def to_list(arr):
    """to_list"""
    return [list(x) for x in arr]


def load_kpoints(fname):
    """load_kpoints"""
    with open(fname, 'r', encoding='ascii') as fi:
        fbuf = fi.read()
        lines = fbuf.split('\n')
    l0 = lines[0]

    if re.match('# crystal', l0):
        nk = np.loadtxt(lines[1:2])
        kpoints_array = np.loadtxt(lines[2:]).reshape(-1, 4)
        # stripping the weights, is thisk ok?
        kpoints = [x[0:3] for x in kpoints_array]
        assert(len(kpoints) == nk)
        return {'vk': kpoints}
    elif re.match('# automatic', l0):
        # return kpoint mesh
        data = np.loadtxt(lines[1:])
        return {'ngridk': [int(x) for x in data[:3]],
                'shiftk': list(data[3:])}
    else:
        raise Exception('unknown k-point specifier: ' + l0)


def load_cell(fname):
    """
    assuming cell is given in Angstrom
    """
    return [list(x) for x in np.loadtxt(fname) / br]


if __name__ == '__main__':
    sirius_json = SIRIUS_JSON

    if len(sys.argv) > 1:
        dirname = sys.argv[1]
    else:
        dirname = './'

    qe_system = read_system(os.path.join(dirname, 'SYSTEM'))
    qe_electrons = read_electrons(os.path.join(dirname, 'ELECTRONS'))
    qe_species = read_species(os.path.join(dirname, 'SPECIES'))
    print('qe_system[\'magnetization\']', qe_system['magnetization'])
    # import ipdb; ipdb.set_trace()

    qe_magnetization = {}
    # combine qe_system and qe_species and store magnetization vector
    for i, data  in enumerate(qe_species):
        # magnetization in z-direction
        if qe_system['magnetization']:
            mag = qe_system['magnetization'][i]
            qe_magnetization[data['type']] = np.array([0, 0, np.sign(mag)])

    # store Atom postions into pandas table (first row = atom type)
    pos = pa.read_csv(os.path.join(dirname, 'POS'),
                      delimiter=r'\s+', skiprows=1, header=None)
    print(pos.dtypes)
    # remove number from atom
    pos[0] = pos[0].apply(lambda x: re.match('[A-Za-z]+', x).group(0))
    print('---pos---')
    print(pos)

    # convert atom positions to atomic units
    pos_dict = {k: to_list(np.array(data[[1, 2, 3]]/br))
                for k, data in pos.groupby(0)}

    if 'magnetization' in qe_system:
        # extend by magnetization
        for atom_type in pos_dict:
            if atom_type in qe_magnetization:
                lpos = [x+list(qe_magnetization[atom_type]) for x in pos_dict[atom_type]]
                pos_dict[atom_type] = lpos
            else:
                # this atom does not have magnetization
                continue
    sirius_json['parameters']['smearing'] = smearing_qe2sirius[qe_system['smearing']]
    sirius_json['parameters']['smearing_width'] = 0.5 * qe_system['degauss']  # rydberg -> hartree
    print('sirius params', sirius_json['parameters'])

    sirius_json['unit_cell']['atom_types'] = list(pos_dict.keys())

    sirius_json['unit_cell']['atoms'] = pos_dict

    sirius_json['unit_cell']['qe_magnetization'] = {x: list(qe_magnetization[x]) for x in qe_magnetization} # just for the record

    sirius_json['unit_cell']['lattice_vectors'] = load_cell(os.path.join(dirname, 'CELL'))

    sirius_json['unit_cell']['atom_files'] = {k: k + '.json' for k in pos_dict.keys()}

    print("natoms: ", sum([len(x) for x in pos_dict.values()]))
    volume = np.abs(np.linalg.det(np.array(sirius_json['unit_cell']['lattice_vectors'])))
    print("volume: %.2f" % volume)

    kpoints_data = load_kpoints(os.path.join(dirname, 'KPOINTS'))
    if 'ngridk' in kpoints_data:
        sirius_json['parameters']['ngridk'] = kpoints_data['ngridk']
        sirius_json['parameters']['shiftk'] = kpoints_data['shiftk']
    else:
        sirius_json['parameters']['vk'] = [list(x[:3]) for x in kpoints_data['vk']]
        # list of kpoints

    #
    if qe_system['nspin'] == 2:
        sirius_json['parameters']['num_mag_dims'] = 1
    elif qe_system['nspin'] == 1:
        sirius_json['parameters']['num_mag_dims'] = 0
    else:
        # default, non-magnetic
        sirius_json['parameters']['num_mag_dims'] = 0
    sirius_json['parameters']['pw_cutoff'] = np.sqrt(qe_system['ecutrho'])
    sirius_json['parameters']['gk_cutoff'] = np.sqrt(qe_system['ecutwfc'])
    # sirius_json['mixer']['beta'] = qe_electrons['mixing_beta']

    # store sirius.json
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
