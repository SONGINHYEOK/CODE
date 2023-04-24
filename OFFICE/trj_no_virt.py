'''
Read in a trajectory, remove all virtual sites, and write out a new
trajectory.

This script creates a PDB file and a trajectory in .XTC format. Use the
`--no-h` option to strip all hydrogen atoms from the output.

Since the order of hydrogen atoms may not be preserved when converting to PDB,
we recommend stripping all hydrogen atoms before proceeding with the analysis.

Proceed with caution.

We cannot guarantee the compatibility of Desmond trajectories with any module
outside of the Schrödinger Suite and we do not provide support for such
external tools.

Copyright Schrodinger, LLC. All rights reserved.
'''

from schrodinger.application.desmond.packages import cui
from schrodinger.application.desmond.packages import topo
from schrodinger.application.desmond.packages import traj


def parse_args():
    '''parse command line arguments'''

    # use `cui.OPTIONAL_TRJ` if a trajectory is optional
    parser = cui.CommandLine(description=__doc__,
                             spec=cui.REQUIRE_MSYS_CMS +
                                  cui.REQUIRE_TRJ +
                                  cui.SLICE_TRJ)

    parser.add_argument('-o', '--output', type=str,
                        default='md_no_virt',
                        help=('Basename for the output CMS and trajectory. '
                              '(Default: md_no_virt)'))

    parser.add_argument('--no-h', action='store_true',
                        help=('Strip hydrogen atoms from the trajectory. '
                              'Note: this will also strip hydrogen atoms '
                              'from water molecules.'))

    return parser.parse_args()


if __name__ == '__main__':

    # parse command line arguments
    args = parse_args()

    (msys_model, cms_model), _ = args.cms

    trj, _ = args.trj

    if args.no_h:
        asl = 'not a.e H'
    else:
        asl = 'all'

    # get the atom IDs (Maestro 1-based)
    atom_ids = cms_model.select_atom(asl)

    atom_ids_w_virtuals = topo.get_aids_with_virtuals(cms_model)
    cui.info(f'This cms/trajectory has {len(atom_ids_w_virtuals)} '
             'virtual sites.')

    # get the global IDs (Desmond's 0-based)
    global_ids = topo.aids2gids(cms_model, atom_ids,
                                include_pseudoatoms=False)

    # extract a structure object from the CMS structure
    new_st = cms_model.extract(atom_ids)

    # create a new trajectory without virtual sites
    new_trj = []
    for i, frame in enumerate(trj):
        cui.info(f'Processing frame {i+1}')
        new_frame = frame.reduce(global_ids, copy=True)
        new_trj.append(new_frame)

    # write new structure as PDB file (has no virtual sites by default)
    new_st.write(f'{args.output}.pdb', format='pdb')

    # write the new trajectory
    traj.write_traj(new_trj, f'{args.output}.xtc')
