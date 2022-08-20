
# Full script
# from ase...

# Teacher
import numpy as np
from ase.lattice.hexagonal import Graphite
from ase.calculators.dftd3 import DFTD3
from ase.constraints import StrainFilter
from ase.optimize.bfgs import BFGS
from ase.io import Trajectory

ccdist = 1.41
layerdist = 3.21
a = ccdist * np.sqrt(3)
c = 2 * layerdist
gra = Graphite('C', latticeconstant={'a': a, 'c': c})

from gpaw import GPAW, PW

for xc in ['LDA', 'PBE', 'DFTD3']:
    if xc == 'DFTD3':
        dft = GPAW(mode=PW(500), kpts=(10, 10, 6), xc='PBE',
                   txt=f'graphite-{xc}.log')
        calc = DFTD3(dft=dft, xc='PBE')
    else:
        calc = GPAW(mode=PW(500), kpts=(10, 10, 6), xc=xc,
                    txt=f'graphite-{xc}.log')

    gra.calc = calc  # Connect system and calculator

    sf = StrainFilter(gra, mask=[1, 1, 1, 0, 0, 0])
    opt = BFGS(sf, trajectory=f'graphite-{xc}.traj')
    # traj = Trajectory(f'graphite-{xc}.traj', 'w', gra)
    # opt.attach(traj)
    opt.run(fmax=0.01)
