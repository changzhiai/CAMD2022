from ase.io import read
from gpaw import GPAW, FermiDirac, PW

atoms = read('chosen_relax.traj')
calc = GPAW(mode=PW(500),
            kpts={'size': (8, 8, 8), 'gamma': True},
            xc='GLLBSC',
            occupations=FermiDirac(width=0.05))

atoms.calc = calc
energy = atoms.get_potential_energy()

# Note! An accurate discontinuity calculation requires a k-point grid that
# gives accurate HOMO/VBM and LUMO/CBM levels (in other words, the k-points of
# the valence band maximum and the conduction band minimum should be
# included in the used k-point grid).
homo, lumo = calc.get_homo_lumo()
response = calc.hamiltonian.xc.response
dxc_pot = response.calculate_discontinuity_potential(homo, lumo)
KS_gap, dxc = response.calculate_discontinuity(dxc_pot)
gap = KS_gap + dxc
print(f"The gap is {gap:.3f} with components: Kohn-Sham gap {KS_gap:.3f} and "
      f"discontinuity gap of {dxc:.3f}.")
