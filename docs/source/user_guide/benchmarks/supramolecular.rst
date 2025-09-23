==============
Supramolecular
==============

S30L
======

Summary
-------

Performance in predicting supramolecular binding energies for the S30L set of
30 large host–guest complexes. Systems contain up to ~200 atoms and feature a
wide range of interaction motifs, including hydrogen and halogen bonding,
π–π stacking, CH–π contacts, nonpolar dispersion, and cation–dipolar
interactions. Net charges span −1 to +4, where 8 of the 30 complexes are charged. These ΔE_emp values serve as the
benchmark dataset for assessing quantum-chemical methods on large noncovalent
complexes.

Metrics
-------

Total binding energy error

For each complex, the binding energy is calculated by taking the difference in energy
between the host-guest complex and the sum of the individual host and guest energies. This is
compared to the reference binding energy, calculated in the same way.

Charged binding energy error

For each charged complex, the binding energy is calculated by taking the difference in energy
between the host-guest complex and the sum of the individual host and guest energies. This is
compared to the reference binding energy, calculated in the same way.

Neutral binding energy error

For each neutral complex, the binding energy is calculated by taking the difference in energy
between the host-guest complex and the sum of the individual host and guest energies. This is
compared to the reference binding energy, calculated in the same way.

Computational cost
------------------

low: tests are likely to take minutes to run on CPU.

Data availability
-----------------

Input structures:

* R. Sure and S. Grimme, ‘S30L - Comprehensive Benchmark of Association (Free) Energies of Realistic Host–Guest Complexes’, J. Chem. Theory Comput., vol. 11, no. 8, pp. 3785–3801, Aug. 2015, doi: 10.1021/acs.jctc.5b00296.
* Stuctures download found in SI
    * TPSS-D3/def2-TZVP geometries of all complexes as Cartesian coordinates.

Reference data:

* The Supporting Information also provides the benchmark binding energies.
    * These are empirical gas-phase reference values (ΔE_emp), obtained by back-correcting the experimental association free energies with theoretical estimates of vibrational and solvation contributions.
