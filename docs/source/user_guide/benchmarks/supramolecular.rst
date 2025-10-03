==============
Supramolecular
==============

PLA15
======

Summary
-------

Performance in predicting protein–ligand active-site interaction energies for the
PLA15 set of 15 complexes. Systems range from 259 to 584 atoms and contain complete
active sites. Ligands contain 37–95 atoms with net charges of −1, 0, or +1 and all contain
aromatic heterocycles. Five ligands contain either divalent of tetrahedral sulfur atoms, and
four and three of them contain F and Cl atoms, respectively.

Metrics
-------

Total MAE

For each complex, the interaction energy is calculated by taking the difference in energy between the protein-ligand complex and the sum of the individual protein and ligand energies. The MAE is computed by comparing predicted interaction energies to reference interaction energies across all 15 systems.

Pearson's r²

The squared Pearson correlation coefficient between predicted and reference interaction energies, measuring the proportion of variance in the reference values explained by the model predictions.

Ion-Ion MAE

For each complex where both protein and ligand fragments have non-zero charges, the interaction energy error is calculated. This metric reports the MAE for these ion-ion interaction systems.

Ion-Neutral MAE

For each complex where one fragment (protein or ligand) has a non-zero charge and the other is neutral, the interaction energy error is calculated. This metric reports the MAE for these ion-neutral interaction systems.


Computational cost
------------------

low: tests are likely to take minutes to run on CPU.

Data availability
-----------------

Input structures:

* K. Kříž and J. Řezáč, ‘protein ligand - Benchmarking of Semiempirical Quantum-Mechanical Methods on Systems Relevant to Computer-Aided Drug Design’, J. Chem. Inf. Model., vol. 60, no. 3, pp. 1453–1460, Mar. 2020, doi: 10.1021/acs.jcim.9b01171.
* Structures download found in SI

Reference data:

* The Supporting Information also provides the interaction energies.
    * The benchmark interaction energies are based on a combination of explicitly correlated MP2-F12 calculations and a DLPNO-CCSD(T) correction
