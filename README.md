**PCMN - Protein Contact Map Network**

Leveraging the "Contact Map Explorer" algorithm by Swenson and Roet, this script compares the contact maps of two related PDB structures. These structures could be, for instance, different conformers extracted from a molecular dynamics trajectory or obtained from crystallography. The script not only calculates the contact map for each input structure but also pinpoints the contacts that are unique to each. Furthermore, by providing a list of specific residues, users can generate a network visualization highlighting their close-range interactions. The default distance threshold for defining a contact is set to 0.35 nanometers, but this value can be modified by the user. Essentially, this script acts as a supplementary utility for the primary Contact Map Explorer tool.

**Installation of dependencies**
'''
pip install matplotlib
pip install argparse
pip install mdtraj
pip install cython
pip install numpy
pip install contact_map
pip install networkx
'''
