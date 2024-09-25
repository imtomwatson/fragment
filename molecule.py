import copy

class Molecule:
    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds

    def __repr__(self):
        return f"Molecule(atoms={self.atoms}, bonds={self.bonds})"


    def grow_molecule(self, atom_idx, bonds, bonded_atom):
        bond = {atom_idx: bonds[atom_idx]}
        for atom in bond[atom_idx]:
            if atom == bonded_atom:
                continue
            bond.update(self.grow_molecule(atom, bonds, atom_idx))
        return bond
    
    
    def get_fragments(self, atom_idx1, atom_idx2):
        atoms = self.atoms
        bonds = self.bonds
    
        if atom_idx2 not in bonds[atom_idx1] or atom_idx1 not in bonds[atom_idx2]:
            raise ValueError("Atoms are not bonded")
    
        # Break the bond in the bonds dict
        bonds[atom_idx1].remove(atom_idx2)
        bonds[atom_idx2].remove(atom_idx1)
    
        molecules = []
        # Grow the molecules from the atoms bonded to the ... bonded atoms
        for idx in [atom_idx1, atom_idx2]:
            bond_dict = {idx: copy.deepcopy(bonds[idx])}
            for atom in bond_dict[idx]:
                bond_dict.update(self.grow_molecule(atom, bonds, idx))

            for atom, bond_list in bond_dict.items():
                _atoms = {i: atoms[i] for i in sorted(list(set([atom] + bond_list)))}

            molecules.append(Molecule(_atoms, bond_dict))
        return molecules


def test_ccn():
    test_data = {
        (0, 1): [Molecule({0: "C"}, {0: []}), Molecule({1: "C", 2: "N"}, {1: [2], 2: [1]})],
        (1, 2): [Molecule({0: "C", 1: "C"}, {0: [1], 1: [0]}), Molecule({2: "N"}, {2: []})]
    }

    for (atom_idx1, atom_idx2), (m1, m2) in test_data.items():
        atoms = {0: "C", 1: "C", 2: "N"}
        bonds = {0: [1], 1: [0,2], 2: [1]}
        molecule = Molecule(atoms, bonds)
        print(f"Breaking bond between {atom_idx1} and {atom_idx2} in {molecule}")

        r1, r2 = molecule.get_fragments(atom_idx1, atom_idx2)
        for a, b in [(r1, m1), (r2, m2)]:
            print(a)
            print(b)
            assert a.atoms == b.atoms, f"Atoms not equal! {a.atoms} != {b.atoms}"
            assert a.bonds == b.bonds, f"Bonds not equal! {a.bonds} != {b.bonds}"

if __name__ == "__main__":
    test_ccn()