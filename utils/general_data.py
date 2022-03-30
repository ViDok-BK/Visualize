ALPHABET = [chr(i) for i in range(65,91)]
REMOTENESS = ["A","B","G","D","E","Z","H"]

PERIODIC_TABLE = {
    "H":[1,1],
    "C":[6,12], "N":[7,14],
    "P":[15,31],
    "S":[16,32], "O":[8,16],
    "F":[9,19], "Cl":[17,35.45], 
    "Br":[35,80], "I":[53,127]
}

CHALCOGEN = ["O","S","P"]
HALOGEN = ["F","Cl","Br","I"]
PNICTOGEN = ["N"]
CARBON = ["C"]
HYDROGEN = ["H"]

ALL_ATOMS = CHALCOGEN + HALOGEN + PNICTOGEN + CARBON + HYDROGEN

DONOR_ACCEPTOR_SYMBOL = ["D","A"]

def MAXIMUM_BOND_LENGTH(element1,element2):
    element1 = element1.strip()
    element2 = element2.strip()
    if element1 == "Cl" or element2 == "Cl":
        return 1.9
    else:
        return 1.78

class AtomClass:
    def __init__(
            self, include_atoms = ALL_ATOMS, 
            exclude_atoms = None, rules = None
            ):
        """
        Rules for atom to fulfill if it is in the class
        """
        self.inclusion = include_atoms
        self.exclusion = [] if exclude_atoms is None else exclude_atoms 
        self.rules = rules 

    def isIn(self,atom):
        element = atom.element.title()
        if element in self.exclusion:
            return False
        if element not in self.inclusion:
            return False
        if self.rules:
            for rule in self.rules:
                if not rule(atom):
                    return False
        return True

#should include rule for further filtering out
H_Acceptor= AtomClass(
            include_atoms = ALL_ATOMS, 
            exclude_atoms = HALOGEN[1:] + CARBON + HYDROGEN
            )
H_Donor_receptor = AtomClass(
            include_atoms = ["H"] )
H_Donor_ligand = AtomClass(
            exclude_atoms = HALOGEN + ["C"], 
            rules = [
                lambda x: not x.element == "O" or x.total_connectivity() < 2,
                lambda x: not x.element == "S" or x.total_connectivity() < 2,
                lambda x: not x.element == "N" or x.total_connectivity() < 3
                ])

#negatively polar atoms, since positively polar atoms are nothing but hydrogen
polar_atoms = AtomClass(
            include_atoms = HALOGEN[1:])

