import numpy as np
import warnings

#from bio_struct import Atom, Molecule, Residue, Protein, Ligand
from utils.structures import Atom, ReceptorAtom, Residue, Protein, Molecule
from utils.general_data import PERIODIC_TABLE, DONOR_ACCEPTOR_SYMBOL, ALL_ATOMS

def safeParseInt(string,begin,length):
    substr = string[begin:begin+length]
    substr = substr.strip()

    try:
        parsedNumber = int(substr)
    except:
        warnings.warn("Warning: cannot turn parsed String into integer")
        return 0
    if parsedNumber - float(parsedNumber) != 0:
        warnings.warn("Warning: Parsed target is not an int, but a float!")
    return parsedNumber

def safeParseFloat(string,begin,length):
    substr = string[begin:begin+length]
    substr = substr.strip()
    try:
        parsedNumber = float(substr)
    except:
        warnings.warn("Warning: cannot turn parsed String into a float, position =", begin)
        return 0

    parsedNumber = float(substr)
    return parsedNumber
    
def safeParseString(string,begin,length):
    if len(string) < begin - 1:
        warnings.warn("Warning: String index out of range")
        return ""

    len1 = length
    if len(string) < begin + length -1:
        len1 = len(string) - begin

    return string[begin:begin+len1]

class BaseReader:
    def __init__(self):
        self.macro = self.macro_class()
        self.macros_list = []
        
    def read(self,file_path):
        with open(file_path,"r") as o:
            #print("open")
            lines = o.readlines()
            for line in lines:
                cue = self.parse(line)
                if cue:
                    self.macros_list.append(self.macro)
                    self.macro = self.macro_class()

            self.finalize()
            
    def parse(self,line,recordName = None):
        #print("parsing")
        if not recordName:
            if len(line) < 6:
                #print("return")
                return
            recordName = safeParseString(line,0,6)
        if recordName == "ATOM  " or recordName == "HETATM":
            self.parseAtom(line)
            return 0
        elif recordName == "REMARK":
            return self.parseRemark(line)
        elif recordName == "TER   ":
            return 1
        else: 
            return self.parseSpecial(line,recordName)

    def parseUnknown(self,recordName):
        raise Exception("Unkonown Record Name: ",recordName)

    def parseRemark(self,line):
        return 0 
        
    def parseAtom(self,line):
        #sample line with index
        #          1         2         3         4         5         6         7
        #0     6     2   67   12       0       8       6   0         0         0     6         
        #ATOM   1768  CA  ARG A 188     -14.327  13.087  76.166  1.00 46.10     0.159 C
        atom = self.atom_class()

        atom.seq_num = safeParseInt(line,6,5)
        #here, node value is used to store the atomic symbols

        element = safeParseString(line,12,4)

        atom.parse_symbol(element)
        altLoc = safeParseString(line,16,1)
        if altLoc != " " and altLoc != "A":
            #print("return",altLoc,"end")
            return

        residue_symbol = safeParseString(line,17,4)
        res_seq_num = safeParseInt(line,22,5)
        #residue.chain_identifier = safeParseString(line,21,1)

        x = safeParseFloat(line,30,8)
        y = safeParseFloat(line,38,8)
        z = safeParseFloat(line,46,8)
        atom.coord = np.array([x,y,z])

        self.macro.append_atom(
            atom,
            residue_symbol = residue_symbol, 
            res_seq_num = res_seq_num)

        #print(atom.element)
        atom.atomic_number = int(PERIODIC_TABLE[atom.element.strip().title()][0])

        return 0

class ReceptorReader(BaseReader):
    def __init__(self):
        self.macro_class = Protein
        self.atom_class = ReceptorAtom
        super().__init__()

    def parseSpecial(self,line,recordName):
        if recordName == "SHEET":
            """
            line: SHEET    1  S1 2 THR A   1  CYS A   4  0
            startChain = A, startResi = 1, endResi = 4
            """
            raise Exception("Not yet implemented: ",recordName) 
        elif recordName =="HELIX":
            raise Exception("Not yet implemented: ",recordName) 
        else:
            self.parseUnknown(recordName)

    def finalize(self):
        pass

class LigandReader(BaseReader):
    def __init__(self):
        self.macro_class = Molecule
        self.atom_class = Atom
        super().__init__()

    def parseSpecial(self,line,recordName):
        if recordName == "ENDROO":
            return 0
        elif recordName == "BRANCH":
            return 0
        elif recordName == "ENDBRA":
            return 0
        elif recordName == "TORSDO":
            return 1
        else:
            self.parseUnknown(recordName)

    def finalize(self):
        pass


class ComplexReader(BaseReader):
    def __init__(self):
        self.macro_class = Molecule
        super().__init__()
        self.receptor_reader = ReceptorReader()
        self.ligand_reader = LigandReader()
        self.reader = None

    def parseRemark(self,line):
        if safeParseString(line,7,4) == "VINA":
            self.reader.finalize()
            self.reader = self.ligand_reader
        return 0

    def parse(self,line):
        if len(line) < 6:
            #print("return")
            return
        recordName = safeParseString(line,0,6)
        if recordName == "REMARK":
            return self.parseRemark(line)
        
        return self.reader.parse(line,recordName)

    def parseComplex(self,complex_path):
        self.reader = self.receptor_reader
        print(complex_path)
        with open(complex_path, "r") as o:
            lines = o.readlines() 
            # print(lines)
            for line in lines:
                cue = self.parse(line)
                if cue:
                    self.reader.macros_list.append(self.reader.macro)
                    self.reader.macro = self.reader.macro_class()
    
    def returnLigandReceptor(self):
        ligand = self.ligand_reader.macros_list
        assert len(ligand) == 1
        ligand = ligand[0]

        receptor = self.receptor_reader.macros_list
        assert len(receptor) == 1
        receptor = receptor[0]

        

        return ligand, receptor
