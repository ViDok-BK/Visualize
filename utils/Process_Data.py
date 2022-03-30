import os,re
import glob

from openbabel import openbabel

from utils.PDBReader import ComplexReader
from utils.SDFReader import SDFReader
from utils.mol_data import LigandGraph

current_dir = os.getcwd()


def convert_pdbqt_to_pdb():
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", "pdb")

    list_files = glob.glob(current_dir + '/data/receptor/*.pdbqt')

    for lf in list_files:
        name = lf.split('/')[-1].split('.')[0]
        path = lf.split('receptor')[0]
        try:
            obConversion.OpenInAndOutFiles(lf,path+'Data_Vidok_Infer/'+name+'.pdb')
            obConversion.Convert()
            obConversion.CloseOutFile()
        except:
            print(lf) 

def get_ranges():
    t = [(24,27) ,(41,49) ,(140, 145), (163, 169) ,(188, 192)]
    ranges = []
    for i in t:
        ranges += [j for j in range(i[0], i[1]+1)]
    return ranges

def get_pocket(f_name_in, f_name_out, ranges):
    lines= open(f_name_in, 'r').readlines()
    t=[[i for i in j.split(' ') if i!=''] for j in lines]
    new_lines = [line for i,line in enumerate(lines) 
            if t[i][0] not in ['ATOM', 'TER'] 
            or re.search('^\d+$', t[i][5]) and int(t[i][5]) in ranges
            or re.search('^\d+$', t[i][4]) and int(t[i][4]) in ranges]

    cnt=1
    for i,line in enumerate(new_lines):
        if t[0]!='ATOM':continue
        t[i][1]=str(cnt)
        cnt+=1
        line = ' '.join(t[i])
        new_lines[i]=line

    with open(f_name_out, 'w') as f:
        f.writelines(new_lines)

def process_pocket():
    ranges = get_ranges()
    list_files = glob.glob(current_dir + '/data/Data_ViDok_Infer/*.pdb')
    for lf in list_files:
        name = lf.split('/')[-1].split('.')[0]
        path = lf.split(name)[0]
        get_pocket(lf, path + '{}_pocket.pdb'.format(name), ranges)

def convert_sdf_after_docking():
    ligand_paths = glob.glob(current_dir+"/data/receptor/*.pdbqt")
    ligand_samples = [path.split(".")[0] for path in ligand_paths][:]
    # ligand_samples = ['/Users/duckhoan/Documents/VS_Code/ViDok/Auto_data/complex_ligand/ZINC000085862539_982']
    # i = 1
    for sample in ligand_samples:
        reader = ComplexReader()
        reader.parseComplex(sample + ".pdbqt")
        lig_graph,receptor = reader.returnLigandReceptor()

        lig_graph = LigandGraph(lig_graph)

        reader = SDFReader()
        start, end = reader.read(sample +".sdf")
        input_graph = LigandGraph(reader.molecule)

        result, mapping = lig_graph.update_connectivity(input_graph)

        # M  CHG  6   3   1   7   1   9   1  11   1  14   1  41   1


        l_end = end.split()

        if l_end[1]!='END':
            k = 3
            for _ in range(int(l_end[2])):
                l_end[k] = str(mapping[int(l_end[k])])
                k+=2
            end = "M  CHG  {}".format(l_end[2])
            for v in l_end[3:]:
                end += "{:>4s}".format(v)
            end += "\nM  END"

        result = "{}{}\n{}".format(start, result, end)

        name = sample.split('/')[-1]
        link = current_dir + '/data/Data_ViDok_Infer/{}'.format(name +".sdf")

        # print("{}: {}".format(i,link))
        # i+=1
        if not os.path.exists(link):
            with open(link, 'w') as f:
                f.write(result) 




def process_data():
    
    convert_pdbqt_to_pdb()
    process_pocket()
    convert_sdf_after_docking()
