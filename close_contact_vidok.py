
import numpy as np
import matplotlib.pyplot as plt
import glob
import csv
import re
import os

from utils.PDBReader import ComplexReader
from mapping_pocket import get_mapping


# current_dir = os.getcwd()

PATH_INTERACT = '../interactions/'
PATH_DATA = '../Vidok_Crawler/complex_ligand/'


def parse_type_interaction(list_interactions):
    interactions = list_interactions.split(';')
    result = []
    for interaction in interactions:
        res = re.findall("[a-zA-Z]+",interaction)
        if len(res)!=0:
            result.append(res[0])

    return result

class BindingSite:
    def __init__(self,close_contact, receptor):
        self.freq = 1
        self.H_donor = 0
        self.H_acceptor = 0
        self.hydrophobic = 0

        self.update_info(close_contact, receptor)

    def update_info(self,close_contact,receptor):
        for res in receptor:
            for atom in res.atoms:
                if atom.seq_num == close_contact[1]:
                    self.atom = atom.element
                    self.atom_seq_num = atom.seq_num
                    self.atom_coord = atom.coord
        
        self.get_contact(close_contact)
        

    def get_contact(self,close_contact):
        if len(close_contact[2])!=0 and len(close_contact[3])!=0:
            if 'Hydropobe' in close_contact[4]:
                self.hydrophobic+=1
            else:
                if 'Acceptor' in close_contact[3] and 'Donor' in close_contact[2]:
                    self.H_acceptor+=1
                elif 'Acceptor' in close_contact[2] and 'Donor' in close_contact[3]:
                    self.H_donor+=1

    
    def merge_close_contact(self,close_contact):
        assert close_contact[1] == self.atom_seq_num
        
        self.freq += 1

        self.get_contact(close_contact)

    def __repr__(self):
        return str(self.atom_seq_num)

class BindingCavity:

    def __init__(self):
        self.binding_sites = {}

    def add_close_contact(self,close_contact, receptor):
        atom_no = close_contact[1]

        if atom_no in self.binding_sites.keys():
            self.binding_sites[atom_no].merge_close_contact(close_contact)
        else:
            bind_site = BindingSite(close_contact, receptor)
            self.binding_sites.update({atom_no:bind_site})

    def get_binding_sites(self):
        Hdonor_sites = []
        Haccep_sites = []
        Hphobic_sites = []
        other_sites = []
        for bs in self.binding_sites.values():
            num_binding = [bs.H_donor, bs.H_acceptor, bs.hydrophobic]
            if np.sum(num_binding)!=0:
                i_max = np.argmax([bs.H_donor, bs.H_acceptor, bs.hydrophobic])
                if i_max==0 :
                    Hdonor_sites.append(bs)
                elif i_max==1:
                    Haccep_sites.append(bs)
                elif i_max==2:
                    Hphobic_sites.append(bs)
                else:
                    other_sites.append(bs)
            else:
                other_sites.append(bs)
                
        return Hdonor_sites, Haccep_sites, Hphobic_sites, other_sites

    def print_binding(self):
        for k,v in self.binding_sites.items():
            print('Atom index: {}---{} {} {}'.format(k,v.H_donor, v.H_acceptor, v.hydrophobic))



         
def visualize():

    interaction_paths = glob.glob(PATH_INTERACT+"*.csv")
    names = sorted([n.split("/")[-1].split('.')[0]for n in interaction_paths])

    cavity = BindingCavity()
    map_pocket = get_mapping()
    
    for name in names:

        #read ligand and receptor from complex.pdbqt
        reader = ComplexReader()
        reader.parseComplex(PATH_DATA+"{}.pdbqt".format(name))
        lig_graph,receptor = reader.returnLigandReceptor()

        #read interaction 
        with open(PATH_INTERACT+'{}.csv'.format(name), mode='r', encoding='UTF8') as f:
            csvreader = csv.reader(f)
            header = next(csvreader)
            for row in csvreader:
                contact = []
                for i in range(len(row)):
                    if i==0:
                        contact.append(int(row[i]))
                    elif i==1:
                        contact.append(map_pocket[int(row[i])])
                    else:
                        contact.append(parse_type_interaction(row[i]))
                #Add and update cavity
                cavity.add_close_contact(contact,receptor)
    
    Hdonor_sites, Haccep_sites, Hphobic_bind_sites, other_sites = cavity.get_binding_sites()

    print("total number of ligands: {:.0f}".format(len(names)))
    print("Number of hydrogen donor sites: {:.0f}".format(len(Hdonor_sites)))
    print("Number of hydrogen acceptor sites: {:.0f}".format(len(Haccep_sites)))
    print("Number of hydrophobic sites: {:.0f}".format(len(Hphobic_bind_sites)))

    fig = plt.figure()
    ax = plt.axes(projection = "3d")
    size = 300

    #visualize sites that has hydrogen bond DONOR
    bind_sites = sorted(Hdonor_sites,key = lambda x: x.atom_seq_num)
    coords = [site.atom_coord for site in bind_sites]
    coords = np.array(coords)

    ax.scatter(
        coords[:,0],
        coords[:,1],
        coords[:,2],
        color = "blue",
        s = size,
        )

    #ax.plot(coords[:,0],coords[:,1],coords[:,2],color = "red")

    #visualize sites that has hydrogen bond ACCEPTOR
    bind_sites = sorted(Haccep_sites,key = lambda x: x.atom_seq_num)
    coords = [site.atom_coord for site in bind_sites]
    coords = np.array(coords)

    ax.scatter(
        coords[:,0],
        coords[:,1],
        coords[:,2],
        color = "red",
        s = size
        )

    #ax.plot(coords[:,0],coords[:,1],coords[:,2],color = "blue")

    #visualize sites that accomodate hydrophobic interaction
    bind_sites = sorted(Hphobic_bind_sites,key = lambda x: x.atom_seq_num)
    coords = [site.atom_coord for site in bind_sites]
    coords = np.array(coords)

    ax.scatter(
        coords[:,0],
        coords[:,1],
        coords[:,2], 
        color = "gray",
        s = size,
        )

    #visualize sites that accomodate unknown interaction
    bind_sites = sorted(other_sites,key = lambda x: x.atom_seq_num)
    coords = [site.atom_coord for site in bind_sites]
    coords = np.array(coords)

    ax.scatter(
        coords[:,0],
        coords[:,1],
        coords[:,2], 
        color = "green",
        s = size,
        )

    #ax.plot(coords[:,0],coords[:,1],coords[:,2], color = "green")
    plt.show()







