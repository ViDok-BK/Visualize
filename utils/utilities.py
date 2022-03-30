import numpy as np

def sqrt_dist(x,x0):
    return np.sum((x - x0)**2)

def dist(x,x0):
    return np.sqrt(sqrt_dist(x,x0))

class Sphere:
    """
    Create a sphere object to check whether or not a given 
    point inside that sphere
    Args:
        + center (numpy array): center of the sphere
        + radius (float): radius of the sphere
    """
    def __init__(self,center,sqrt_radius):
        if isinstance(center,list):
            center = np.array(center)
        self.center = center
        self.sqrt_radius = sqrt_radius

    def isIn(self,coord):
        """
        list of np.array -> bool
        return True if the point of given coordinate is inside the 
        sphere and False other wise
        """
        if (sqrt_dist(coord,self.center) < self.sqrt_radius):
            return True
        return False
        
def spherical_molecule(mol):
    center = np.zeros(3)
    for atom in mol.atoms:
        center += atom.coord
    center /= len(mol.atoms)
    
    sqrt_radius = 0
    for atom in mol.atoms:
        sqrt_radius = max(sqrt_radius,sqrt_dist(center,atom.coord))

    return Sphere(center = center, sqrt_radius = sqrt_radius)


