###############################################################################
# Encoding utf-8                                                              #
# F. Madeira and L. Krippahl, 2012                                            #
# This code is part of Pycoevol distribution.                                 #
# This work is public domain.                                                 #
###############################################################################

"""Adapted from Surface Area (ASA) - (C) Bosco Ho 
http://boscoh.com/protein/calculating-the-solvent-accessible-surface-area-asa

Calculates the Solvente Accessible Surface Area (SASA) using the classic
'rolling ball' algorithm - A. Shrake & J. A. Rupley.
Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin.
J Mol Biol. 79 (1973) 351-371.
"""

from src.UTILS import radii
from math import pi, sqrt, cos, sin

SMALL = 1E-6
two_char_elements = [el for el, r in radii.items() if len(el) == 2]

def SASA(input, output):
    mol = Molecule(input)
    atoms = mol.atoms()
    add_radii(atoms)
    
    n_sphere = 960
    asas = calculateSASA(atoms, 1.4, n_sphere)
    
    for asa, atom in zip(asas, atoms):
        atom.bfactor = asa
    mol.write_pdb(output)
    return

def generateSpherePoints(n):
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.
    """
    points = []
    inc = pi * (3 - sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = sqrt(1 - y * y)
        phi = k * inc
        points.append([cos(phi) * r, y, sin(phi) * r])
    return points


def findNeighborIndices(atoms, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + probe + probe
    indices = range(k)
    indices.extend(range(k + 1, len(atoms)))
    for i in indices:
        atom_i = atoms[i]
        dist = pos_distance(atom_k.pos, atom_i.pos)
        if dist < radius + atom_i.radius:
            neighbor_indices.append(i)
    return neighbor_indices


def calculateSASA(atoms, probe, n_sphere_point=960):
    """
    Returns list of accessible surface areas of the atoms, using the probe
    and atom radius to define the surface.
    """
    sphere_points = generateSpherePoints(n_sphere_point)

    const = 4.0 * pi / len(sphere_points)
    test_point = Vector3d()
    areas = []
    for i, atom_i in enumerate(atoms):
        neighbor_indices = findNeighborIndices(atoms, probe, i)
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0
        radius = probe + atom_i.radius

        n_accessible_point = 0
        for point in sphere_points:
            is_accessible = True

            test_point.x = point[0] * radius + atom_i.pos.x
            test_point.y = point[1] * radius + atom_i.pos.y
            test_point.z = point[2] * radius + atom_i.pos.z

            cycled_indices = range(j_closest_neighbor, n_neighbor)
            cycled_indices.extend(range(j_closest_neighbor))

            for j in cycled_indices:
                atom_j = atoms[neighbor_indices[j]]
                r = atom_j.radius + probe
                diff_sq = pos_distance_sq(atom_j.pos, test_point)
                if diff_sq < r * r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1

        area = const * n_accessible_point * radius * radius 
        areas.append(area)
    return areas


def add_radii(atoms):
    for atom in atoms:
        if atom.element in radii:
            atom.radius = radii[atom.element]
        else:
            atom.radius = radii['.']

def pos_distance_sq(p1, p2):
    x = p1.x - p2.x
    y = p1.y - p2.y
    z = p1.z - p2.z
    return x * x + y * y + z * z;

def pos_distance(p1, p2):
    return sqrt(pos_distance_sq(p2, p1))

class Molecule:
    def __init__(self, pdb=""):
        self.id = ''
        self._atoms = []
        if pdb:
            self.read_pdb(pdb)

    def n_atom(self):
        return len(self._atoms)

    def atoms(self):
        return self._atoms

    def atom(self, i):
        return self._atoms[i]
        
    def clear(self):
        for atom in self._atoms:
            del atom
        del self._atoms[:]

    def transform(self, matrix):
        for atom in self._atoms:
            atom.pos.transform(matrix)

    def insert_atom(self, atom):
        self._atoms.append(atom)
        
    def erase_atom(self, atom_type):
        for atom in self._atoms:
            if atom.type == atom_type:
                self._atoms.remove(atom)
                del atom
                return

    def read_pdb(self, fname):
        self.clear()
        for line in open(fname, 'r').readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = AtomFromPdbLine(line);
                if len(self._atoms) == 1:
                    self.id = atom.chain_id
                self.insert_atom(atom)
            if line.startswith("ENDMDL"):
                return

    def write_pdb(self, pdb):
        f = open(pdb, 'w')
        n_atom = 0
        for atom in sorted(self._atoms, cmp=cmp_atom):
            n_atom += 1
            atom.num = n_atom
            f.write(atom.pdb_str() + '\n')
        f.close()

def AtomFromPdbLine(line):
    """Returns an Atom object from an atom line in a pdb file."""
    atom = Atom()
    if line.startswith('HETATM'):
        atom.is_hetatm = True
    else:
        atom.is_hetatm = False
    atom.num = int(line[6:11])
    atom.type = line[12:16].strip(" ")
    element = ''
    for c in line[12:15]:
        if not c.isdigit() and c != " ":
            element += c
    if element[:2] in two_char_elements:
        atom.element = element[:2]
    else:
        atom.element = element[0]
    atom.res_type = line[17:20]
    atom.chain_id = line[21]
    atom.res_num = int(line[22:26])
    atom.res_insert = line[26]
    if atom.res_insert == " ":
        atom.res_insert = ""
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom.pos.set(x, y, z)
    try:
        atom.occupancy = float(line[54:60])
    except:
        atom.occupancy = 100.0
    try:
        atom.bfactor = float(line[60:66])
    except:
        atom.bfactor = 0.0
    return atom
    
    
def cmp_atom(a1, a2):
    if a1.num < a2.num:
        return -1
    else:
        return 0

def pad_atom_type(in_atom_type):
    atom_type = in_atom_type
    if len(atom_type) == 1:
        atom_type = " %s    " % atom_type
    elif len(atom_type) == 2:
        atom_type = " %s " % atom_type
    elif len(atom_type) == 3:
        if atom_type[0].isdigit():
            atom_type = "%s " % atom_type
        else:
            atom_type = " %s" % atom_type
    return atom_type

class Atom:
    def __init__(self):
        self.is_hetatm = False
        self.pos = Vector3d()
        self.vel = Vector3d()
        self.mass = 0.0
        self.type = ""
        self.element = ""
        self.chain_id = " "
        self.res_type = ""
        self.res_num = ""
        self.res_insert = ""
        self.bfactor = 0.0
        self.occupancy = 0.0
        self.num = 0
    
    def pdb_str(self):
        return str(self.chain_id) + "\t" + str(self.res_type) + "\t" + \
            str(self.res_num) + "\t" + str(self.bfactor)
                             
    def __str__(self):
        return "%s%s-%s (% .1f % .1f % .1f)" \
                        % (self.res_type, self.res_num,
                                self.type, self.pos.x,
                                self.pos.y, self.pos.z)
        
class Vector3d:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, rhs):
        return Vector3d(rhs.x + self.x, rhs.y + self.y, rhs.z + self.z)

    def __sub__(self, rhs):
        return Vector3d(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)

    def __neg__(self):
        return Vector3d(-self.x, -self.y, -self.z)

    def __pos__(self):
        return Vector3d(self.x, self.y, self.z)

    def __eq__(self, rhs):
        return (is_near_zero(self.x - rhs.x) and \
                        is_near_zero(self.y - rhs.y) and \
                        is_near_zero(self.z - rhs.z))

    def __str__(self):
        return "(% .2f, % .2f, % .2f)" % (self.x, self.y, self.z)

    def __repr__(self):
        return "Vector3d(%f, %f, %f)" % (self.x, self.y, self.z)

    def set(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def copy(self):
        return Vector3d(self.x, self.y, self.z)

    def length_sq(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    def length(self):
        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def scale(self, scale):
        self.x *= scale
        self.y *= scale
        self.z *= scale

    def normalize(self):
        self.scale(1.0 / self.length())

    def scaled_vec(self, scale):
        v = self.copy()
        v.scale(scale)
        return v

    def normal_vec(self):
        return self.scaled_vec(1.0 / self.length())

    def parallel_vec(self, axis):
        axis_len = axis.length()
        if is_near_zero(axis_len):
            result = self
        else:
            result = axis.scaled_vec(dot(self, axis) 
                             / axis.length() / axis.length())
        return result

    def perpendicular_vec(self, axis):
        return self - self.parallel_vec(axis)

    def transform(self, matrix):
        x = matrix.elem00 * self.x + \
                matrix.elem10 * self.y + \
                matrix.elem20 * self.z + \
                matrix.elem30
        y = matrix.elem01 * self.x + \
                matrix.elem11 * self.y + \
                matrix.elem21 * self.z + \
                matrix.elem31
        z = matrix.elem02 * self.x + \
                matrix.elem12 * self.y + \
                matrix.elem22 * self.z + \
                matrix.elem32
        self.x, self.y, self.z = x, y, z

def is_near_zero(a):
    return a < SMALL
  
def dot(a, b):
    return a.x * b.x + a.y * b.y + a.z * b.z

