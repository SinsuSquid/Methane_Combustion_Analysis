import sys
import numpy as np

RCUT = 2
DUMP_FREQ = 1000

atominfo = {}
coord = {}
pbc = {}

def checkFile():
    if (len(sys.argv) < 2):
        print("USAGE : python3 simple_cluster.py XX.lammpstrj")
        exit(1)

def readTimestep():
    with open(sys.argv[1],'r') as fp:
        while(True):
            line = fp.readline().strip()
            if (line == "ITEM: TIMESTEP"): 
                timestep = int(fp.readline())
                if (timestep % 100000 == 0): print(f"Now reading STEP : {timestep} ...")
            elif (line == "ITEM: NUMBER OF ATOMS"):
                numAtoms = int(fp.readline())
            elif (line == "ITEM: BOX BOUNDS pp pp pp"):
                checkPBC(fp, timestep)
            elif (line == "ITEM: ATOMS id type x y z"):
                getTraj(fp, numAtoms, timestep)
            elif (line == ''): break;

def checkPBC(fp, timestep):
    xmin, xmax = fp.readline().split(' ')
    ymin, ymax = fp.readline().split(' ')
    zmin, zmax = fp.readline().split(' ')
    data = np.array([[xmin, xmax],
                     [ymin, ymax],
                     [zmin, zmax]], dtype = float)
    pbc[timestep] = data[:,1] - data[:,0]

def getTraj(fp, numAtoms, timestep):
    idType = []
    xyz = []
    for i in range(numAtoms):
        tmp = fp.readline().split()
        idType.append([tmp[0], tmp[1]])
        xyz.append([tmp[2], tmp[3], tmp[4]])
    idType = np.array(idType, dtype = int)
    xyz = np.array(xyz, dtype = float)
    atominfo[timestep] = idType
    coord[timestep] = xyz

def checkCluster():
    global atominfo; global coord; global pbc; global RCUT; global DUMP_FREQ

    for i in range(len(coord)):
        if (i % 10 == 0): print(f"Now cheking clusters for STEP : {i} ...")
        clusterList = []
        coord_thisStep = coord[i * DUMP_FREQ]
        atominfo_thisStep = atominfo[i * DUMP_FREQ]

        for j in range(len(coord_thisStep)):
            clusterSet = set()
            clusterSet.add(atominfo_thisStep[j][0])
            ref = coord_thisStep[j]
            for k in range(j + 1, len(coord_thisStep)):
                cmp = coord_thisStep[k]
                length = getLength(ref, cmp, pbc[i * DUMP_FREQ])
                if (length < RCUT):
                    clusterSet.add(atominfo_thisStep[k][0])
            print(clusterSet)


def getLength(ref, cmp, pbc):
    length = 0.0
    for i in range(3): # for x,y,z
        perAxis = np.power(np.power(ref[i] - cmp[i], 2),0.5)
        if ((perAxis - 0.5 * pbc[i]) > 0): 
            length += np.power((perAxis - pbc[i]),2)
        else : 
            length += np.power(perAxis,2)
    return np.power(length, 0.5)

if __name__ == "__main__":
    checkFile()
    readTimestep()

    checkCluster()
