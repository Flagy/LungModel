from math import pi
from matplotlib import pyplot as plt
import numpy as np
import json
import logging
import math

logging.basicConfig(level=logging.DEBUG)

class EasyBroncho(object):
    """
    Description of the class:
    This class is giving a basic model of the broncho thought as a parallel of a resistive
    element and a compliance.
    All the physical quantities are expressed in the international system.
    """
    eta = 1.81e-5 #Pa*s https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    def __init__(self, generationNumber=False, paramsFromJson = False, length = None, diameter = None, resistance = None, compliance = None):
        if (generationNumber is False):
            self.length = length
            self.diameter = diameter
            self.resistance = resistance 
            self.compliance = compliance
        else:
            self.generationNumber = generationNumber
            if (paramsFromJson is False):
                self.length = length
                self.diameter = 0.018*0.5**(self.generationNumber/3.0) # 1.8 centimeters is 0.018 m
            else:
                # Reading from json file, this file is static for now. But it would be pointless (at least for now to make it dynamic)
                self.from_json()

            self.resistance = (8*EasyBroncho.eta*self.length)/(pi*(self.diameter/2)**4) # Poiseuille resistance   
            self.compliance = compliance

    def getTau(self):
        return self.resistance * self.compliance

    def from_json(self):
        with open("WeibelsModel_parameters.json") as json_file:
            d = json.load(json_file)['generation'][self.generationNumber]
        self.from_dict(d)

    def from_dict(self, d):
        for key, value in d.items():
            if type(value) is dict:
                value = EasyBroncho(value)
            self.__dict__[key] = value
    
    def from_reducedSkeleton_off(self, filename):
        ## Read the off file
        with open(filename, "r") as f:
            # Check the off header
            if f.readline().strip() != 'OFF':
                raise("Not a valid OFF file. --> Check header.")
            # Store info regarding number of vertices, faces and edges
            vertices, faces, edges = f.readline().strip().split(' ')
            ## Couple vertice index with coordinates
            vert = []
            edge = []
            for _ in range(int(vertices)):
                # everything is stored as a string, need conversion.
                vert.append(tuple([float(cord) for cord in f.readline().strip().split(' ')]))
            for _ in range(int(edges)):
                # Supposing the first number for now is always 2
                edge.append(tuple([int(num) for num in f.readline().strip().split(' ')[1:]]))
                
        ## Each edge corresponds to a resistance. Calculate it.
        # Supposing the first number for now is always 2
        # Calculate the distance between vert[0] and vert[1]
        terminal = []
        el1List = [e[0] for e in edge]
        el2List = [e[1] for e in edge]
        for el in el2List:
            terminal.append(el not in el1List)
        pass

        tree = extractTreeStructure(edge.copy(), terminal.copy())
        gen = max([len(el) for el in tree]) # create an empty matrix of the right dimensions.
        tree = createCompleteTree(tree, vert, gen)
        Rmatrix = getRmatrixFromTree(tree)
        return Rmatrix

def createCompleteTree(tree, vert, gen):
    """This function reconstruct the tree in a way to have all the terminals at the last generation
    The indexes won't be changed and"""
    for idx, path in enumerate(tree):
        if len(path) == gen: # this is the case of deeper path
            # overwriting the tree with the resistance value of each edge.
            tree[idx] = [getRfromEdge(vert, edge, currGen)[1] for currGen, edge in enumerate(path)] 
        else: # this is when the path is incomplete
            idxPath = 0
            line = []
            while(idxPath < len(path)):
                L, R = getRfromEdge(vert, path[idxPath], idxPath)
                line.append(R)
                idxPath = idxPath + 1
            while(idxPath < gen):
                previousLength = L
                L, R = getPhysioRfromGen(idxPath, previousLength)
                line.append(R)
                idxPath = idxPath + 1
            tree[idx] = line
    return tree

def getPhysioRfromGen(currGen, prevL):
    logger = logging.getLogger("getPhysioRfromGen")
    mu = 20e-9
    h = 0.5**(1/3)
    lengthEdge = prevL*h 
    radiusEdge = (0.018*0.5**(currGen/3.0))/2.0
    logger.debug(f"lengthEdge: {lengthEdge}")
    return (lengthEdge, 8*mu*lengthEdge/(np.pi*radiusEdge**4))
    
def getRfromEdge(vert, singleEdge, currGen):
    logger = logging.getLogger("getRfromEdge")
    mu = 20e-9
    idx1, idx2 = singleEdge
    lengthEdge = np.linalg.norm(np.asarray(vert[idx1]) - np.asarray(vert[idx2])) 
    radiusEdge = (0.018*0.5**(currGen/3.0))/2.0
    logger.debug(f"lengthEdge: {lengthEdge}")
    return (lengthEdge, 8*mu*lengthEdge/(np.pi*radiusEdge**4))

def extractTreeStructure(edge, terminal):
    """ 
    Example of input: np.array((0, 1) , (1, 2), (1, 3), (2, 4), (4, 7))
    Ex of output: [[(0, 1), (1, 2), (2,4), (4, 7), (7, 11)], [(0,1), (1,3), (3,5), (5, 8), (8, 12)]]
    """
    logger = logging.getLogger("extractTreeStructure")
    newTree = [] #root tree: (0, 1) toujours.
    branch = [(0, 1)]
    for _ in range(len(edge)):
        for j in range(len(edge)):
            if edge[j][0] == branch[-1][1]:
                branch.append(edge[j])
                if terminal[j]:
                    list2 = [el[1] for el in edge]
                    if branch[-1][0] in list2:
                        branch2 = list(branch)
                        last = branch2.pop()
                        newTree.append(branch)
                        list1 = [el[0] for el in edge]
                        list2 = [el[1] for el in edge]
                        el2toadd = last[1] + 1
                        try:
                            if terminal[list2.index(el2toadd)]:
                                branch2.append((last[0], el2toadd))
                                newTree.append(branch2)
                        except ValueError:
                            pass
                        for b1, b2 in zip(branch[::-1], branch2[::-1]):
                            list1 = [el[0] for el in edge]
                            if b1[1] not in list1:
                                terminal.pop(edge.index(b1))
                                edge.remove(b1)
                            if b1[1] != b2[1]:
                                if b2[1] not in list1:
                                    terminal.pop(edge.index(b2))
                                    edge.remove(b2)
                        branch = [(0, 1)]
                        break
                    else:
                        assert("Attention!")
    
    return newTree
            
def get_gen(ordered_tree):
    return max([len(el) for el in ordered_tree])

def p(t):
    t = t % 4
    if 0 <= t and t< 0.6:
        return -10 * np.sin(np.pi * t / 1.2)
    elif 0.6 <= t and t< 0.8:
        return -10
    elif 0.8 <= t and t < 1:
        return -10 * np.sin(np.pi * (t - 0.6) / 0.4)
    elif 1 <= t and t < 4:
        return 0
    else:
        raise ValueError('Check t in definition of p(t).')

if __name__ == "__main__":
    
    a = EasyBroncho()
    orderedTree = a.from_reducedSkeleton_off("skel-edge-reduced.off")
    print("Fine")
    """   
    R00 = EasyBroncho(0, paramsFromJson=True).resistance
    R10 = R11 = EasyBroncho(1, paramsFromJson=True).resistance / 2
    R20 = R21 = R22 = R23 = EasyBroncho(1, paramsFromJson=True).resistance / 4
    A = np.array([[R00+R10+R20, R00+R10, R00, R00], [R00+R10, R00+R10+R21, R00, R00], [R00, R00, R00+R11+R22, R00+R11], [R00, R00, R00+R11, R00+R11+R23]])
    Ainv = np.linalg.inv(A)
    tt = np.linspace(0, 10, 1000)
    P = []
    Q = []
    for t in tt:
        P.append(np.linspace(p(t), 0, 4))   
        Q.append(Ainv @ P[-1])
    Q = np.array(Q)
    print(P[0])
    
    for numberOfGen in range(len(Q[0])):
        Pressure = np.array([i[numberOfGen] for i in P])
        Flow = np.array([i[numberOfGen] for i in Q])
        # print(Flow)
        plt.plot(tt, Pressure, label="Pressure", color="red")
        plt.plot(tt, Flow*100, label="Flow", color="blue")
        plt.legend()
        plt.grid(linestyle='--')
        plt.title(f"Generation {numberOfGen+1}")
        plt.show()
    """


