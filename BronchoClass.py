from math import sqrt

class EasyBroncho(object):
    """
    Description of the class:
    This class is giving a basic model of the broncho thought as a parallel of a resistive
    element and a compliance.
    """
    eta = #1.81 * 10^-5 Pa*s https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    def __init__(self, generationNumber=False, length = 0.0, diameter = 0.0, resistance = 0.0, compliance = 0.0):
        if (generationNumber is False):
            self.length = length
            self.diameter = diameter
            self.resistance = resistance 
            self.compliance = compliance
        else:
            self.generationNumber = generationNumber
            self.length = length
            self.diameter = 1.8*0.5**(self.generationNumber/3.0)
            self.resistance = (8*eta*self.length)/(pi*(self.diameter/2)**4) #Poiseuille resistance


if __name__ == "__main__":
    trachea = EasyBroncho()
    print ("Trachea diameter is: %d" % trachea.diameter)