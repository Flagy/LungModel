from math import pi
import json

class EasyBroncho(object):
    """
    Description of the class:
    This class is giving a basic model of the broncho thought as a parallel of a resistive
    element and a compliance.
    All the physical quantities are expressed in the international system.
    """
    eta = 1.81e-5 #Pa*s https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    def __init__(self, generationNumber=False, paramsFromJson = False, length = 0.0, diameter = 0.0, resistance = 0.0, compliance = 0.0):
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


if __name__ == "__main__":
    bifurcation = EasyBroncho(1, paramsFromJson=True)
    print("Bifurcation generation is: %d" % bifurcation.generationNumber)
    print("Bifurcation diameter is: %f m" % bifurcation.diameter)
    print("Bifurcation length is: %f m" % bifurcation.length)
    print("Bifurcation resistance is: %f (Pa*s)/m^3" % bifurcation.resistance)