from math import sqrt

class EasyBroncho(object):
    """
    Description of the class:
    This class is giving a basic model of the broncho thought as a parallel of a resistive
    element and a compliance.
    """
    def __init__(self):
        self.length = 0.0
        self.diameter = 0.0
        self.R = 0.0 
        self.C = 0.0


if __name__ == "__main__":
    trachea = EasyBroncho()
    print ("Trachea diameter is: %d" % trachea.diameter)