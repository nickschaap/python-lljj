import sys
import math

########## CLASS DEFINITIONS ##########
# Class object that represents one event
class Event:
    def __init__(self, header):
        self.eventNum = int(header[1])
        self.particles = []

    def addFinalStateParticle(self, particle) :
        self.particles.append(particle)

    def selectionCut(eventchar) :
        return

# Class object for a final state particle
class Particle:
    def __init__(self, event, data):
        self.event = event
        self.num = int(data[0])
        self.typ = int(data[1])
        self.eta = float(data[2])
        self.phi = float(data[3])
        self.pt = float(data[4])
        self.jmas = float(data[5])
        self.ntrk = float(data[6])
        self.btag = float(data[7])
        self.hadem = float(data[8])
        self.dum1 = float(data[9])
        self.dum2 = float(data[10])
        self.px = float(self.pt * math.cos(self.phi))
        self.py = float(self.pt * math.sin(self.phi))
        self.theta = float(2 * math.atan(math.exp(-self.eta)))
        self.pz = float(self.pt / math.tan(self.theta))
        self.p = float(math.sqrt(math.pow(self.px, 2) + math.pow(self.py, 2) + math.pow(self.pz, 2)))
    
    def ptVec(self) :
        return (self.px, self.py)

    def fourVector(self) :
        return (self.jmas, self.px, self.py, self.pz)

def fourDotProduct(particle1, particle2) :
    fourVec1 = particle1.fourVector()
    fourVec2 = particle2.fourVector()
    return fourVec1[0] * fourVec2[0] - (fourVec1[1]*fourVec2[1] + fourVec1[2]*fourVec2[2] + fourVec1[3]*fourVec2[3])

########## MAIN FUNCTIONS ##########
# Exits the program if no lhco file is passed in as
# the second argument. Or if there are more than 2
# arguments.
if(len(sys.argv) != 2) :
    print('Usage: python analysis.py <filename>')
    sys.exit()



# Reads in LHCO File and creates Event and Particle
# objects accordingly. The events list is updated accordingly.
events = []
file = open(sys.argv[1], "r")
currEvent = None
for line in file:
    data = line.split()
    if(data[0] == "0") :
        if(currEvent != None) :
            events.append(currEvent)
        currEvent = Event(data)
    elif(data[0] != "#") :
        particle = Particle(currEvent, data)
        currEvent.addFinalStateParticle(particle)

print(events[0].particles[0].px)

