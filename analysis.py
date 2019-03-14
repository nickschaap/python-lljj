import sys
import math

########## PID DICTIONARY ##########
pids = { "PIDu": 2, "PIDubar" : -2, "PIDd" : 1, "PIDdbar" : -1, "PIDs" : 3, "PIDsbar" : -3, "PIDc" : 4,
    "PIDcbar" : -4, "PIDb" : 5, "PIDbbar" : -5, "PIDt" : 6, "PIDtbar" : -6, "PIDeminus" : 11, "PIDeplus" : -11,
    "PIDnue" : 12, "PIDnuebar" : -12, "PIDmuminus" : 13, "PIDmuplus" : -13, "PIDnumu" : 14, "PIDnumubar" : -14,
    "PIDtauminus" : 15, "PIDtauplus" : -15, "PIDnutau" : 16, "PIDnutaubar" : -16, "PIDgluon" : 21, "PIDphoton" : 22,
    "PIDZ" : 23, "PIDWplus" : 24, "PIDWminus" : -24, "PIDh" : 25, "PIDVz" : 9000001, "PIDWLp" : 9000002, "PIDWLm" : -9000002,
    "PIDWRp" : 9000003, "PIDWRm" : -9000003, "PIDZprime" : 9000004, "PIDAprime" : 9000005, "PIDeRm" : 9000011,
    "PIDeRp" : -9000011, "PIDmuRm" : 9000013, "PIDmuRp" : -9000013, "PIDtauRm" : 9000015, "PIDtauRp" : -9000015,
    "PIDNRe" : 9000012, "PIDNRebar" : -9000012, "PIDNRmu" : 9000014, "PIDNRmubar" : -9000014, "PIDNRtau" : 9000016,
    "PIDNRtaubar" : -9000016 }

lightleptonPIDlist = [-13, 13, -11, 11]
alljetPIDlist = [21, 5, 15, -15]

########## CLASS DEFINITIONS ##########
# Class object that represents one event
class Event:
    def __init__(self, header):
        self.eventNum = int(header[1])
        self.particles = []

    def addFinalStateParticle(self, particle) :
        self.particles.append(particle)

    def selectionCut(self) :
        self.particles.sort(key=lambda x: x.pt)

        METfind = list(filter(lambda x: x.pid == 12, self.particles))
        lepfind = list(filter(lambda x: x.pid in lightleptonPIDlist, self.particles))
        jetfind = list(filter(lambda x: x.pid in alljetPIDlist, self.particles))

        print("MET\n")
        for particle in METfind :
            print(particle.pid, "\n")

        print("lep\n")
        for particle in lepfind :
            print(particle.pid, "\n")

        print("jet\n")
        for particle in jetfind :
            print(particle.pid, "\n")


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
        pid = None
        if (self.typ == 0) :
            pid = 22
        elif(self.typ == 1 and self.ntrk == 1) :
            pid = -11
        elif(self.typ == 1 and self.ntrk == -1) :
            pid = 11
        elif(self.typ == 2 and self.ntrk == 1) :
            pid = -13
        elif(self.typ == 2 and self.ntrk == -1) :
            pid = 13
        elif(self.typ == 3 and self.ntrk > 0) :
            pid = -15
        elif(self.typ == 3 and self.ntrk < 0) :
            pid = 15
        elif(self.typ == 4 and self.btag > 0) :
            pid = 5
        elif(self.typ == 4 and self.btag == 0) :
            pid = 21
        elif(self.typ == 6) :
            pid = 12
        self.pid = pid
    
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

events[0].selectionCut()

