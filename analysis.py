import sys
import math
import itertools
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1F
from ROOT import gROOT

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
invisibleList = [12, -12, 14, -14, 16, -16]

########## CLASS DEFINITIONS ##########
# Class object that represents one event
class Event:
    def __init__(self, header):
        self.eventNum = int(header[1])
        self.particles = []
        self.keep = True

    def addFinalStateParticle(self, particle) :
        self.particles.append(particle)

    def totalVisiblepT() : 
        sum = [0.0, 0.0]
        for particle in self.particles :
            if(not(particle.pid in invisibleList)) :
                sum[0] += particle.px
                sum[1] += particle.py
        return sum

    def totalVisiblep() : 
        sum = [0.0, 0.0, 0.0]
        for particle in self.particles :
            if(not(particle.pid in invisibleList)) :
                sum[0] += particle.px
                sum[1] += particle.py
                sum[2] += particle.pz
        return sum
    
    def selectionCut(self) :
        self.particles.sort(key=lambda x: x.pt)

        ptj = 20
        pta = 10
        ptl = 10
        etaj = 3
        etal = 2.5
        etaa = 2.5
        drjj = 0.0
        drll = 0.4
        drla = 0.4

        METfind = list(filter(lambda x: x.pid == 12, self.particles))
        lepfind = list(filter(lambda x: x.pid in lightleptonPIDlist, self.particles))
        jetfind = list(filter(lambda x: x.pid in alljetPIDlist, self.particles))

        self.METfind = METfind
        self.lepfind = lepfind
        self.jetfind = jetfind
        # Visible find?

        if(not(len(lepfind) >= 2 and len(jetfind) >= 2 and len(METfind) > 0 and jetfind[0].pt >= ptj and lepfind[0].pt >= ptl)) :
            self.keep = False
            return
        
        twojets = set(itertools.combinations(jetfind, 2))
        twoleps = set(itertools.combinations(lepfind, 2))
        jetandlep = set(itertools.product(jetfind, lepfind))

        for pair in twojets :
            if(not(abs(etaOf(pair[0])) < etaj and abs(etaOf(pair[1])) < etaj and deltaR(pair[0], pair[1]) > drjj)) :
                self.keep = False
                return
        
        for pair in twoleps :
            if(not(abs(etaOf(pair[0])) < etal and abs(etaOf(pair[1])) < etal and deltaR(pair[0], pair[1]) > drll)) :
                self.keep = False
                return

        for pair in jetandlep :
            if(not(deltaR(pair[0], pair[1]) > drla)) :
                self.keep = False
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

        if(self.px > 0 and self.py < 0) :
            self.phi = math.atan(self.py/self.px)
        elif((self.px < 0 and self.py > 0) or (self.px < 0 and self.py < 0)) :
            self.phi = math.atan(self.py/self.px) + math.pi
        elif(self.px > 0 and self.py < 0) :
            self.phi = math.atan(self.py/self.px) + 2*math.pi
        
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

    def threeVector(self) :
        return (self.px, self.py, self.pz)

def fourDotProduct(particle1, particle2) :
    fourVec1 = particle1.fourVector()
    fourVec2 = particle2.fourVector()
    return fourVec1[0] * fourVec2[0] - (fourVec1[1]*fourVec2[1] + fourVec1[2]*fourVec2[2] + fourVec1[3]*fourVec2[3])

def deltaR(particle1, particle2) :
    return math.sqrt(math.pow(etaOf(particle1) - etaOf(particle2), 2) + math.pow(deltaPhi(particle1,particle2),2))

def deltaPhi(particle1, particle2) :
    difference = particle1.phi - particle2.phi
    return min(difference, 2*math.pi - difference, 2*math.pi+difference, key=abs)

def rapOf(particle) :
    return 0.5*math.log10((particle.jmas + particle.pz)/(particle.jmas - particle.pz))

def ptOf(particle) : 
    return particle.pt

def etaOf(particle) :
    return -1*math.log10(math.tan(math.acos(particle.pz/particle.p)/2))

def fourLength(particle) :
    return math.sqrt(math.pow(particle.jmas, 2) - math.pow(particle.pz, 2) - math.pow(particle.py, 2) - math.pow(particle.px, 2))

def fourLengthSq(particle) :
    return math.pow(particle.jmas, 2) - math.pow(particle.pz, 2) - math.pow(particle.py, 2) - math.pow(particle.px, 2)

########## MAIN FUNCTIONS ##########
# Exits the program if no lhco file is passed in as
# the second argument. Or if there are more than 2
# arguments.
if(len(sys.argv) != 3) :
    print('Usage: python analysis.py <filename> <filename>')
    sys.exit()



# Reads in LHCO File and creates Event and Particle
# objects accordingly. The events list is updated accordingly.
events = []
file = open(sys.argv[1], "r")
currEvent = None

for line in file :
    data = line.split()
    if(data[0] == "0") :
        if(currEvent != None) :
            events.append(currEvent)
        currEvent = Event(data)
    elif(data[0] != "#") :
        particle = Particle(currEvent, data)
        currEvent.addFinalStateParticle(particle)
events.append(currEvent)

for event in events :
    event.selectionCut()

bkevents = []
bkfile = open(sys.argv[2], "r")
bkcurrEvent = None

for line in bkfile :
    data = line.split()
    if(data[0] == "0") :
        if(bkcurrEvent != None) :
            bkevents.append(bkcurrEvent)
        bkcurrEvent = Event(data)
    elif(data[0] != "#") :
        particle = Particle(bkcurrEvent, data)
        bkcurrEvent.addFinalStateParticle(particle)
bkevents.append(bkcurrEvent)


for event in bkevents :
    event.selectionCut()
canvas = TCanvas('METs')
example = TFile( 'analysis.root', 'RECREATE' )
histo1 = TH1F('MET', 'MET', 100, 0, 1000)
histo2 = TH1F('METBK', 'MET', 100, 0, 1000)

for event in events :
    if(event.keep and event.METfind[0].pt > 150) :
        histo1.Fill(event.METfind[0].pt)
norm = histo1.GetEntries()
histo1.Scale(1/norm)

for event in bkevents :
    if(event.keep and event.METfind[0].pt > 150) :
        histo2.Fill(event.METfind[0].pt, 1)
norm = histo2.GetEntries()
histo2.Scale(1/norm)

histo1.SetLineColor(2)
histo1.Draw()

histo2.SetLineColor(1)
histo2.Draw('same')
canvas.Update()
example.Write()