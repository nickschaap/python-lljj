import sys
import math
import itertools
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1F, THStack, TImage
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
        self.cut = False

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
        self.p2 = float(math.pow(self.px, 2) + math.pow(self.py, 2) + math.pow(self.pz, 2))
        self.p = float(math.sqrt(self.p2))

        self.jmas = math.sqrt(math.pow(self.jmas,2) + self.p2)
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

def fourLengthVec(vector) :
    return math.sqrt(math.pow(vector[0], 2) - math.pow(vector[1], 2) - math.pow(vector[2], 2) - math.pow(vector[3], 2))

def fourLengthSq(particle) :
    return math.pow(particle.jmas, 2) - math.pow(particle.pz, 2) - math.pow(particle.py, 2) - math.pow(particle.px, 2)


def eventCount(events) :
    return len(list(filter(lambda x: x.keep and not x.cut, events)))

########## MAIN FUNCTIONS ##########
# Exits the program if no lhco file is passed in as
# the second argument. Or if there are more than 2
# arguments.
if(len(sys.argv) != 4) :
    print('Usage: python analysis.py <signal> <background> <background2>')
    sys.exit()



# Reads in LHCO File and creates Event and Particle
# objects accordingly. The events list is updated accordingly.
def readFileAndPresel(file) :
    events = []
    file = open(file, "r")
    currEvent = None
    for line in file :
        data = line.split()
        if(data[0] == "0") :
            if(currEvent != None) :
                currEvent.selectionCut()
                events.append(currEvent)
            currEvent = Event(data)
        elif(data[0] != "#") :
            particle = Particle(currEvent, data)
            currEvent.addFinalStateParticle(particle)
    currEvent.selectionCut()
    events.append(currEvent)
    return events

# def showEvents(signal, background, cutVar, cutFunc) :
#     canvas = TCanvas(cutVar,cutVar,10,10,500,500)
#     histo1 = TH1F(cutVar + '_SIGNAL', cutVar, 100, 0, 1000)
#     histo2 = TH1F(cutVar + '_BACKGROUND', cutVar, 100, 0, 1000)
#     for event in signal :
#         if(event.keep) :
#             histo1.Fill(cutFunc(event))
#     for event in background :
#         if(event.keep) :
#             histo2.Fill(cutFunc(event))
#     histo1.SetLineColor(2)
#     histo2.SetLineColor(1)
#     histo2.DrawNormalized()
#     histo1.DrawNormalized('same')
#     canvas.Update()
#     canvas.Write()

def showEvents(eventLists, cutVar, cutFunc, showAll) :
    canvas = TCanvas(cutVar,cutVar,10,10,500,500)
    currhisto = None
    count = 1
    for eventList in eventLists :
        currhisto = TH1F(cutVar + '_' + str(count), cutVar, 100, 0, 1000)
        for event in eventList :
            if(event.keep and (showAll or not event.cut)) :
                currhisto.Fill(cutFunc(event))
        count = count + 1
        currhisto.SetLineColor(count)
        if (count == 1) :
            currhisto.DrawNormalized()
        else :
            currhisto.DrawNormalized('same')
    canvas.Update()
    canvas.Write()

def addCut(eventLists, cutFunc, min, max) :
    for eventList in eventLists :
        for event in eventList :
            if(event.keep) :
                cut = cutFunc(event)
                if(cut >= min and cut <= max) :
                    event.cut = False
                else :
                    event.cut = True
            else :
                event.cut = True

signal = readFileAndPresel(sys.argv[1])
background = readFileAndPresel(sys.argv[2])
background2 = readFileAndPresel(sys.argv[3])
allEvents = [signal, background, background2]

def invJJ(event) :
    p1 = event.jetfind[0].fourVector()
    p2 = event.jetfind[1].fourVector()
    p3 = (p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2], p1[3]+p2[3])
    print(p3)
    return fourLengthVec(p3)


rootfile = TFile( 'analysis.root', 'RECREATE' )
# addCut(signal, background + background2, lambda x: invJJ(x), 0, 150)
# showEvents(allEvents, 'Mjj', lambda x: invJJ(x), False)


# MET Cut
addCut(allEvents, lambda x: x.METfind[0].pt, 175, 100000)
showEvents(allEvents, 'MET', lambda x: x.METfind[0].pt, True)

# Mll Cut
def invLL(event) :
    p1 = event.lepfind[0].fourVector()
    p2 = event.lepfind[1].fourVector()
    p3 = (p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2], p1[3]+p2[3])
    return fourLengthVec(p3)


addCut(allEvents, lambda x: invLL(x), 175, 100000)
showEvents(allEvents, 'Mll', lambda x: invLL(x), True)

# Ptl Cut
addCut(allEvents, lambda x: ptOf(x.lepfind[0]), 200, 100000)
showEvents(allEvents, 'Ptl', lambda x: ptOf(x.lepfind[0]), True)


def getSignificance(signal, background) :
    weightSig = 1.01*math.pow(10,-4)*3000*math.pow(10,3)/10000
    weightBK = 2.7*math.pow(10,-3)*3000*math.pow(10,3)/10000
    weightBK2 = 0.0147*3000*math.pow(10,3)/10000
    S = eventCount(signal) * weightSig
    print(S)
    B = eventCount(background) * weightBK
    print(B)
    B2 = eventCount(background2) * weightBK2
    print(B2)
    print('Significance: ', S/math.sqrt(S+ B + B2))

print('Signal: ', eventCount(signal))
print('Background: ', eventCount(background))
print('Background 2: ', eventCount(background2))
getSignificance(signal, background)
rootfile.Write()