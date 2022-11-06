from ROOT import TFile,TH1F,TH2F
import numpy as np
from matplotlib import pyplot as plt
import uproot
from skhep.visual import MplPlotter as skh_plt
import glob
import pdb



file100 = "../crab100keV_unif.root"
file2 = "../crab2MeV.root"
#f = uproot.open(file2)["ntuple/S2"]
#S2 = f.arrays()
f = TFile(file2)

Nsteps = f.ntuple.S2.GetEntries()
Lgain = TH1F("LG","LEM-gain",40,50.,150.)

oldEvent = 0
ecnt = 0
ecntevt = []
opcnt = 0
opcntevt = []
ypt = 0
yevt = []
print("S2 array is " + str(Nsteps) + " long.")
for step in range(Nsteps):
    if not step%1E6:
        print(str(step/1E6) + " millionth step.")
    S2 = f.ntuple.S2.GetEntry(step)

    event = f.ntuple.S2.Event
    if event == oldEvent:
        if f.ntuple.S2.PID == 11:
            ecnt += 1
        if f.ntuple.S2.PID == -22:
            opcnt += 1
        ypt = (f.ntuple.S2.Y>580) or (f.ntuple.S2.Y<20) or (np.sqrt(f.ntuple.S2.Z**2.+f.ntuple.S2.X**2.) > 130.)
#        if ypt:
#         pdb.set_trace()
        ym = f.ntuple.S2.Y
    else:
        oldEvent = event
#        print (" ypt " + str(ypt))
        ypt = 0
        if not ypt:
            ecntevt.append(ecnt/2000./32.) # 100, per keV
            opcntevt.append(opcnt/2000./32.)
            yevt.append(ym)
#            pdb.set_trace()
#        if event==20:
#            pdb.set_trace()
        ecnt = 0
        opcnt = 0
        ypt = 0


ent = np.array(ecntevt)
opnt = np.array(opcntevt)
ynt = np.array(yevt)

both = []
#both.append(ent)
both.append(opnt)


labelv = np.array(("S2 optPhotons"))
binsz = 2.0
cts, be, er = skh_plt.hist(both,bins=np.arange(50.,150.,binsz),errorbars=False, histtype='step',label=labelv)#,stacked='true')
hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er
pdb.set_trace()

fileout = "S2gain_2MeV_fid"
plt.legend()
plt.title(fileout)
plt.xlabel('optphots/2000/32')
plt.ylabel('entries')
#plt.yscale('log')

plt.savefig(fileout+'.png')
plt.close()
