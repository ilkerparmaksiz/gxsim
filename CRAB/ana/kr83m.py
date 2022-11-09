import numpy as np
from ROOT import TFile,TH1F,TCanvas
from matplotlib import pyplot as plt
from skhep.visual import MplPlotter as skh_plt
import pdb

filename = "../gam_41keV.root"

f = TFile(filename)

f.ntuple.S1.Draw("Event>>hS1b(100,0.,100.)")

histo = []
for ii in range(f.hS1b.GetNbinsX()):   
    histo.append(f.hS1b.GetBinContent(ii))

#pdb.set_trace()
cts, be, er = skh_plt.hist(histo[1:],bins=np.arange(0.,1950.,10.),errorbars=False, histtype='step')



plt.text(100,20.,"mean/stdev, meanEnergy [keV]: " + str(np.round(np.mean(histo[1:]),2)) + "/" + str(np.round(np.std(histo[1:]),3)) + ", "  + str(np.round(np.mean(histo[1:]) *22.6/1000.,3)))
fileout = "S1_41keV"
#plt.legend()
plt.title(fileout)
plt.xlabel('S1 thermales')
plt.ylabel('entries')
#plt.yscale('log')                                                                                                                                         

plt.savefig(fileout+'.png')
plt.close()


f.ntuple.S2.Draw("Event>>hS2b(100,0.,100.)")

histo = []
for ii in range(f.hS2b.GetNbinsX()):   
    histo.append(f.hS2b.GetBinContent(ii))


cts, be, er = skh_plt.hist(histo[1:],bins=np.arange(190000.,210000.,1000.),errorbars=False, histtype='step')


plt.text(190000,20.,"mean/stdev, meanEnergy [keV]: " + str(np.round(np.mean(histo[1:]),2)) + "/" + str(np.round(np.std(histo[1:]),3)) + ", "  + str(np.round(np.mean(histo[1:]) /110*22.6/1000.,3)))
fileout = "S2_41keV"
#plt.legend()
plt.title(fileout)
plt.xlabel('S2 thermales')
plt.ylabel('entries')
#plt.yscale('log')                                                                                                                                         

plt.savefig(fileout+'.png')
plt.close()


c1 = TCanvas()
f.ntuple.S2.Draw("X:Z","Event<3","colz")

c1.SaveAs(fileout+"S2-plane"+".png")

