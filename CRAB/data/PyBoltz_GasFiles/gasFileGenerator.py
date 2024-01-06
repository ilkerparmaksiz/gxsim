import time
import argparse
from PyBoltz.OdieRun import *
import numpy as np
print("")
print("PyBoltz, adapted by B. Al. Atoum, A.D. McDonald, and B.J.P. Jones")
print("  from the original FORTRAN code MagBoltz by S. Biagi")

print("The PyBoltz to Garfield++ converter, Odie, written by A. B. Cudd.")
print("Parser Added by Ilker Parmaksiz on Nov 3 2023")
print("-----------------------------------------------------------------")

parser=argparse.ArgumentParser(prog='PyBoltzGasFileCreator',description="Generate Gas Files to be used with garfieldpp",epilog="https://github.com/UTA-REST/PyBoltz.git",fromfile_prefix_chars='@')

## Handling Arguments
parser.add_argument('-g',nargs="+",type=str,help="Add Gasess you want such as Ar,Xe,CH4 such as G1 G2 G3")
parser.add_argument('-f',nargs="+",type=float,help="add Fraction of the gasses such as F1 F2 F3")
parser.add_argument('-o',type=str,help="Output File Name")
parser.add_argument('-p',type=float,help="Pressure in bars")
parser.add_argument('-e',type=int,nargs="+",help="Electrical field Min Max PointSize Ex 100 1000 100 (100 points between 100 and 1000")
parser.add_argument('-c',type=int,help="Adjust Max amount of collision by changing c where MaxCollision is c*4e6 ")
args=parser.parse_args()
print(args)
if(args.g==None or args.f==None  or args.p==None or args.e==None or args.c==None):
    print("One of the variables is empty \n")
    print(f" -g {args.g} \n -f {args.f} \n -c {args.c} \n -p {args.p} \n -e {args.e}")
    exit(0)

if(args.o==None):
    if(len(args.g)>1):
        FileName=args.g[0]+"_"+args.g[1]+"_"+str(args.f[0])+"_"+str(args.f[1])+"_"+str(args.p)+".gas"
    else:
        FileName=args.g[0]+"_"+str(args.f[0])+"_"+str(args.p)+".gas"

else:
    FileName=args.o
## Electrical Fields 
if(args.e[0]>args.e[1] or args.e[2]>args.e[1]):
    print("Electrical Fields are not Assigned Correctly\n")
    print(f" -e {args.e} \n")
    print("Set [LowerLimit HigherLimit PointSize] in V/cm without units ")
    exit(0)





#Set up helper object
Odie = OdieRun()

# Configure settings for our simulation
base_settings = {
    'Gases':args.g,
    'Fractions': args.f,
    'Max_collisions': 4e6*args.c,
    'EField_Vcm': 100,
    'Max_electron_energy': 0,
    'Temperature_C': 23,
    'Pressure_Torr': 750.062*args.p,
    'BField_Tesla': 0,
    'BField_angle': 0,
    'Angular_dist_model': 2,
    'Enable_penning': 0,
    'Enable_thermal_motion': 1,
    'ConsoleOutputFlag': 1
}

t_start = time.time()

#Load base settings into Odie, this is required for
#running on a grid.
Odie.LoadSettings(base_settings, PrintSettings=True)

#Odie also can read in the settings from a JSON file
#Odie.LoadSettings(JSONFileName="input.json", PrintSettings=True)

#Generate output for a grid of possible EFields
#in four steps from 50 V/cm to 200 V/cm on a linear scale
GridOutput = Odie.GenerateGasGrid(args.e[0], args.e[1], args.e[2], LogScale=False)

#Generate output for a grid of possible EFields
#in four steps from 50 V/cm to 200 V/cm on a log scale,
#and in three steps from 0 to 1 Telsa for the BField (linear scale)
#GridOutput = Odie.GenerateGasGrid(50, 200, 4, LogScale=True, 0, 1, 3)

#Finally, a grid of EFields, BFields, and E-B Angles can be generated
#GridOutput = Odie.GenerateGasGrid(
#    minE=50, maxE=200, nE=4, LogScale=False, minB=0, maxB=1, nB=3, minA=0, maxA=90, nA=3
#)

t_end = time.time()

#Write output to Garfield++ style gas file
Odie.WriteGasFile(FileName, GridOutput)

print("Simulation time: {}\n".format(t_end - t_start))

#The output for each grid point is indexed by a tuple of the EField,
#E-B Angle, and BField. In general this loop will iterate over all points.
print("Printing (some) gas properties...")
for e in Odie.GridSettings['EFields']:
    for b in Odie.GridSettings['BFields']:
        for a in Odie.GridSettings['EBAngles']:
            Output = GridOutput[(e, b, a)]

            print("\nE={} V/cm, B={} T, A={}".format(e, b, a))
            print(
                "Vz: {:.3f} +/- {:.3f}".format(
                    Output['Drift_vel'].val[2], Output['Drift_vel'].err[2]
                )
            )
            print("DL: {:.3f} +/- {:.3f}".format(Output['DL'].val, Output['DL'].err))
            print("DT: {:.3f} +/- {:.3f}".format(Output['DT'].val, Output['DT'].err))
