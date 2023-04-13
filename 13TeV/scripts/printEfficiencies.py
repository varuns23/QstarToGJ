#!/usr/bin/env python

import sys, os, subprocess, string, re
from ROOT import *
from array import array

masses = array('d')
Efficiencies = array('d')
Efficiencies_sigma = array('d')

mass_start = 1.0
mass_step = 0.05
nsteps = 88

for imass in range(0,nsteps+1):
  mass = mass_start + imass * mass_step
  masses.append(mass)

print "mass[] ="
count_lines = 0; 
for i in range(0,len(masses)):
  count_lines = count_lines+1;
  if count_lines > 10:
    print ""
    count_lines = 1;
  print masses[i],",",

print "\n"
print "Efficiencies[] ="

inputEff = array('d', [0.488178, 0.582729, 0.598829, 0.600764, 0.598127, 0.583467])

count_lines = 0;
for ieff in range(0,len(inputEff)-1):
  diff = (inputEff[ieff+1] - inputEff[ieff])/20.0;
  for jeff in range(0,20):
    count_lines = count_lines + 1;
    if count_lines > 89:
       break;
    Efficiencies.append(inputEff[ieff] + jeff * diff)
    Efficiencies_sigma.append(inputEff[ieff] + jeff * diff)

for ieff in range(0, len(Efficiencies)):
  if ieff % 10 == 0:
    print "";
  print  "%0.6f" % Efficiencies[ieff],", ",

print "\n\n"

for i in range(0, len(Efficiencies)):
  Efficiencies_sigma.append( Efficiencies[ len(Efficiencies) - i - 1])

for i in range(0, len( Efficiencies_sigma)):  
  print "%0.6f" % Efficiencies_sigma[i],", ",

