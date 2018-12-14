import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import copy

def localEnergyIsing2D(config, J, B, row, col):
  energyB = 0
  energyJ = 0
  energyB+=B*config[row][col]
  if row > 0:
    energyJ += J*config[row-1][col]*config[row][col]
  if row < len(config)-1:
    energyJ += J*config[row+1][col]*config[row][col]
  if col > 0:
    energyJ += J*config[row][col-1]*config[row][col]
  if col < len(config[0])-1:
    energyJ += J*config[row][col+1]*config[row][col]
  return energyB+energyJ


def energyIsing2D(config, J, B):
  energyB = 0
  energyJ = 0
  for row in range(len(config)-1):
    for col in range(len(config[0])-1):
      if row > 0:
        energyJ += J*config[row-1][col]*config[row][col]
      if row < len(config)-1:
        energyJ += J*config[row+1][col]*config[row][col]
      if col > 0:
        energyJ += J*config[row][col-1]*config[row][col]
      if col < len(config[0])-1:
        energyJ += J*config[row][col+1]*config[row][col]
    return energyB+(energyJ/2)

def MCstepConstNv3(kT, J, B, config):
  config = np.array(config) #explicitly ensure array object
  newConfig = copy.deepcopy(config)
  configDidChange = False
  [rowA,colA,rowB,colB] = np.random.randint(0,len(config),size=4)
  spinA=config[rowA,colA]
  spinB=config[rowB,colB]
  if spinA != spinB:
    newConfig[rowA][colA]=spinB
    newConfig[rowB][colB]=spinA
    
    if ((rowA-rowB)**2)+((colA-colB)**2) < 5:
      minRow = max(min(rowA,rowB)-1,0)
      maxRow = min(max(rowA,rowB)+1,len(config)-1)
      minCol = max(min(colA,colB)-1,0)
      maxCol = min(max(colA,colB)+1,len(config)-1)
      Estart = energyIsing2D(config[minRow:maxRow+1,minCol:maxCol+1], J, B)
      Eend = energyIsing2D(newConfig[minRow:maxRow+1,minCol:maxCol+1], J, B)
    else:
      Estart = localEnergyIsing2D(config,J,B,rowA,colA) + localEnergyIsing2D(config,J,B,rowB,colB)
      Eend = localEnergyIsing2D(newConfig,J,B,rowA,colA) + localEnergyIsing2D(newConfig,J,B,rowB,colB)
    
    if Eend < Estart:
      configDidChange=True
    else:
      if np.random.random() <= math.exp(-(Eend-Estart)/kT):
        configDidChange=True
  if configDidChange==True:
    return newConfig
  else:
    return config
    
def insertionAverage(kT,Eint,Eads,config,nTrials):
  theSum = 0
  nInsert = 0
  for i in range(nTrials):
    [row,col] = np.random.randint(0,len(config),size=2)
    if config[row][col]==0:
      nInsert+=1
      newConfig = copy.deepcopy(config)
      newConfig[row][col]=1 
      theSum+=math.exp(-localEnergyIsing2D(newConfig,Eint,Eads,row,col)/kT)
  print("we inserted", nInsert)
  return theSum/nTrials

  
def runWidomInsertion(kT, Eint, Eads, nEquil, nSwap, nInsertionTrials, nRepeats, config, MCStepFunction, insertionFunction):
  newConfig = copy.deepcopy(config)
  theSum = 0
  for i in range(nEquil):
    newConfig = MCStepFunction(kT,Eint,Eads,newConfig)
  for i in range(nRepeats):
    for j in range(nSwap):
      newConfig = MCStepFunction(kT,Eint,Eads,newConfig)
    theSum+=insertionFunction(kT,Eint,Eads,newConfig,nInsertionTrials)
  print("the average is ", theSum/nRepeats)
  print(newConfig)
  return theSum/nRepeats
  
def generateConfig(occupancy, length):
  output = np.zeros((length,length))
  for i in range(length):
    for j in range(length):
      if (i*length)+(j+1) <= occupancy*(length**2):
        output[i][j] = 1
  return output

def langmuir(x, K):
  return K*x/(1+K*x)

def genAdsorptionData(kT, Eint, Eads, boxlength):
  xs = []
  ys = []
  samples = [0.01, 0.05, 0.1, 0.2, 0.5, 0.55, 0.75, 0.8, 0.81, 0.91]
  for occupancy in samples:
    config = generateConfig(occupancy, boxlength)
    xs.append(occupancy/runWidomInsertion(kT, Eint, Eads, 10000, 1000, 200, 200, config, MCstepConstNv3, insertionAverage))
    ys.append(occupancy)
  data = [xs,ys, samples]
  popt, pcov = curve_fit(langmuir, xs, ys)
  print("K for temperature ", kT, "is ", popt[0])
  print(data)
  return data
  
def plotAdsorption(kT):
  data = genAdsorptionData(kT,0.0,-1.0,10)
  fig = plt.figure(3)
  ax = fig.add_subplot(111)
  ax.plot(data[0],data[1], label = str(kT))
  plt.legend(loc = "center right", borderaxespad = 0.)
  plt.savefig("adsoroption(T)legend")
    
    
def plotIsotherms(temps):
  legendList = []
  for temp in temps:
    plotAdsorption(temp)
    legendList.append(temp)

plotIsotherms([ 0.5])

############################13-19###################################


def plotEint(Eint):
  data = genAdsorptionData(0.5, Eint, 1, 10)
  fig2 = plt.figure(2)
  ax2 = fig2.add_subplot(111)
  ax2.plot(data[0],data[1])
  plt.savefig("adsorption(Eint)")
  
def plotEints(Eints):
  for Eint in Eints:
    plotEint(Eint)
    
#plotEints([-5, -4, -3.5, -3, -2.5, -2, -1, -0.5, -0.1])