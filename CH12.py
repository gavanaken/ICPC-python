import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

###########################################12-5############################################################
def makeSixbySix():
  sixbysix = [[1,1,1,-1,-1,-1],[1,1,1,-1,-1,-1],[1,1,1,-1,-1,-1],[1,1,1,-1,-1,-1],[1,1,1,-1,-1,-1],[1,1,1,-1,-1,-1],]

  return(np.array(sixbysix))
  
def energyIsing2D(config, J, B):
  energyB = 0
  energyJ = 0
  for row in range(len(config)-1):
    for col in range(len(config)-1):
      if row > 0:
        energyJ += J*config[row-1][col]*config[row][col]
      if row < len(config)-1:
        energyJ += J*config[row+1][col]*config[row][col]
      if col > 1:
        energyJ += J*config[row][col-1]*config[row][col]
      if col < len(config)-1:
        energyJ += J*config[row][col+1]*config[row][col]
    return energyB+(energyJ/2)
  
def netMagnetizationPerSpin(config):
  sum = 0
  for i in range(len(config)-1):
    for j in range(len(config)-1):
      sum += config[i][j]
  return sum/(len(config)**2)
  

def MCstep(kT, J, B, config):
  newConfig = config
  [row,col] = np.random.randint(0,high=6,size=2)
  Estart = config[row][col]*B
  if row > 0:
    Estart += J*config[row-1][col]*config[row][col]
  if row < len(config)-1:
    Estart += J*config[row+1][col]*config[row][col]
  if col > 0:
    Estart += J*config[row][col-1]*config[row][col]
  if col < len(config)-1:
    Estart += J*config[row][col+1]*config[row][col]
    
  Eend = -Estart
  
  if Eend < Estart:
    newConfig[row][col]*=-1
  else:
    if np.random.random() <= math.exp(-(Eend-Estart)/kT):
      newConfig[row][col] *=-1
      
  return newConfig

def runMC(kT, J, B, nEquil, nDataCol, sampleInterval, config, MCstepFunction, energyFunction):
  Esamples = []
  Msamples = []
  newConfig = config
  for i in range(nEquil):
    newConfig = MCstepFunction(kT, J, B, newConfig)
  for i in range(nDataCol):
    newConfig = MCstepFunction(kT,J,B,newConfig)
    if i%sampleInterval == 0:
      Esamples.append(energyFunction(newConfig,J,B))
      Msamples.append(abs(netMagnetizationPerSpin(newConfig)))
  return [np.mean(Esamples)/(len(config)**2), np.var(Esamples)/(len(config)**2), np.mean(Msamples), newConfig]
      
def genData():
  DataList = []
  for i in range(1, 11):
    DataList.append(runMC(i, -1, 0, 10000, 200000, 1000, makeSixbySix(), MCstep, energyIsing2D))
  return DataList
  
Data = genData()
print(Data)

def graphMag(Data):
  xs = []
  ys = []
  for i in range(1,11):
    xs.append(i)
  for datum in Data:
    ys.append(datum[2])
  fig1 = plt.figure(1)
  ax = fig1.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("magGraph")
  
graphMag(Data)

def graphEn(Data):
  xs = []
  ys = []
  for i in range(1,11):
    xs.append(i)
  for datum in Data:
    ys.append(datum[0])
  fig2 = plt.figure(2)
  ax = fig2.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("enGraph")

graphEn(Data)
  
def graphHeatCap(Data):
  xs = []
  ys = []
  for i in range(1,11):
    xs.append(i)
  for datum in Data:
    ys.append(datum[1])
  fig3 = plt.figure(3)
  ax = fig3.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("heatcapGraph")
  
graphHeatCap(Data)

##################################################12-6###########################################################
def gennewData():
  DataList = []
  for i in range(1, 11):
    DataList.append(runMC(i, 1, 0, 10000, 200000, 1000, makeSixbySix(), MCstep, energyIsing2D))
  return DataList
  
newData = gennewData()

graphMag(newData)
graphEn(newData)
graphHeatCap(newData)



##################################################12-7###########################################################
def gennewDataB(T):
  DataList = []
  i = 0
  while i <= 5:
    DataList.append(runMC(T, 1, i, 10000, 200000, 1000, makeSixbySix(), MCstep, energyIsing2D))
    i += 0.5
  return DataList

def graphMagB(Data):
  xs = []
  ys = []
  i = 0
  while i <= 5:
    xs.append(i)
    i += 0.5
  for datum in Data:
    ys.append(datum[2])
  fig1 = plt.figure(4)
  ax = fig1.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("magGraphB")

def graphEnB(Data):
  xs = []
  ys = []
  i = 0
  while i <= 5:
    xs.append(i)
    i += 0.5
  for datum in Data:
    ys.append(datum[0])
  fig2 = plt.figure(5)
  ax = fig2.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("enGraphB")
  
def graphHeatCapB(Data):
  xs = []
  ys = []
  i = 0
  while i <= 5:
    xs.append(i)
    i += 0.5
  for datum in Data:
    ys.append(datum[1])
  fig3 = plt.figure(6)
  ax = fig3.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("heatcapGraphB")

newDataB = gennewDataB(1)
graphMagB(newDataB)
graphEnB(newDataB)
graphHeatCapB(newDataB)

newDataB10 = gennewDataB(10)
graphMagB(newDataB10)
graphEnB(newDataB10)
graphHeatCapB(newDataB10)

##################################################12-12###########################################################

# note: idk how to do this... but i figured out that the periodic boundary index for i is given by 
#(len(config) + (i%len(config)))% len(config) wait this only works in one direction.... 


##################################################12-13###########################################################

def make36():
  l36 = []
  for i in range(36):
    if i%2 == 0:
      l36.append(1)
    else:
      l36.append(-1)
  return l36
  

      
def energyIsing1D(config, J, B):
  energyB = 0
  energyJ = 0
  for i in range(len(config)-1):
      if i > 0:
        energyJ += J*config[i-1]*config[i]
      if i < len(config)-1:
        energyJ += J*config[i+1]*config[i]
  return energyB+(energyJ/2)
  
def netMagnetizationPerSpin1D(config):
  sum = 0
  for i in range(len(config)-1):
    sum += config[i]
  return sum/(len(config)**2)
  

def MCstep1D(kT, J, B, config):
  newConfig = config
  i = np.random.randint(0,high=6)
  Estart = config[i]*B
  if i > 0:
    Estart += J*config[i-1]*config[i]
  if i < len(config)-1:
    Estart += J*config[i+1]*config[i]
  Eend = -Estart
  
  if Eend < Estart:
    newConfig[i]*=-1
  else:
    if np.random.random() <= math.exp(-(Eend-Estart)/kT):
      newConfig[i] *=-1
      
  return newConfig

def runMC1D(kT, J, B, nEquil, nDataCol, sampleInterval, config, MCstepFunction, energyFunction):
  Esamples = []
  Msamples = []
  newConfig = config
  for i in range(nEquil):
    newConfig = MCstepFunction(kT, J, B, newConfig)
  for i in range(nDataCol):
    newConfig = MCstepFunction(kT,J,B,newConfig)
    if i%sampleInterval == 0:
      Esamples.append(energyFunction(newConfig,J,B))
      Msamples.append(abs(netMagnetizationPerSpin1D(newConfig)))
  return [np.mean(Esamples)/(len(config)**2), np.var(Esamples)/(len(config)**2), np.mean(Msamples), newConfig]
      
def genData1D():
  DataList = []
  for i in range(1, 11):
    DataList.append(runMC1D(i, -1, 0, 10000, 200000, 1000, make36(), MCstep1D, energyIsing1D))
  return DataList
  
data1D = genData1D()

def graphMag1D(Data):
  xs = []
  ys = []
  for i in range(1,11):
    xs.append(i)
  for datum in Data:
    ys.append(datum[2])
  fig1 = plt.figure(7)
  ax = fig1.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("magGraph1D")
  
graphMag1D(data1D)

def graphEn1D(Data):
  xs = []
  ys = []
  for i in range(1,11):
    xs.append(i)
  for datum in Data:
    ys.append(datum[0])
  fig2 = plt.figure(8)
  ax = fig2.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("enGraph1D")

graphEn1D(data1D)
  
def graphHeatCap1D(Data):
  xs = []
  ys = []
  for i in range(1,11):
    xs.append(i)
  for datum in Data:
    ys.append(datum[1])
  fig3 = plt.figure(9)
  ax = fig3.add_subplot(111)
  ax.plot(xs, ys)
  plt.savefig("heatcapGraph1D")
  
graphHeatCap1D(data1D)
