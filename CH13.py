import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

# 13-3:

def genCheckConfig(N):
  npoutput = np.zeros((N, N))
  output = npoutput.tolist()
  for i in range(N):
    for j in range(N):
      if i%2 == 0:
        if j%2 ==0:
          output[i][j] = +1
        else:
          output[i][j] = -1
      else:
        if j%2 == 0:
          output[i][j] = -1
        else:
          output[i][j] = +1
      
  return(output)
  
def energyIsing2DPBC(config, J, B):
  energyB = 0
  energyJ = 0
  for row in range(len(config)-1):
    for col in range(len(config)-1):
      energyB+=B*config[row][col]
      if row > 0:
        energyJ += J*config[row-1][col]*config[row][col]
        energyJ += J*config[len(config)-1][col]*config[row][col]
      if row < len(config)-1:
        energyJ += J*config[row+1][col]*config[row][col]
        energyJ += J*config[0][col]*config[row][col]
      if col > 0:
        energyJ += J*config[row][col-1]*config[row][col]
        energyJ += J*config[0][col]*config[row][col]
      if col < len(config)-1:
        energyJ += J*config[row][col+1]*config[row][col]
        energyJ += J*config[row][len(config)-1]*config[row][col]
  return energyB+(energyJ/2)
  
def netMagnetizationPerSpin(config):
  sum = 0
  for i in range(len(config)-1):
    for j in range(len(config)-1):
      sum += config[i][j]
  return sum/(len(config)**2)
  

def MCstepPBC(kT, J, B, config):
  newConfig = config
  [row,col] = np.random.randint(0,len(config),size=2)
  Estart = config[row][col]*B
  if row > 0:
    Estart += J*config[row-1][col]*config[row][col]
    Estart += J*config[len(config)-1][col]*config[row][col]
  if row < len(config)-1:
    Estart += J*config[row+1][col]*config[row][col]
    Estart += J*config[0][col]*config[row][col]
  if col > 0:
    Estart += J*config[row][col-1]*config[row][col]
    Estart += J*config[0][col]*config[row][col]
  if col < len(config)-1:
    Estart += J*config[row][col+1]*config[row][col]
    Estart += J*config[row][len(config)-1]*config[row][col]
    
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
  
  
def testSizes(sizes):
  data = []
  for size in sizes:
    thisSize = []
    config = genCheckConfig(size)
    kT = 1
    while kT <= 4:
      print ("starting size", size, "with temperature", kT)
      meanEPerSpin, heatCapacity, netMagPerSpin, finalConfig = runMC(kT, -1, 0, 10*(size**2), 500*(size**2), 30*(size**2), config, MCstepPBC, energyIsing2DPBC)
      thisSize.append([size, kT, meanEPerSpin, heatCapacity, netMagPerSpin])
      kT += 0.25
    data.append(thisSize)
  return data



def graphSizes(sizes):
  data = testSizes(sizes)
  f, (enax, heatax, magax) = plt.subplots(3)
  enax.set_ylabel("mean energy/spin")
  magax.set_ylabel("abs mean net mag/spin")
  heatax.set_ylabel("heatCap")
  enax.set_xlabel("temperature")
  heatax.set_xlabel("temperature")
  magax.set_xlabel("temperature")
  for size in data:
    en = []
    heat = []
    mag = []
    kT = []
    legendList = []
    for temp in size:
      en.append(temp[2])
      heat.append(temp[3])
      mag.append(temp[4])
      kT.append(temp[1])
    enax.plot(kT, en, label=str(size[0][0]))
    heatax.plot(kT, heat, label=str(size[0][0]))
    magax.plot(kT, mag, label=str(size[0][0]))
    legendList.append(str(size)+"x"+str(size))
  plt.legend(bbox_to_anchor=(0.01, 1.7), loc="center left", borderaxespad=0.)
  plt.savefig("13-3.png")
      
  
  
#graphSizes([4,10,20,30])

# 13-5:

def runMCSus(kT, J, B, nEquil, nDataCol, sampleInterval, config, MCstepFunction, energyFunction):
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
  MsamplesSus = Msamples
  for Msample in MsamplesSus:
    Msample = Msample**2
  return [np.mean(Esamples)/(len(config)**2), np.var(Esamples)/(len(config)**2), np.mean(Msamples), np.mean(MsamplesSus), newConfig]

def testSizesSus(sizes):
  data = []
  for size in sizes:
    thisSize = []
    config = genCheckConfig(size)
    kT = 1
    while kT <= 4:
      print ("starting size", size, "with temperature", kT)
      meanEPerSpin, heatCapacity, netMagPerSpin, netMagPerSpinSq, finalConfig = runMCSus(kT, -1, 0, 10*(size**2), 500*(size**2), 30*(size**2), config, MCstepPBC, energyIsing2DPBC)
      thisSize.append([size, kT, meanEPerSpin, heatCapacity, netMagPerSpin, netMagPerSpinSq])
      kT += 0.25
    data.append(thisSize)
  return data

def findMagSus(data):
  N = data[0]
  kT = data[1]
  M = data[4]
  Msq = data[5]
  return (N**2/kT)*(Msq - M**2)
  

  
def graphSizesSus(sizes):
  f, (ax) = plt.subplots(1)
  data = testSizesSus(sizes)
  for size in data:
    magSus = []
    kT = []
    for temp in size:
      magSus.append(findMagSus(temp))
      kT.append(temp[1])
    ax.plot(kT, magSus, label = size[0][0])
  plt.legend(bbox_to_anchor=(0.9, 0.1), loc="lower right", borderaxespad=0.)
  plt.savefig("13-5.png")
  
  
#graphSizesSus([4,10,20,30])


  
# 13-11 :

def genMixConfig(N):
  npoutput = np.zeros((N, N))
  output = npoutput.tolist()
  for i in range(N):
    for j in range(N):
      if j <= (N/2) - 1:
        output[i][j] = 1
      else:
        output[i][j] = -1
      
  return(output)

def genDataMix(kT, size):
  config = genMixConfig(10)
  meanEPerSpin, heatCapacity, netMagPerSpin, finalConfig = runMC(kT, -1, 0, 10*(size**2), 500*(size**2), 30*(size**2), config, MCstepPBC, energyIsing2DPBC)
  return meanEPerSpin, heatCapacity, netMagPerSpin, finalConfig
  
def graphMix(size):
  f, (enax, heatax, magax) = plt.subplots(3)
  enax.set_ylabel("mean energy/spin")
  magax.set_ylabel("abs mean net mag/spin")
  heatax.set_ylabel("heatCap")
  enax.set_xlabel("temperature")
  heatax.set_xlabel("temperature")
  magax.set_xlabel("temperature")
  kT = 0.5
  en = []
  heat = []
  mag = []
  temps = []
  legendList = []
  while kT <= 4:
    meanEPerSpin, heatCapacity, netMagPerSpin, finalConfig = genDataMix(kT, size)
    en.append(meanEPerSpin)
    heat.append(heatCapacity)
    mag.append(netMagPerSpin)
    temps.append(kT)
    kT += 0.5
  enax.plot(temps, en)
  heatax.plot(temps, heat)
  magax.plot(temps, mag)
  print(np.array(finalConfig))
  plt.savefig("13-11.png")
  
graphMix(6)
  
  