import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math

# 10-2
# 10-3
# 10-6
# 11-1
# 11-2
# 11-6
# 11-10

# *** 10-2 ***
'''
V = molar V

[ I typed these up in python format because you-never-know ]:
Generally, at T_B, 2nd virial coefficient goes to 0, so let's write the R-K equation in terms of virial coefficients:

Z = PV/RT = 1 + B_2V/V + B_3V/(V**2) + ........

R-K: P = (R*T/(V-b))-(a/((T**0.5)*V*(V+b)))
  rewrite: P = ((R*T)/V)*((1/1-b/V)-(a/((T**3/2)*R*(V+b))))
  bionomial expansion, 1/(1-x) == 1+x+x**2+...
  Therefore, P =((R*T)/V)*((1+(b-(a/(R*T)))/V)+...)


let's write the van der waals equation in terms of virial coefficients:

Z = PV/RT = 1 + B_2V/V + B_3V/(V**2) + ........

vdw: P = (R*T/(V-b)) - (a/V**2)
  rewrite: P = ((R*T)/V)*((1/(1-(b/V)-(1/(R*T*V)))))
  bionomial expansion, 1/(1-x) == 1+x+x**2+...
  Therefore, P = ((R*T)/V)*((1+(b-(a/(R*T)))/V)+...)
  
We can see that in both cases, B_2V is equal to b-(a/(R*T))

T_B, B_2V -> 0: b-(a/R*T) = 0
  solve for T_B = a/(b*R)
'''

# *** 10-3 ***
'''
V = molar V

alpha = (1/V)*((dV/dT)_P)

for ideal gas: P = RT/V ... V = RT/P

dV/dT (ideal) = R/P

alpha = (1/V)*(R/P)_P

alpha = R/(V*P)

Note that (V*P)/R = T for ideal gas, so...

alpha = 1/T

'''


# *** 10-6 ***
'''
I got my constants from McQuarrie...
We are using bars

Methane:
a = 2.3026
b = 0.043067

Nitrogen
a = 1.3661
b = 0.038577

'''

def vdw_P(a, b, V, T):
  R = 0.08314
  P = ((R*T)/(V-b)) - (a/(V**2))
  return P
  
def plot_vdw_P_V():
  fig_M = plt.figure(1)
  ax = fig_M.add_subplot(111)
  
  xs_M = []
  ys_M = []
  T_M = 200
  while T_M <= 500:
    xs_M.append(T_M)
    ys_M.append(vdw_P(2.2762, 0.043067, 1, T_M))
    T_M += 0.01
  
  xs_N = []
  ys_N = []
  T_N = 200
  while T_N <= 500:
    xs_N.append(T_N)
    ys_N.append(vdw_P(1.3661, 0.038577, 1, T_N))
    T_N += 0.01
  
  ax.set_title('P(T) van der waals (Methane and Nitrogen)')
  ax.plot(xs_M, ys_M)
  ax.plot(xs_N, ys_N)
  ax.set_xlabel('Temperature (Kelvin)')
  ax.set_ylabel('Pressure (Bar)')
  plt.savefig('graph_methane_nitrogen_vdw.png')

plot_vdw_P_V()


'''
vdw: P = (R*T/(V-b)) - (a/V**2)

or: V**3 - (b + ((R*T)/P))*(V**2) + (a/P)*V - (a*b)/P = 0

'''

def vdw_V(a, b, P, T):
  R = 0.08314
  coeff = [1.0, -(P*b+R*T)/P, a/P,  -a*b/P]
  ans = np.roots(coeff).real
  imag = np.roots(coeff).imag
  absimag = np.absolute(imag)
  index = np.argmin(absimag)
  return (ans[index])

print ('10-6 (b) the density question')
print(1/vdw_V(2.2762, 0.043067, 500, 500))

def calc_Z(P, V, T):
  R = 0.08314
  return (P*V)/(R*T)

def plot_Z(a, b):
  listsOfPlots = []
  tempList = [180, 190, 200, 220, 250]
  for T in tempList:
  
    P = 1
    Zs = []
    Ps = []
    legendList = []
    while (P <= 200):
      V = vdw_V(a, b, P, T)
      Z = calc_Z(P, V, T)
      Zs.append(Z)
      Ps.append(P)
      P+=1
    
    listsOfPlots.append([Ps, Zs])
    legendList.append(str(T)+"K")
  
  fig2 = plt.figure(2)
  ax = fig2.add_subplot(111)
  
  for i in range(len(listsOfPlots)):
    ax.plot(listsOfPlots[i][0], listsOfPlots[i][1], label = (str(tempList[i])+" K"))
  
  ax.set_title("Isothermal Z vs. P")  
  ax.set_xlabel('Pressure (bar)')
  ax.set_ylabel('Z')
  ax.legend(bbox_to_anchor=(0.9, 0.1), loc="lower right", borderaxespad=0.)
  
  plt.savefig('Z vs P.png')
  

plot_Z(2.2762, 0.043067)

# *** 11-1 ***

pH1 = 0.75
pH2 = 0.5

probabilityDict = {
  '[0, 0]' : pH1*pH2,
  '[0, 1]' : pH1*(1-pH2),
  '[1, 0]' : (1-pH1)*(pH2),
  '[1, 1]' : (1-pH1)*(1-pH2)
}

def Probability(coins):
  return probabilityDict[str(coins)]


# *** 11-2 ***

def energyVibDiatomic(hbarOmega, v):
  return hbarOmega*(v+0.5)

def MCstep(kT, hbarOmega, v):
  vprime = max(v + np.random.choice([-1,1]), 0)
  deltaE = energyVibDiatomic(hbarOmega, vprime)-energyVibDiatomic(hbarOmega, v)
  if deltaE <= 0:
    newV = vprime
  else:
    if np.random.random() <= math.exp(-deltaE/kT):
      newV = vprime
    else:
      newV = v
  return newV
  
def runMC(kT, hbarOmega, nEquil, nDataCol):
  v = 0 
  vSamples = []
  for i in range (nEquil):
    v = MCstep(kT,hbarOmega,v)
  for i in range(nDataCol):
    v = MCstep(kT,hbarOmega,v)
    if i%100==0:
      vSamples.append(v)
  return vSamples
  
def MCProbResult(vSamples,vTarget):
  count = 0
  for sample in vSamples:
    if sample == vTarget:
      count +=1
  return float(count/len(vSamples))
  
def convergenceTest():
  kT = 1
  hbarOmega = 1
  nRuns = []
  probs = []
  ntoAdd = 2000
  while ntoAdd < 100000:
    nRuns.append(ntoAdd)
    ntoAdd += 2000
  
  for runs in nRuns:
    probs.append(MCProbResult(runMC(kT, hbarOmega, 100, runs), 1))
    
  fig = plt.figure(3)
  ax = fig.add_subplot(111)
  ax.plot(nRuns, probs, zorder=1, color = "black")
  ax.scatter(nRuns, probs, s=3, zorder=2, color = "red")
  ax.set_xlabel('Number of Runs')
  ax.set_ylabel('p(v=1) at kT/ħω=1')
  ax.set_title("Convergence of Probability as Number of Runs ↑")
  plt.savefig("convergenceTest.png")

convergenceTest()

# *** 11-6 ***
def MCstep2(kT, E, energies):
  Eprime = np.random.choice(energies)
  deltaE = Eprime - E
  if deltaE <= 0:
    newE = Eprime
  else:
    if np.random.random() <= math.exp(-deltaE/kT):
      newE = Eprime
    else:
      newE = E
  return newE
    
  
def runMC2(kT, energies, nEquil, nDataCol):
  E = 0
  Esamples = []
  for i in range(nEquil):
    E = MCstep2(kT, E, energies)
  for i in range (nDataCol):
    E = MCstep2(kT, E, energies)
    if i%100 == 0:
      Esamples.append(E)
  return Esamples
  
def MCProbResult2(Esamples, targetE):
  count = 0
  for sample in Esamples:
    if sample == targetE:
      count +=1
  return count/len(Esamples)


energiesList = [0,0,0,0,0.9427,0.9427]
EsamplesList = runMC2(8.617e-5*5000, energiesList, 100, 100000)
print("\nthe probability of finding Iodine in the excited state:")
print(MCProbResult2(EsamplesList, 0.9427))
  
# *** 11-10 ***

def harmE(k, x):
  return (k*(x**2))/2
  
def MCstep3(kT, k, x):
  xprime = max(x+np.random.uniform(-1,1),0)
  deltaE = harmE(k, xprime) - harmE(k, x)
  if deltaE <= 0:
    newx = xprime
    accepted = True
  else:
    if np.random.random() <= math.exp(-deltaE/kT):
      newx = xprime
      accepted = True
    else:
      newx = x
      accepted = False
  return newx, accepted

def runMC3(kT, k, nEquil, nDataCol):
  x = 0
  xsamples = []
  for i in range(nEquil):
    x, accepted = MCstep3(kT, k, x)
  numaccepted = 0
  for i in range (nDataCol):
    x, accepted = MCstep3(kT, k, x)
    if accepted:
      numaccepted += 1
    if i%100 == 0:
      xsamples.append(x)
  print('\nacceptance fraction is: ' + str(numaccepted/nDataCol))
  return xsamples

def MC3Hist(kT, k, nEquil, nDataCol):
  xsamples = runMC3(kT, k, nEquil, nDataCol)
  npxsamples = np.array(xsamples)
  print('The mean position for the harmonic particle with kT = ' + str(kT) +' is:')
  print(np.mean(npxsamples))
  print('The standard deviation for the harmonic particle is:')
  print(np.std(npxsamples))
  fig = plt.figure(4)
  ax = fig.add_subplot(111)
  ax.hist(xsamples, label = 'kT = ' + str(kT))
  ax.legend(bbox_to_anchor=(0.9, 0.1), loc="lower right", borderaxespad=0.)
  ax.set_title("distribution of particle positions")
  ax.set_ylabel("frequency")
  ax.set_xlabel("x")
  plt.savefig('harmonic particle positions for kT = ' + str(kT) + '.png')
  
MC3Hist(1, 1, 100, 1000) 
MC3Hist(0.1, 1, 100, 1000)
MC3Hist(10, 1, 100, 1000)