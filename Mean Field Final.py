import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.lines import Line2D
from random import gauss, random, choice
from math import pi, sin, exp, cos
import csv

'''****************************Enter Input Here******************************'''
num = 50 # The number of oscillators in the system (N)
kMu = 1.25 # The mean coupling constant (K)
kSig = 0 # The standard dev for the coupling constant
omegaMu = 5 # The mean angular frequency
omegaSig = 0.5 # The standard dev for the angular frequency
nFrame = 400 # How many frames it will run (1 frame = 1/30 seconds)
animate = False
save = False
'''**************************************************************************'''

# This is the class to represent a system of coupled oscillators
class Osc(object):
            
    num #number of oscillators
    kMu #mean coupling constant
    kSig #scoupling constant standard deviation
    omegaMu #mean angular frequency (rad/sec)
    omegaSig #angular frequency standard deviation
    kArray = np.zeros([num,num]) #coupling array for individual coupling
    oscArray = np.zeros([5,num]) #data array
    timeArray = [] # holds theta(t) values
    rCo = complex(0,0) #coherence/phase
    aPhi = 0
    aOmega = 0
    kTot = 0
    pTot = 0

    def __init__(self,num = 25,kMu = 1,kSig = 0,omegaMu = 5, omegaSig = 0.5):
        self.num = num
        self.kMu = kMu
        self.kSig = kSig
        self.omegaMu = omegaMu
        self.omegaSig = omegaSig
        
    for i in range(num): #filling in data
        oscArray[0][i] = random()*2*pi  #Theta (rad)
        #oscArray[0][i] = choice([4.25,4.75])
        oscArray[1][i] = gauss(omegaMu,omegaSig) #Omega (rad/sec)
        oscArray[2][i] = oscArray[1][i] #Natural Omega (rad/sec)
        oscArray[3][i] = 0.001 #Mass (kg)
        oscArray[4][i] = 0.01 #radius (m)
        timeArray.append([oscArray[0][i]]) #Fill in time array with lists
        
        for j in range(num): #fill the coupling array with gaussian distribution
            if i != j:
                kArray[i][j] = gauss(kMu,kSig)
                
    def inRad(self,rad, inc): #method to check period boundaries
        nRad = rad+inc
        if nRad > 2*pi:
            nRad -= 2*pi
        elif nRad < 0:
            nRad += 2*pi
        return nRad
    
    def runTime(self, timeStep): #runs the kuramoto equation once
        coherence = 0
        noise = gauss(0,0.1)
        sinT = 0
        cosT = 0
        self.kTot = 0
        self.pTot = 0
        tOmega = 0
        if self.aPhi > 2*pi:
            self.aPhi-=2*pi
        elif self.aPhi < 0:
            self.aPhi+=2*pi
        for i in range(self.num):
            couple = 0
            tOmega += self.oscArray[1][i]
            coherence += np.exp(complex(0,self.oscArray[0][i]))
            sinT += sin(self.oscArray[0][i])
            cosT += cos(self.oscArray[0][i])
            self.timeArray[i].append(self.oscArray[0][i])
            self.kTot += 0.5* self.oscArray[3][i]*self.oscArray[4][i]**2*self.oscArray[1][i]**2
            for j in range(self.num):
                couple += sin(self.oscArray[0][j]-self.oscArray[0][i])
                if i != j:
                    self.pTot += 1-cos(self.oscArray[0][i]-self.oscArray[0][j])
            self.oscArray[1][i] = self.oscArray[2][i]+ (self.kArray[i][j]/self.num)*couple
            self.oscArray[0][i] = self.inRad(self.oscArray[0][i],
                                        self.oscArray[1][i]*timeStep)
        self.rCo = coherence/self.num
        self.aPhi = np.arctan2(sinT,cosT)
        self.pTot = self.pTot * self.kMu/(2*self.num)
        self.aOmega = tOmega/self.num
            
    def getTheta(self):
        theta = self.oscArray[[0],:]
        return theta
    
    def getOmega(self):
        omega = self.oscArray.copy()[[1],:]
        return omega,

    def getOstd(self):
        return np.std(self.oscArray[1])

    def getR(self):
        return self.rCo

    def getE(self):
        return self.kTot +self.pTot

if animate:
    osc = Osc(num,kMu,kSig,omegaMu,omegaSig) #optional to fill in params here (num, kMu, kSig, omegaMu, omegaSig)
    rArray = []
    eArray = []
    oArray = []
    osArray = []
    
    fig = plt.figure(1)
    ax = plt.axes(xlim=(0,num),ylim=(0,2*pi))
    plot, = ax.plot([],[],'bo')
    plt.ylabel('Theta (radians)')
    plt.title('Kuramoto System')

    def init():
        plot.set_data([],[])
        return plot,

    def animate(i):
        print(i)
        osc.runTime(1/30)
        rCo = osc.getR()
        r= rCo/np.exp(complex(0,osc.aPhi))
        rArray.append(r.real)
        #eArray.append(osc.getE())
        #oArray.append(osc.aOmega)
        osArray.append(osc.getOstd())
        x = range(osc.num)
        y = osc.getTheta()
        plot.set_data(x,y)
        return plot,

    anim = animation.FuncAnimation(fig,animate, init_func=init, frames = nFrame,
                             interval=20,blit=True)
    plt.show()
    if save:
        anim.save('kuramoto.mp4', fps=30,extra_args=['-vcodec', 'libx264'])

    
if __name__ == '__main__':
    cor = []
    std = []
    osc = Osc(num,kMu,kSig,omegaMu,omegaSig)
    for i in range(nFrame):
        osc.runTime(1/30)
        r = osc.getR()/np.exp(complex(0,osc.aPhi))
        cor.append(r.real)
        std.append(osc.getOstd())
    '''with open('n{}choice45.csv'.format(num,omegaSig),'w') as file:
        writer = csv.writer(file)
        writer.writerow(cor)
        writer.writerow(std)
        file.close()'''
    plt.plot(cor)
    plt.plot(std)
    plt.show()
