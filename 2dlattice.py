import numpy as np
import matplotlib.pyplot as plt
from math import pi
import imageio
import scipy.stats as scp
from random import random, gauss
from Tools import acf, getCircle
import csv
'''**************************************************************************'''
xdim = 64
ydim = 64
omegaMu = 5
omegaSig = 1
k = 5
tStep = 1/15
nFrame = 500
save = False
noise = 0
'''**************************************************************************'''

omegaN = np.random.normal(omegaMu,omegaSig,(xdim,ydim))
omegaC = omegaN.copy()
theta = np.random.uniform(0,2*pi,(xdim,ydim))

save_omega = []
save_theta = []
cor = []

def incTheta(i,j):
    theta[i][j]+= omegaC[i][j]*tStep
    if theta[i][j] > 2*pi:
        theta[i][j]-=2*pi
    if theta[i][j] < 0:
        theta[i][j]+=2*pi

def getAdjacent(lattice,i,j): # gets all adjacent coords to (i,j)
    coords = []
    if i != 0:
        coords.append(lattice[i-1][j])
    else:
        coords.append(lattice[xdim-1][j])
    if i != xdim-1:
        coords.append(lattice[i+1][j])
    else:
        coords.append(lattice[0][j])
    if j != 0:
        coords.append(lattice[i][j-1])
    else:
        coords.append(lattice[i][ydim-1])
    if j != ydim-1:
        coords.append(lattice[i][j+1])
    else:
        coords.append(lattice[i][0])
    return coords

def runAS(time):
    crr = 0
    tsin = 0
    tcos = 0
    tact = np.zeros((xdim,ydim))
    for i in range(xdim):
        for j in range(ydim):
            couple = 0
            ac = 0
            temp = theta[i][j]
            crr+=np.exp(complex(0,theta[i][j]))
            tsin += np.sin(theta[i][j])
            tcos += np.cos(theta[i][j])
            for a in getAdjacent(theta,i,j):
                couple += np.sin(a-temp) + omegaSig*gauss(0,noise)
            omegaC[i][j] = omegaN[i][j] + (k/4)*couple
    for i in range(xdim):
        for j in range(ydim):
            incTheta(i,j)
    save_omega.append(omegaC.copy())
    save_theta.append(theta.copy())
    aphi = np.arctan2(tsin,tcos)
    Co = (crr/(xdim*ydim))/np.exp(complex(0,aphi))
    cor.append(Co.real)

def runIS(time):
    for i in range(xdim):
        for j in range(ydim):
            temp = theta[i][j]
            couple = 0
            for k in range(xdim):
                for l in range(ydim):
                    r2 = (i-k)**2 + (j-l)**2
                    if r2!=0:
                        couple+= np.sin(theta[k][l]-temp)/r2
            omegaC[i][j] = omegaN[i][j]+ (k/(xdim*ydim))*couple
    for i in range(xdim):
        for j in range(ydim):
            incTheta(i,j)
    save_omega.append(omegaC.copy())
    save_theta.append(theta.copy())

def getS(lat):
    d = list(np.reshape(lat.copy(),xdim*ydim))
    d.sort()
    x = np.linspace(d[0]+((d[-1]-d[0])/251),d[-1],250)
    i = 0
    pd = np.ones(250)
    oMin = d[0]
    oMax = d[-1]
    for c in x:
        while d[0] < c:
            pd[i]+=1
            temp = d.pop(0)
        i+=1
    Z = np.sum(pd)
    tS = np.zeros((xdim,ydim))
    for i in range(xdim):
        for j in range(ydim):
            tS[i][j] = -1*(pd[int(249*(lat[i][j]-oMin)/(oMax-oMin))]/Z)*np.log(pd[int(249*(lat[i][j]-oMin)/(oMax-oMin))]/Z)
    return tS,pd

def stdR(lat):
    stdd = []
    for i in range(xdim):
        for j in range(ydim):
            temp = []
            for r in range(min(xdim//2,ydim//2)):
                coords = getCircle(xdim//2,ydim//2,r,bounds = (0,xdim-1))
                tot = 0
                for c in coords:
                    tot+=(lat[c]-lat[i][j])**2
                temp.append(tot**0.5/len(coords))
            stdd.append(temp)
    std = []
    for i in range(len(stdd[0])):
        temp = []
        for dd in stdd:
            temp.append(dd[i])
        std.append(np.average(temp))
    return std

save_s = []
save_stdT = []
save_stdO = []
for s in range(1):
    omegaN = np.random.normal(omegaMu,omegaSig,(xdim,ydim))
    omegaC = omegaN.copy()
    theta = np.random.uniform(0,2*pi,(xdim,ydim))
    save_theta.clear()
    save_omega.clear()
    for i in range(nFrame):
        print(i)
        runAS(3)
    entropy = []
    for OC in save_omega:
        tS,pd = getS(OC)
        entropy.append(np.sum(tS))
    save_s.append(entropy)
    std_theta = []
    for s in save_theta:
        std_theta.append(scp.circstd(s))
    save_stdT.append(std_theta)
    std_omega = []
    for s in save_omega:
        std_omega.append(np.std(s))
    save_stdO.append(std_omega)

ent = []
tnt = []
ont = []
for i in range(nFrame):
    stemp = []
    ttemp = []
    otemp = []
    for e in save_s:
        stemp.append(e[i])
    ent.append(np.average(stemp))
    for t in save_stdT:
        ttemp.append(t[i])
    tnt.append(np.average(ttemp))
    for o in save_stdO:
        otemp.append(o[i])
    ont.append(np.average(otemp))

'''with open('2d{}entropy.csv'.format(k),'w') as file:
    writer = csv.writer(file)
    writer.writerow(ent)
    file.close()
with open('2d{}stdTheta.csv'.format(k),'w') as file:
    writer = csv.writer(file)
    writer.writerow(tnt)
    file.close()
with open('2d{}stdOmega.csv'.format(k),'w') as file:
    writer = csv.writer(file)
    writer.writerow(ont)
    file.close()'''

plt.plot(ent)
plt.show()

std_omega = []

for s in save_omega:
    std_omega.append(np.std(s))

def saveFrame(lattice,f,saveList):
    lwidth,lheight = np.shape(lattice[f])
    temp = np.zeros((512,512),dtype=np.uint8)
    maxV = np.max(lattice)
    cpixel = 0
    for i in range(lwidth):
        for j in range(lheight):
            cpixel = lattice[f][i][j]/maxV
            for k in range(512//lwidth):
                for l in range(512//lheight):
                    temp[i*(512//lwidth)+k][j*(512//lheight)+l] = (cpixel*255)
                 
    saveList.append(temp)

if save:
    print('saving...')
    for i in range(nFrame):
        print(i)
        saveFrame(save_theta,i,final_theta)
        saveFrame(save_omega,i,final_omega)
    imageio.mimwrite('{} by {} k {} Theta.mp4'.format(xdim,ydim,k),final_theta,fps=15)
    imageio.mimwrite('{} by {} k {} Omega.mp4'.format(xdim,ydim,k),final_omega,fps=15)



plt.plot(entropy)
plt.show()
tS,pd = getS(omegaC)
std = stdR(omegaC)
ocor = acf(tS)
plt.plot(ocor)
plt.plot(std)

