#Details:
#Particles have different sizes (Req=1+a*rand()) a~0.1
#Particles reaching the right edge are reinserted in the first left eighth with the same
#speed v, direction n, and y-position.
#good noise value 1.0
#good density 1.0
#good valeu for wall_osc~0.01
#input_data = sys.argv[0:4]
#tau=float(input_data[1])
#division_time=int(input_data[2])
#initial_state=int(input_data[3]) 0-start from t=0; 1-continue    

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import math as math
import random as rand
import os
import sys
from numba import jit
from mydefs import initial_0, initial_1, outfile_init, save_state, make_graph, state_init

#Particle class definition
class particle:
   # noise=0.6 #original value
    #noise=0.6
#    noise_T=10.0
    v0 = 0.1
    mu=1.0
    Frep=30.0
    #Fadh=0.75 #original value
#    Fadh=10.0
    #Fadh=0.1
    #Req=5./6. # original work by Szabo
    #R0=1.0
#    Req=1.0
    R0=3.0
    def __init__(self, x, y, vx, vy, ident, Raio_equilibrio,fadh,noise):
        self.r = np.array([x,y])
        self.v =  np.array([vx*self.v0,vy*self.v0])
        self.theta = np.arctan2(vy,vx)
        self.n = np.array([np.cos(self.theta),np.sin(self.theta)])
        self.ident = ident
        self.Mybox = int((self.r[0]+L[0])/lbox)+nb[0]*int((self.r[1]+L[1])/lbox)
        self.Force =np.array([0.,0.])
        self.Req = 2+Raio_equilibrio
        self.cross = 0
        self.cross2 = 0
        self.Fadh = fadh
        self.noise = noise

    def mybox(self): #Each particle calculates the box it is in
        if np.isnan(self.r[0]) == True or np.isnan(self.r[1]) == True :
            print(self.ident, self.r, self.Force)
            exit()
        j=int((self.r[0]+L[0])/lbox)+nb[0]*int((self.r[1]+L[1])/lbox)
        return j

    def set_new_born_neigh(self):
        newbox=self.Mybox
        #Create lists with the indexes of the neighboring boxes to change the particles there too
        #This particle will be in the lists of three boxes above its box and one in the left
        new_neigh_box_list=[]
        if newbox%nb[0] != 0 :                              #test if not first column
            new_neigh_box_list.append(newbox-1)             #add to left box list  
        if newbox%nb[0] != 0 and newbox > nb[0] :           #test if not first column and not first line
            new_neigh_box_list.append(newbox-nb[0]-1)       #add left box above list
        if newbox > nb[0] :                                 #test if not first line
            new_neigh_box_list.append(newbox-nb[0])         #add central box above list
        if newbox > nb[0] and newbox%nb[0] != nb[0]-1:      #test if not frst line and last column
            new_neigh_box_list.append(newbox-nb[0]+1)       #add right box list above
        #put the particle in the new_neighborlists
        for k in new_neigh_box_list :
            box[k].neighboxlist.append(self.ident)

    def changebox(self):
        if self.Mybox != -1 :
            oldbox=self.Mybox
            newbox=self.mybox()
            if newbox!=oldbox : #verify particle box change
                box[oldbox].mylist.remove(self.ident) #this is safe since particles ident is unique
                box[newbox].mylist.append(self.ident)
                self.Mybox=newbox
                #Create lists with the indexes of the neighboring boxes to change the particles there too
                #This particle was in the lists of three boxes above its box and one in the left
                old_neigh_box_list=[]
                new_neigh_box_list=[]
                if oldbox%nb[0] != 0 :                              #test if not first column
                    old_neigh_box_list.append(oldbox-1)             #index of the box on the left 
                if newbox%nb[0] != 0 :                              #test if not first column
                    new_neigh_box_list.append(newbox-1)             #add to left box list  
                if oldbox%nb[0] != 0 and oldbox > nb[0] :           #test if not first column and not first line
                    old_neigh_box_list.append(oldbox-nb[0]-1)       #add left box above list
                if newbox%nb[0] != 0 and newbox > nb[0] :           #test if not first column and not first line
                    new_neigh_box_list.append(newbox-nb[0]-1)       #add left box above list
                if oldbox  > nb[0] :                                #test if not first line
                    old_neigh_box_list.append(oldbox-nb[0])         #add central box above list   
                if newbox > nb[0] :                                 #test if not first line
                    new_neigh_box_list.append(newbox-nb[0])         #add central box above list
                if oldbox > nb[0]  and oldbox%nb[0] != nb[0]-1:     #test if not first line and last column
                    old_neigh_box_list.append(oldbox-nb[0]+1)       #add right box list above
                if newbox > nb[0] and newbox%nb[0] != nb[0]-1:      #test if not frst line and last column
                    new_neigh_box_list.append(newbox-nb[0]+1)       #add right box list above
                #take the particle off the old_neighborlists
                for k in old_neigh_box_list :
                    # print(oldbox,newbox,self.ident,k)
                    # print(box[k].neighboxlist)
                    box[k].neighboxlist.remove(self.ident)
                #put the particle in the new_neighborlists
                for k in new_neigh_box_list :
                    box[k].neighboxlist.append(self.ident)
                    # print(oldbox,newbox,self.ident,k)
                    # print(box[k].neighboxlist)
        return 

    def mov(self): #Particle moviment
        if self.Mybox != -1 :
            self.v=self.v0*self.n+self.mu*self.Force
            dr=self.v*dt
            self.r+=dr
            self.autovelchange(self.v,self.n)
            self.theta+=self.dtheta*dt/tau+self.noise*(rand.random()-0.5)*np.sqrt(dt)
            self.n[0]=np.cos(self.theta)
            self.n[1]=np.sin(self.theta)
            self.contour()
        return self.r,self.v, self.n

    
    def contour(self):
        if self.r[0]<-L[0]:
            self.n[0]=np.abs(self.n[0])
            self.v[0]=np.abs(self.v[0])
            self.theta=(rand.random()-0.5)*3.14
            self.r[0]=-L[0]+wall_osc#*rand.random()*self.Req
        if self.r[0]>L[0]: #will enter the death_list
            death_list.append(self.ident)
            live_list.remove(self.ident)
            box[self.Mybox].mylist.remove(self.ident)
            self.Mybox=-1
        if self.r[1]<-L[1]:
            self.n[1]=np.abs(self.n[1])
            self.v[1]=np.abs(self.v[1])
            self.theta=-self.theta
            self.r[1]=-L[1]+wall_osc#*self.Req*rand.random()
        if self.r[1]>L[1]:
            self.n[1]=-np.abs(self.n[1])
            self.v[1]=-np.abs(self.v[1])
            self.theta=-self.theta
            self.r[1]=L[1]-wall_osc#*self.Req*rand.random()
        normr=np.linalg.norm(self.r)
        #Particles hitting the cylinder loose their velocity component perpendicular to the cylinder if invading it
        if normr < cylinder_radius:
            vaux=self.v
            vpar=np.dot(self.v,self.r)*self.r/np.linalg.norm(self.r)**2 #velocity component parallel to cylinder radius
            vper=self.v-vpar #perpendicular component to cylinder radius (parallel to cylinder contour)
            crit = np.dot(vpar,self.r)
            if crit < 0 :
                self.v=vper
                self.r-=vaux*dt
            self.r+=self.v*dt 
            self.n=self.v/np.linalg.norm(self.v)
            self.theta = np.arctan2(self.v[1],self.v[0])
            # if self.r[1]<0:
            #     self.theta=math.atan2(self.r[0],-self.r[1])
            #     self.n = np.array([np.cos(self.theta),np.sin(self.theta)])
            #     self.v = self.v0*self.n
            # if self.r[1]>0:
            #     self.theta=math.atan2(-self.r[0],self.r[1])
            #     self.n = np.array([np.cos(self.theta),np.sin(self.theta)])
            #     self.v =  self.v0*self.n
            # self.r=self.r*cylinder_radius/normr
        return self.r, self.theta, self.n, self.v

    def autovelchange(self,v,n):
        vnorm=np.linalg.norm(v)
        u=self.v/vnorm
        self.dtheta=np.arcsin(self.n[0]*u[1]-self.n[1]*u[0])
        return self.dtheta

    def forces_between_particles(self):
        if self.Mybox != -1 :
            def force(self,dr,Req):
                normdr=np.linalg.norm(dr)
                if(normdr<Req):
                    if dr.any() < 10**(-8):
                        print(self.ident,dr)
                        exit()
                    f=self.Frep*(1/self.Req-1/normdr)*dr
                else:
                    if(normdr<self.R0):
                        f=self.Fadh*dr*(1-Req/normdr)/(self.R0-Req)
                    else:
                        f=0.0
                return f
            for i in box[self.Mybox].mylist:
                if self.ident != i :
                    dr=part[i].r-self.r
                    if np.linalg.norm(dr) < 10**(-8):
                        self.r[0]+=0.00001*rand.random()
                        dr=part[i].r-self.r
                    Req=(part[i].Req+self.Req)/2.
                    self.Force+=force(self,dr,Req)
            for i in box[self.Mybox].neighboxlist:
                Req=(part[i].Req+self.Req)/2.
                dr=part[i].r-self.r
                f=force(self,dr,Req)
                self.Force+=f
                part[i].Force-=f

    def zeroforce(self):
        self.Force=np.array([0.,0.])

    def zerocross(self): #testing the number of particles in the last cell column
        if L[0]-self.r[0] > self.R0 and L[0]-self.r[0] < 2*self.R0:
            self.cross=1
        else :
            self.cross=0

    def zerocross2(self): #testing the number of particles in the last cell column
        if L[0]+self.r[0] > self.R0 and L[0]+self.r[0] < 2*self.R0:
            self.cross2=1
        else :
            self.cross2=0

    def out_of_cylinder(self,cylinder_radius):
        if np.linalg.norm(self.r)<cylinder_radius :
            self.r[0]=cylinder_radius+(L[0]-cylinder_radius)*rand.random()
            self.r[1]=(2*rand.random()-1)*L[1]

#Box class definition
class boite:
    def __init__(self,index):
        self.index = index
        self.mylist = []
        self.neighboxlist = []

    def neighbor_list(self):
        self.neighboxlist=[]
        if self.index%nb[0]!=nb[0]-1:                                #test if not last column
            self.neighboxlist.extend(box[self.index+1].mylist)       #add box list on the right
        if self.index%nb[0]!=0 and self.index+nb[0]<nb2:             #test if not first column and not last line
            self.neighboxlist.extend(box[self.index+nb[0]-1].mylist) #add left box list below 
        if self.index+nb[0]<nb2:                                     #test if not last line
            self.neighboxlist.extend(box[self.index+nb[0]].mylist)   #add central box list below
        if self.index+nb[0]<nb2 and self.index%nb[0]!=nb[0]-1:       #test if not last line and last column
            self.neighboxlist.extend(box[self.index+nb[0]+1].mylist) #add right box list below
#Main program                                  
#global variables
global N,L,lbox,nb,dt,nb2,t,cylinder_radius,death_list,live_list
wall_osc = 0.001  #distance from wall of a reinjected boid
output_file_name="output_simu_szabo.dat"
output_counter_name = 'last_column_counter.txt'
state_file_name = "state_simu_szabo.dat"

passos=20000
save_time = 1000
input_data = sys.argv[0:7]
if len(sys.argv) !=7 :
    print("Need 6 arguments:tau,division_time,fadh,noise,size_disp, initial state\nFor example: python boidscontinous_stokes.py 4 100 0.1 1 0.2 0 ")
    exit()
tau=float(input_data[1])
division_time=int(input_data[2])
fadh=float(input_data[3])
noise=float(input_data[4])
size_disp = float(input_data[5]) #particle size Req dispersion
initial_state=int(input_data[6])  #if 0, launch initial conditions, if 1, read from file


if initial_state== 0 :
    rand.seed(0.1)
    N,L,lbox,nb,nb2,dt,exit_fig,cylinder_radius,t,figindex,death_list,live_list=initial_0()
    #initialize N particles
    part=list(particle(-L[0]+L[0]*rand.random()/4.,L[1]*2*(rand.random()-0.5), 0.1*(rand.random()-0.5),0.1*(rand.random()-0.5), i, size_disp*(rand.random()-0.5),fadh,noise) for i in range(N))
    #avoiding cells to enter de cylinder
    list(map(lambda i:i.out_of_cylinder(cylinder_radius), part)) 
    output_counter_name = 'last_column_counter.txt'
    output_counter_file = open(output_counter_name,'w')
    #save parameters, initialize outfile
    outfile_init(output_file_name,N,lbox,exit_fig,cylinder_radius,L,dt,size_disp,division_time,part[0].noise,part[0].v0,part[0].mu,part[0].Frep,part[0].Fadh,part[0].R0,tau)

if initial_state == 1:
    state_file_name="state_simu_szabo.dat"
    output_counter_name = 'last_column_counter.txt'
    output_counter_file = open(output_counter_name,'a')
    death_list, x, y, vx, vy, N, lbox, exit_fig, cylinder_radius, L, nb, nb2, dt, t, noise, Fadh, figindex, tau, size_disp, division_time = initial_1(tau,division_time,state_file_name)
    part=[]
    live_list=[]
    for i in range(N):
        part.append(particle(x[i],y[i],vx[i],vy[i],i,size_disp*(rand.random()-0.5),fadh,noise))
        if i not in death_list:
            live_list.append(i)

#create the boxes
box=list(boite(i) for i in range(nb2))

# Construct list of particles in each box
for i in range(N):
    if part[i].r[0] > L[0]:
        part[i].Mybox=-1
        # delta=part[i].r[0]-L[0]
        # part[i].r[0]-=delta*1.5
    else:
        part[i].Mybox=part[i].mybox()
        #print(i,N,part[i].Mybox,nb2, part[i].r)
        box[part[i].Mybox].mylist.append(i)

# Construct list of particles in neighboring boxes

list(map(lambda i:i.neighbor_list(), box))


#System evolution
intt=0
while(t<passos*dt):
    #Calculate the forces
    list(map(lambda i:i.forces_between_particles(), part))
    #Move all particles
    t+=dt #update time
    list(map(lambda i:i.mov(), part))
    #Find the newboxes and the neighboring boxes
    list(map(lambda i:i.changebox(), part))
    #Construct the list of particles in neighboring boxes
    #list(map(lambda i:i.neighbor_list(), box))
    #Reset forces
    list(map(lambda i:i.zeroforce(), part))

    #Count those crossing the end within a cell width 
    if intt%100 == 0 : 
        cross,cross2=0,0
        list(map(lambda i:i.zerocross(), part))
        list(map(lambda i:i.zerocross2(), part))
        for i in part:
            cross += i.cross
            cross2 += i.cross2
        output_counter_file.write("%f %d %d\n"%(t,cross,cross2))
    if intt%division_time == 0 : #division
        if death_list : #take a death to revival
            mother=rand.choice(live_list)
            while part[mother].r[0]>-3*L[0]/4:
                mother=rand.choice(live_list)
                print(mother,part[mother].r[0])
            new=death_list.pop(0)
            live_list.append(new)
            print("new_zombie=",new,part[mother].r[0])
            index = 1
            
        else:
            mother=rand.choice(live_list) #create a new cell
            while part[mother].r[0]>-3*L[0]/4:
                mother=rand.choice(live_list)
            new = N
            print(mother, "not zombie",part[mother].r[0])
            live_list.append(N)
            N+=1
            index =0
            #print(new,mother)
        #initial conditions to the new born
        x=part[mother].r[0]+0.01*(np.random.rand()-0.5)
        y=part[mother].r[1]+0.01*(np.random.rand()-0.5)
        if np.abs(x)>L[0]:x=x-np.sign(x)*0.005
        if np.abs(y)>L[1]:y=y-np.sign(y)*0.005
        vx,vy=-part[mother].v
        if index == 0:
        #create the ne w born and put it on the box
            part.append(particle(x,y,vx,vy,new,size_disp*(rand.random()-0.5),fadh,noise))
        if index == 1:
            part[new].r=np.array([x,y])
            part[new].v=np.array([vx,vy])
        newbox=part[new].mybox()
        part[new].Mybox=newbox
        if newbox > 19999:
            print(newbox,x,y)
        box[newbox].mylist.append(new)
        part[new].set_new_born_neigh()
#        if new < N : print(newbox)
    # Make a scatter graph
    delta=5.
    sizes=6.0
#    if intt > 10000 : exit_fig = 10
    if(intt%exit_fig==0):
#        print(t,cross,cross2)
        figindex=make_graph(L,live_list,part,figindex,cylinder_radius)
    if intt%save_time == 0:
        state_init(state_file_name,N,lbox,exit_fig,cylinder_radius,L,dt,t,size_disp,division_time,noise,part[0].v0,part[0].mu,part[0].Frep,part[0].Fadh,part[0].R0,tau,death_list,figindex)
        save_state(state_file_name,part)
        save_state(output_file_name,part)
    intt+=1

