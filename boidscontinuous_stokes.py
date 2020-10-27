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
#function to save system state each 10000 steps

def save_state(N,lbox,exit_fig,cylinder_radius,L,dt,t,figindex,tau,part):
    state_file_name="state_simu_szabo.dat"
    state_file = open(state_file_name,'w')
    state_file.write("Number_of_particles: %d\n"%N)
    state_file.write("Box_size: %d\n"%lbox)
    state_file.write("Steps_between_images: %d \n"%exit_fig)
    state_file.write("Radius: %f \n"%cylinder_radius)
    state_file.write("Obst_position: 0. 0.\n")
    state_file.write("Dimensions: %d %d \n"%(L[0],L[1]))
    state_file.write("Max-dist: 4 \n")
    state_file.write("dt: %f \n"%dt)
    state_file.write("t: %f \n"%t)
    state_file.write("noise: %f\n"%part[0].noise)
    state_file.write("v0: %f\n"%part[0].v0)
    state_file.write("mu: %f\n"%part[0].mu)
    state_file.write("Frep: %f\n"%part[0].Frep)
    state_file.write("Fadh: %f\n"%part[0].Fadh)
    state_file.write("R0: %f\n"%part[0].R0)
    state_file.write("figindex: %d\n"%figindex)
    state_file.write("tau: %f\n"%tau)
    state_file.write("size_disp: %f\n"%size_disp)           
    state_file.write("division_time: %d \n"%division_time)
    state_file.write("death_list: ")
    for i in death_list:
        state_file.write("%d "%i)
    state_file.write("\n")
    state_file.write("x y vx vy \n")
    x,y,vx,vy=[],[],[],[]
    map(lambda i:x.append(i.r[0]), part)
    map(lambda i:y.append(i.r[1]), part)
    map(lambda i:vx.append(i.v[0]), part)
    map(lambda i:vy.append(i.v[1]), part)
    for i in range(len(x)):
        state_file.write("%d %f %f %f %f \n"%(i,x[i],y[i],vx[i],vy[i]))


#Particle class definition
class particle:
   # noise=0.6 #original value
    noise=0.1
#    noise_T=10.0
    v0 = 0.1
    mu=1.0
    Frep=30.0
    #Fadh=0.75 #original value
#    Fadh=10.0
    Fadh=10.0
    #Req=5./6. # original work by Szabo
    #R0=1.0
#    Req=1.0
    R0=3.0
    def __init__(self, x, y, vx, vy, ident, Raio_equilibrio):
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
    
    def mybox(self): #Each particle calculates the box it is in
        if np.isnan(self.r[0]) == True or np.isnan(self.r[1]) == True :
            print self.ident, self.r, self.Force
            exit()
        j=int((self.r[0]+L[0])/lbox)+nb[0]*int((self.r[1]+L[1])/lbox)
        return j

    def changebox(self):
        if self.Mybox != -1 :
            newbox=self.mybox()
            if newbox!=self.Mybox : #verify particle box change
                box[self.Mybox].mylist.remove(self.ident) #this is safe since particles ident is unique
                box[newbox].mylist.append(self.ident)
#                if self.ident not in box[self.Mybox].mylist:
#                    print box[self.Mybox].mylist,self.ident,self.r, newbox, self.Mybox,t
#                    exit()
                #box[self.Mybox].mylist.remove(self.ident)
                self.Mybox=newbox
        return self.Mybox

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
                        print self.ident,dr
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
initial_state=0 #if 0 launch initial conditions if 1 read from file
wall_osc = 0.001  #distance from wall of a reinjected boid
size_disp = 0.1  #particle size Req dispersion
output_file_name="output_simu_szabo.dat"
output_counter_name = 'last_column_counter.txt'
passos=1800000
save_time = 1000
input_data = sys.argv[0:4]
tau=float(input_data[1])
division_time=int(input_data[2])
initial_state=int(input_data[3])    
if initial_state== 0 :
    N=10
    L=np.array([100,50])
    lbox=1
    nb=2*L#/lbox)
    nb2=nb[1]*nb[0]
    dt=0.01
    exit_fig=1000
    rand.seed(0.1)
    cylinder_radius=15
    t=0
    figindex=0
    death_list=[]
    live_list=list(range(N))
    output_file = open(output_file_name,'w')
    output_counter_file = open(output_counter_name,'w')
    output_file.write("Number_of_particles: %d\n"%N)
    output_file.write("Box_size: %d\n"%lbox)
    output_file.write("Steps_between_images: %d \n"%exit_fig)
    output_file.write("Radius: %f \n"%cylinder_radius)
    output_file.write("Obst_position: 0. 0.\n")
    output_file.write("Dimensions: %d %d \n"%(L[0],L[1]))
    output_file.write("Max-dist: 4 \n")
    output_file.write("dt: %f \n"%dt)
    output_file.write("size-dispersion: %f \n"%size_disp)
    output_file.write("division_time: %f \n"%division_time)

#initialize N particles
    #part=list(particle(-L[0]+L[0]/2*rand.random(),L[1]*2*(rand.random()-0.5), rand.random()-0.5,rand.random()-0.5, i) for i in range(N))
    #part=list(particle(-L[0]+3*L[0]/4.*rand.random(),L[1]*2*(rand.random()-0.5), rand.random()-0.5,rand.random()-0.5, i, size_disp*(rand.random()-0.5) ) for i in range(N))
    #part=list(particle(L[0]*2*(rand.random()-0.5),L[1]*2*(rand.random()-0.5), rand.random()-0.5,rand.random()-0.5, i, size_disp*(rand.random()-0.5) ) for i in range(N))
    #part=list(particle(-L[0]+3*L[0]/4.*rand.random(),L[1]*2*(rand.random()-0.5), 1.,0., i, size_disp*(rand.random()-0.5) ) for i in range(N))
    part=list(particle(-L[0]+L[0]*rand.random()/4.,L[1]*2*(rand.random()-0.5), 1.,0., i, size_disp*(rand.random()-0.5) ) for i in range(N))
    map(lambda i:i.out_of_cylinder(cylinder_radius), part) #avoiding cells to enter de cylinder
    output_file.write("noise: %f\n"%part[0].noise)
    output_file.write("v0: %f\n"%part[0].v0)
    output_file.write("mu: %f\n"%part[0].mu)
    output_file.write("Frep: %f\n"%part[0].Frep)
    output_file.write("Fadh: %f\n"%part[0].Fadh)
    output_file.write("R0: %f\n"%part[0].R0)
    output_file.write("tau: %f\n"%tau)                  

if initial_state == 1:
    state_file_name = "state_simu_szabo.dat"
    state_file = open(state_file_name)
    output_file=open(output_file_name,"a")
    output_counter_file = open(output_counter_name,'a')
    line_splitted = [1]
    part=[]
    live_list=[]
    part_counter = 0
    
    while line_splitted[0] != "x":
        line = state_file.readline()
        if not line :
            break
        line_splitted = line.split()
        if line_splitted[0] == "Number_of_particles:" :
            N = int(line_splitted[1])
        if line_splitted[0] == "Box_size:" :
            lbox = int(line_splitted[1])
        if line_splitted[0] ==  "Steps_between_images:":
            exit_fig = int(line_splitted[1])
        if line_splitted[0] ==  "Radius:":
            cylinder_radius = float(line_splitted[1])
        if line_splitted[0] == "Dimensions:":
            L0,L1           = int(line_splitted[1]),int(line_splitted[2])
            L=np.array([L0,L1])
            nb=2*L#/lbox)
            nb2=nb[1]*nb[0]
        if line_splitted[0] ==  "dt:":
            dt = float(line_splitted[1])
        if line_splitted[0] ==  "t:":
            t = float(line_splitted[1])
        if line_splitted[0] == "noise:":
            noise = float(line_splitted[1])
        if line_splitted[0] == "v0:":
            v0 = float(line_splitted[1])
        if line_splitted[0] == "mu:":
            mu = float(line_splitted[1])
        if line_splitted[0] == "Frep:":
            Frep            = float(line_splitted[1])
        if line_splitted[0] == "Fadh:":
            Fadh = line_splitted[1]
        if line_splitted[0] == "R0:":
            R0              = float(line_splitted[1])
        if line_splitted[0] == "figindex:":
            figindex        = int(line_splitted[1])
        if line_splitted[0] == "tau:":
            taul             = float(line_splitted[1])
            if taul != tau :
                print "Attention! Previous simu performed with other tau"
                exit()
        if line_splitted[0] == "size_disp:":
            size_disp             = float(line_splitted[1])
        if line_splitted[0] == "division_time:":
            division_timel     = int(line_splitted[1])
            if division_timel != division_time :
                print "Attention! Previous simu performed with other division_time"
                exit()
        if line_splitted[0] == "death_list:":
            if len(line_splitted) == 1 :
                death_list = []
            else:
                death_list=map(int,line_splitted[1:])
        
            #            state_file.write("division_time: %f \n"%division_time)
    while 1 :
        line = state_file.readline()
        if not line :
            break
        line_splitted = line.split()
        live_list.append(int(line_splitted[0]))
        part.append(particle(float(line_splitted[1]),float(line_splitted[2]),float(line_splitted[3]),float(line_splitted[4]),int(line_splitted[0]),size_disp*(rand.random()-0.5)))
box=list(boite(i) for i in range(nb2))

# Construct list of particles in each box
#print len(part),N

for i in range(N):
    if part[i].r[0] > L[0]:
        delta=part[i].r[0]-L[0]
        part[i].r[0]-=delta*1.5
    part[i].Mybox=part[i].mybox()
    #print i,N,part[i].Mybox,nb2, part[i].r
    box[part[i].Mybox].mylist.append(i)
    # Construct list of particles in neighboring boxes

map(lambda i:i.neighbor_list(), box)

#System evolution
intt=0
while(t<passos*dt):
    #Calculate the forces
    map(lambda i:i.forces_between_particles(), part)
    #Move all particles
    t+=dt #update time
    map(lambda i:i.mov(), part)
    #Find the newboxes
    map(lambda i:i.changebox(), part)
    #Construct the list of particles in neighboring boxes
    map(lambda i:i.neighbor_list(), box)
    #Reset forces
    map(lambda i:i.zeroforce(), part)
    #Count those crossing the end within a cell width 
    if intt%100 == 0 : 
        cross,cross2=0,0
        map(lambda i:i.zerocross(), part)
        map(lambda i:i.zerocross2(), part)
        for i in part:
            cross += i.cross
            cross2 += i.cross2
        output_counter_file.write("%f %d %d\n"%(t,cross,cross2))
    if intt%division_time == 0 : #division
        if death_list : #take a death to revival
            mother=rand.choice(live_list)
            while part[mother].r[0]>-3*L[0]/4:
                mother=rand.choice(live_list)
                print mother,part[mother].r[0]
            new=death_list.pop(0)
            live_list.append(new)
            print "new_zombie=",new,part[mother].r[0]
            index = 1
            
        else:
            mother=rand.choice(live_list) #create a new cell
            while part[mother].r[0]>-3*L[0]/4:
                mother=rand.choice(live_list)
            new = N
            print mother, "not zombie",part[mother].r[0]
            live_list.append(N)
            N+=1
            index =0
            #print new,mother
        #initial conditions to the new born
        x=part[mother].r[0]+0.01*(np.random.rand()-0.5)
        y=part[mother].r[1]+0.01*(np.random.rand()-0.5)
        if np.abs(x)>L[0]:x=x-np.sign(x)*0.005
        if np.abs(y)>L[1]:y=y-np.sign(y)*0.005
        vx,vy=-part[mother].v
        if index == 0:
        #create the ne w born and put it on the box
            part.append(particle(x,y,vx,vy,new,size_disp*(rand.random()-0.5)))
        if index == 1:
            part[new].r=np.array([x,y])
            part[new].v=np.array([vx,vy])
        newbox=part[new].mybox()
        part[new].Mybox=newbox
        if newbox > 19999:
            print newbox,x,y
        box[newbox].mylist.append(new)
#        if new < N : print newbox
    # Make a scatter graph
    delta=5.
    sizes=6.0
#    if intt > 10000 : exit_fig = 10
    if(intt%exit_fig==0):
#        print t,cross,cross2
        output_file.write("x y vx vy \n")
        plt.axis([-L[0]-delta,L[0]+delta,-L[1]-delta,L[1]+delta])
        plt.axes().set_aspect(1.0)
        circle=plt.Circle((0.,0.),radius=cylinder_radius-0.5,color='r')
        x,y,vx,vy=[],[],[],[]
        for i in live_list:
            x.append(part[i].r[0])
            y.append(part[i].r[1])
            vx.append(part[i].v[0])
            vy.append(part[i].v[1])
        for i in range(len(x)):
            output_file.write("%f %f %f %f \n"%(x[i],y[i],vx[i],vy[i]))

        plt.scatter(x,y,s=sizes,alpha=0.3)
        name=str(figindex)+".png"
        fig = plt.gcf()
        plt.rc("savefig",dpi=300)
        ax = fig.gca()
        ax.add_artist(circle)
        fig.savefig(name,bbox_inches='tight')
        figindex+=1
        fig.clf()
    if intt%save_time == 0:
        save_state(N,lbox,exit_fig,cylinder_radius,L,dt,t,figindex,tau,part)
    intt+=1

