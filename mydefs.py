def outfile_init(output_file_name,N,lbox,exit_fig,cylinder_radius,L,dt,size_disp,division_time,noise,v0,mu,Frep,Fadh,R0,tau):
    output_file = open(output_file_name,'w')
    output_file.write("Number_of_particles: %d\n"%N)
    output_file.write("Box_size: %d\n"%lbox)
    output_file.write("Steps_between_images: %d \n"%exit_fig)
    output_file.write("Radius: %f \n"%cylinder_radius)
    output_file.write("Obst_position: 0. 0.\n")
    output_file.write("Dimensions: %d %d \n"%(L[0],L[1]))
    output_file.write("Max-dist: 4 \n")
    output_file.write("dt: %f \n"%dt)
    output_file.write("size-dispersion: %f \n"%size_disp)
    output_file.write("division_time: %d \n"%division_time)

    output_file.write("noise: %f\n"%noise)
    output_file.write("v0: %f\n"%v0)
    output_file.write("mu: %f\n"%mu)
    output_file.write("Frep: %f\n"%Frep)
    output_file.write("Fadh: %f\n"%Fadh)
    output_file.write("R0: %f\n"%R0)
    output_file.write("tau: %f\n"%tau)
    return

def state_init(output_file_name,N,lbox,exit_fig,cylinder_radius,L,dt,t,size_disp,division_time,noise,v0,mu,Frep,Fadh,R0,tau,death_list,figindex):
    output_file = open(output_file_name,'w')
    output_file.write("Number_of_particles: %d\n"%N)
    output_file.write("Box_size: %d\n"%lbox)
    output_file.write("Steps_between_images: %d \n"%exit_fig)
    output_file.write("Radius: %f \n"%cylinder_radius)
    output_file.write("Obst_position: 0. 0.\n")
    output_file.write("Dimensions: %d %d \n"%(L[0],L[1]))
    output_file.write("Max-dist: 4 \n")
    output_file.write("dt: %f \n"%dt)
    output_file.write("t: %f \n"%t)
    output_file.write("size-disp: %f \n"%size_disp)
    output_file.write("noise: %f\n"%noise)
    output_file.write("v0: %f\n"%v0)
    output_file.write("mu: %f\n"%mu)
    output_file.write("Frep: %f\n"%Frep)
    output_file.write("Fadh: %f\n"%Fadh)
    output_file.write("R0: %f\n"%R0)
    output_file.write("tau: %f\n"%tau)
    output_file.write("figindex: %d\n"%figindex)
    output_file.write("death_list: ")
    for i in death_list:
        output_file.write("%d "%i)
    output_file.write("\n")
    return

def save_state(output_file_name,part):
    output_file = open(output_file_name,'a')
    output_file.write("x y vx vy \n")
    for i in part:
        output_file.write("%f %f %f %f \n"%(i.r[0],i.r[1],i.v[0],i.v[1]))

def make_graph(L,live_list,part,figindex,cylinder_radius):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    delta=5.
    sizes=6.0
    fig=plt.figure()
    plt.axis([-L[0]-delta,L[0]+delta,-L[1]-delta,L[1]+delta])
    plt.axes().set_aspect(1.0)
    fig.patch.set_facecolor('#E0E0E0')
    fig.patch.set_alpha(0.7)
    circle=plt.Circle((0.,0.),radius=cylinder_radius-0.5,color='r')
    x,y,vx,vy=[],[],[],[]
    for i in live_list:
        x.append(part[i].r[0])
        y.append(part[i].r[1])
        vx.append(part[i].v[0])
        vy.append(part[i].v[1])

    plt.scatter(x,y,s=sizes,alpha=0.5)
    name=str(figindex)+".png"
    fig = plt.gcf()
    plt.rc("savefig",dpi=300)
    ax = fig.gca()
    ax.patch.set_facecolor('#E0E0E0')
    ax.patch.set_alpha(0.5)
    ax.add_artist(circle)
    fig.savefig(name,bbox_inches='tight',facecolor=fig.get_facecolor())
    figindex+=1
    fig.clf()
    return figindex

def initial_0():
    import numpy as np
    N=10
    L=np.array([100,50])
    lbox=2
    nbf=2*L/lbox
    nb=nbf.astype(int)
    nb2=nb[1]*nb[0]
    dt=0.01
    exit_fig=1000
    cylinder_radius=15
    t=0
    figindex=0
    death_list=[]
    live_list=list(range(N))
    return N,L,lbox,nb,nb2,dt,exit_fig,cylinder_radius,t,figindex,death_list,live_list

def initial_1(tau,division_time,v0,state_file_name):
    import numpy as np
    import random as rand
    state_file = open(state_file_name)
    line_splitted = [1]
    
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
            nbf=2*L/lbox
            nb=nbf.astype(int)
            nb2=nb[1]*nb[0]
        if line_splitted[0] ==  "dt:":
            dt = float(line_splitted[1])
        if line_splitted[0] ==  "t:":
            t = float(line_splitted[1])
        if line_splitted[0] == "noise:":
            noise = float(line_splitted[1])
        if line_splitted[0] == "v0:":
            v0l = float(line_splitted[1])
            print(v0l,v0)
            delta_v0=np.abs(v0l-v0)
            if delta_v0 > 1.e-6 :
                print("Attention! Previous simu performed with other v0\n Exiting!!")
                exit()

        if line_splitted[0] == "mu:":
            mu = float(line_splitted[1])
        if line_splitted[0] == "Frep:":
            frep            = float(line_splitted[1])
        if line_splitted[0] == "R0:":
            R0              = float(line_splitted[1])
        if line_splitted[0] == "Fadh:":
            fadh = line_splitted[1]
        if line_splitted[0] == "figindex:":
            figindex        = int(line_splitted[1])
        if line_splitted[0] == "tau:":
            taul             = float(line_splitted[1])
            delta_tau=np.abs(taul-tau)
            if delta_tau > 1.e-6 :
                print("Attention! Previous simu performed with other tau!\n Exiting!!")
                exit()
        if line_splitted[0] == "size-disp:":
            size_disp             = float(line_splitted[1])
        if line_splitted[0] == "division_time:":
            division_timel     = int(line_splitted[1])
            delta_div=np.abs(division_timel-division_time)
            if delta_div > 0 :
                print("Attention! Previous simu performed with other division_time")
                exit()
        if line_splitted[0] == "death_list:":
            if len(line_splitted) == 1 :
                death_list = []
            else:
                death_list=list(map(int,line_splitted[1:]))
        
            #            state_file.write("division_time: %f \n"%division_time)
    live_list,x,y,vx,vy=[],[],[],[],[]
    while 1 :
        line = state_file.readline()
        if not line :
            break
        line_splitted = line.split()
        x.append(float(line_splitted[0]))
        y.append(float(line_splitted[1]))
        vx.append(float(line_splitted[2]))
        vy.append(float(line_splitted[3]))


    return death_list, x, y, vx, vy, N, lbox, exit_fig, cylinder_radius, L, nb, nb2, dt, t, v0, fadh, figindex, tau, size_disp, division_time
