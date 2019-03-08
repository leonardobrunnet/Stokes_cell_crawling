import copy
import math as math
import os
import sys
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt

#************************************************************
# diagonalization of texture matrices. Input: texture_box
# output: axis_a, axis_b and angle in degrees
def axis_angle(M):
    w,v        = np.linalg.eig(M)
    index_max  = np.argmax(w)
    index_min  = np.argmin(w)
    axis_a     = w[index_max]
    axis_b     = w[index_min]
    ang_elipse = math.atan2(v[1,index_max],v[0,index_max])*180/math.pi
    return axis_a,axis_b,ang_elipse
            
# function mapping delaunay out index to particle index

def map_focus_region_to_part(points,list_neighbors,index_particle):
    for i,w in enumerate(points):
        aux=index_particle[i]
        part[aux].r=np.array(w)
        part[aux].list_neigh=[]
        for j in list_neighbors[i]:
            part[aux].list_neigh.append(index_particle[j])

############Particle class definition##############
class particle:

    def __init__(self,ident):
        self.r=np.array([0,0])
        self.r_old=np.array([0,0])
        self.ident=ident  #indice geral da particula
        self.list_neigh=[]
        self.list_neigh_old=[]
        self.M=np.zeros((2,2))
        self.m_list=[]
        self.dm_list=[]
        self.lc_list=[] 
        self.la_list=[]
        self.ld_list=[]
        self.C=np.zeros((2,2))
        self.CT=np.zeros((2,2))
        self.B=np.zeros((2,2))
        self.T=np.zeros((2,2))

    def texture(self):
        self.M=np.zeros((2,2))
        n=len(self.list_neigh)
        if n > 0:
            for i in self.list_neigh:
                self.M+=self.mat(i)
            self.M/=n


        
    def mat(self,i):
        l=self.r-part[i].r
        m=np.outer(l,l)
        return m
    
    def copy_to_old(self):
        self.list_neigh_old=copy.deepcopy(self.list_neigh)
        self.r_old=copy.deepcopy(self.r)

    def l_av_dl(self,i):
        l=self.r-part[i].r
        l_old=self.r_old-part[i].r_old
        lav=(l+l_old)/2.
        dl=l-l_old
        return lav,dl
        
    def UT(self):
        list_c=list(set(self.list_neigh).intersection(self.list_neigh_old))
        list_a=list(set(self.list_neigh).difference(self.list_neigh_old))
        list_d=list(set(self.list_neigh_old).difference(self.list_neigh))
        self.zeros()
        Nc=len(list_c)
        Na=len(list_a)
        Nd=len(list_d)
        Ntot=Nc+Na+Nd
        for i in list_c:
            lav,dl=self.l_av_dl(i)
            c=np.outer(lav,dl)
            ct=np.outer(dl,lav)
            self.C+=c
            self.CT+=ct
        if Nc > 0 :
            self.C/=Ntot
            self.CT/=Ntot
            self.B=self.C+self.CT
        for i in list_a:
            ma=self.mat(i)
            self.T+=ma
        for i in list_d:
            md=self.mat(i)
            self.T-=md
        self.T/=Ntot
            

              
    def zeros(self):
        self.average_m=np.zeros((2,2))
        self.m_list=[]
        self.dm_list=[]
        self.C=np.zeros((2,2))
        self.CT=np.zeros((2,2))
        self.B=np.zeros((2,2))
        self.T=np.zeros((2,2))
###############Particle class definition ends here###########



def delaunay(points):
    tri = Delaunay(points)
    x,y=[],[]
    z,zz=[],[]
    for i,w in enumerate(points):
        if i%50 == 0:
            x.append(w[0])
            y.append(w[1])
        else :
            z.append(w[0])
            zz.append(w[1])
    fig=plt.scatter(x,y,s=30,c='b')
    fig=plt.scatter(z,zz,s=30,c='g')
    plt.show()
    list_neigh = [ [] for i in range(len(points)) ]
    for i in tri.simplices:
        for j in i:
            for l in i:
                if l != j :
                    if np.linalg.norm(points[j]-points[l]) < 10 : #
                        if l not in list_neigh[j]:
                            list_neigh[j].append(l)
                        if j not in list_neigh[l]:
                            list_neigh[l].append(j)
    x,y=[],[]
    for i,w in enumerate(list_neigh) :
        if i%50==0 :
            for j in w :
                x.append(points[j][0])
                y.append(points[j][1])
    fig=plt.scatter(x,y,s=30,c='r')
    fig=plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
    plt.show()
    exit()

    return list_neigh

def create_gnu_script(arrow_size, box_per_line_x, box_per_column_y, vel_win_file_name, dens_win_file_name, path):
    proportion_x, proportion_y             = 1.0, 0.7
    grid_x, grid_y, levels                 = 200, 200, 4
    image_resolution_x, image_resolution_y = 1300, 1300
    name_output_map            = "densidade-velocidade.png"
    file_script_den_vel        = open(path+"/scriptdenvel.gnu","w")
   
    file_script_den_vel.write("set size %1.2f,%1.2f \n"% (proportion_x, proportion_y))
    file_script_den_vel.write("set palette defined ( 0 '#000000',\\\n")
    file_script_den_vel.write("                      1 '#0000ff',\\\n")
    file_script_den_vel.write("                      2 '#00ffff',\\\n")
    file_script_den_vel.write("                      3 '#00ff00',\\\n")
    file_script_den_vel.write("                      4 '#ffff00',\\\n")
    file_script_den_vel.write("                      5 '#ff0000')\n")
    file_script_den_vel.write("set nokey \n")
    file_script_den_vel.write("set dgrid3d %d,%d,%d \n"% (grid_x, grid_y, levels))
    file_script_den_vel.write("set pm3d explicit \n")
    file_script_den_vel.write("set output \"| head -n -2 > toto.dat\" \n")
    file_script_den_vel.write("set table \n")
    file_script_den_vel.write("splot \"%s\" using 1:2:3 \n"% dens_win_file_name)
    file_script_den_vel.write("unset table \n")
    file_script_den_vel.write("unset dgrid3d \n")
    file_script_den_vel.write("mtf = %f \n"% arrow_size)
    file_script_den_vel.write(" \n")
    file_script_den_vel.write(" \n")
    file_script_den_vel.write("set pm3d map \n")
    file_script_den_vel.write("splot [%d:%d][%d:%d] \"toto.dat\" \n"% (0, box_per_line_x, 0, box_per_column_y))
    file_script_den_vel.write("replot \"%s\" u($1):($2):(0.0):(mtf*$3):(mtf*$4):(0.0) with vectors head size 1.5,20,60 lt rgb \"black\" \n"% vel_win_file_name)
    file_script_den_vel.write("pause -1 \n")
    file_script_den_vel.write("set terminal png large size %d,%d \n"% (image_resolution_x, image_resolution_y)) 
    file_script_den_vel.write("set output \"%s\" \n"% name_output_map)
    file_script_den_vel.write("replot \n")  


def read_param_vic_greg(file_par_simu) :
    while 1 :
        line = file_par_simu.readline()
        if not line:
            break #EOF
        if line.replace( '\r', '' ) == '\n' : #identifies blank lines
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = file_par_simu.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        if line_splitted[1] == 'Lx' :
            Lx = int(line_splitted[2])
        if line_splitted[1] == 'Ly' :
            Ly = int(line_splitted[2])
        if line_splitted[1] == 'R_ESFERA' :
            R_OBST = float(line_splitted[2])
        if line_splitted[1] == 'L_CENTRO_X' :
            X_OBST = float(line_splitted[2])
        if line_splitted[1] == 'L_CENTRO_Y' :
            Y_OBST = float(line_splitted[2])
        if line_splitted[1] == 'R_MAX' :
            box_size = float(line_splitted[2])
    return Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size

def read_param(file_input_parameter) :
    while 1 :
        line = file_input_parameter.readline()
        if not line:
            break #EOF
        if line.replace( '\r', '' ) == '\n' : #identifies blank lines
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = file_input_parameter.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        if line_splitted[0] == 'window' :
            window_size = int(line_splitted[1])
        if line_splitted[0] == 'time_0' :
            time_0 = int(line_splitted[1])
        if line_splitted[0] == 'time_f' :
            time_f = int(line_splitted[1])
        if line_splitted[0] == 'voronoi' :
            voronoi = line_splitted[1]
        if line_splitted[0] == 'obstacle' :
            obstacle = line_splitted[1]
        if line_splitted[0] == 'x0' :
            x0 = float(line_splitted[1])
        if line_splitted[0] == 'xf' :
            xf = float(line_splitted[1])
        if line_splitted[0] == 'y0' :
            y0 = float(line_splitted[1])
        if line_splitted[0] == 'yf' :
            yf = float(line_splitted[1])
        if line_splitted[0] == 'file' :
            filename = line_splitted[1]
    return window_size, time_0, time_f, obstacle, voronoi, x0, xf, y0, yf, filename



def box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf) :
    box_total   = box_per_column_y * box_per_line_x
    vx_now         = list(0. for i in range(box_total))
    vy_now         = list(0. for i in range(box_total))
    density_now    = list(0  for i in range(box_total))
    vx_tot         = list(0. for i in range(box_total))
    vy_tot         = list(0. for i in range(box_total))
    density_tot    = list(0  for i in range(box_total))
    vx_win         = list(0. for i in range(box_total))
    vy_win         = list(0. for i in range(box_total))
    density_win    = list(0  for i in range(box_total))
    axis_a_tot     = list(0. for i in range(box_total))
    axis_b_tot     = list(0. for i in range(box_total))
    ang_elipse_tot = list(0. for i in range(box_total))
    axis_a_win     = list(0. for i in range(box_total))
    axis_b_win     = list(0. for i in range(box_total))
    ang_elipse_win = list(0. for i in range(box_total))
    #    ratio=float(box_per_column_y)/box_per_line_x
    ratio = float((yf-y0)) / (xf-x0)
    vid_def.write("set size ratio %f  \n" % ratio)
    vid_def.write("set size ratio %f  \n" % ratio)
    vid_def.write("set xrange [0:%f]  \n" % box_per_line_x)
    vid_def.write("set yrange [0:%f]  \n" % box_per_column_y)    
    vel_win.write("set size ratio %f  \n" % ratio)
    vel_win.write("arrow=2.5\n")
    vid_veloc_dens.write("set size ratio %f  \n" % ratio)
    vid_veloc_dens.write("arrow=1.\n")
    vid_veloc_dens.write("unset key \n")
    vid_veloc_dens.write("set cbrange [0:1] \n")
    vid_veloc_dens.write("set palette defined ( 0 '#0000ff',\\\n")
    vid_veloc_dens.write("                      1 '#00ffff',\\\n")
    vid_veloc_dens.write("                      2 '#00ff00',\\\n")
    vid_veloc_dens.write("                      3 '#ffff00',\\\n")
    vid_veloc_dens.write("                      4 '#ff0000')\n")

    return box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, density_tot, vx_win, vy_win, \
           density_win, axis_a_tot, axis_b_tot, ang_elipse_tot, axis_a_win, axis_b_win, ang_elipse_win

def box_variables_definition_experiment(box_per_column_y, box_per_line_x):
    box_total   = box_per_column_y * box_per_line_x
    ratio          = float(box_per_column_y) / box_per_line_x
    vx_tot         = list(0. for i in range(box_total))
    vy_tot         = list(0. for i in range(box_total))
    density_tot    = list(0  for i in range(box_total))
    axis_a_tot     = list(0. for i in range(box_total))
    axis_b_tot     = list(0. for i in range(box_total))
    ang_elipse_tot = list(0. for i in range(box_total))
    vid_def.write("set size ratio %f  \n" % ratio)
    vid_veloc_dens.write("set size ratio %f  \n" % ratio)
    vid_veloc_dens.write("arrow=1.\n")
    vid_veloc_dens.write("unset key \n")
    vid_veloc_dens.write("set cbrange [0:1] \n")
    vid_veloc_dens.write("set palette defined ( 0 '#0000ff',\\\n")
    vid_veloc_dens.write("                      1 '#00ffff',\\\n")
    vid_veloc_dens.write("                      2 '#00ff00',\\\n")
    vid_veloc_dens.write("                      3 '#ffff00',\\\n")
    vid_veloc_dens.write("                      4 '#ff0000')\n")
    vel_win.write("set size ratio %f  \n" % ratio)
    vel_win.write("arrow=1.\n")

    return box_total, ratio, vx_tot, vy_tot, density_tot, axis_a_tot, axis_b_tot, ang_elipse_tot


def velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0, x0, y0, xf, yf):
    #Here we write each image to the gnuplot velocity-density movie script
    vid_veloc_dens.write("plot [%f:%f] [%f:%f] \'-\' u ($1):($2):(arrow*$3):(arrow*$4):($5) with vectors head size  0.6,20,60  filled palette title \"%d\"\n" % \
    (0, xf-x0, 0, yf-y0, image))
    if system_type == "experiment":
        v0      = 1 #this should be the real velocity, if we can measure...
        density = 1 #this should be changed if we measure real density (density_now)
        for i in range(box_total):
            module = math.sqrt(vx_now[i]**2+vy_now[i]**2)
            if module > 0. :
                vid_veloc_dens.write('%i %i %f %f %f %f\n' % (x[i], y[i], vx_now[i]/module, vy_now[i]/module, module, density)) #density_now should be used case we have it
    else :
        for box in range(box_total) :
            if density_now[box] != 0 :
                x, y = box % box_per_line_x, box / box_per_line_x
                module = math.sqrt(vx_now[box]**2 + vy_now[box]**2)
                if module > 0. :
                    vid_veloc_dens.write('%i %i  %f %f %f %f %d\n' % (x, y, vx_now[box]/module, vy_now[box]/module, module/v0, density_now[box], box))
    vid_veloc_dens.write("e \n")
    vid_veloc_dens.write("pause .1 \n")

def deformation_elipsis_script(x, y, axis_b, axis_a, ang_elipse, system_type) :
    #Deformation elipsis gnuplot script
    if system_type == 'experiment' :
        for i in range(box_total) :
            vid_def.write("set object %i ellipse at %i,%i size %f,0.5 angle %f \n" % (i+1, x[i], y[i], axis_b[i] / (2*axis_a[i]), ang_elipse[i]))
        vid_def.write("plot \'-\' w d notitle\n")
        for i in range(box_total) :
            vid_def.write("%i %i \n" % (x[i], y[i]))
        vid_def.write("e \n")
        vid_def.write("pause .1 \n")
        vid_def.write("unset for [i=1:%i] object i \n" % (box_total+1))

        
def deformation_elipsis_script_simu(box_per_line_x, box_total, axis_b, axis_a, ang_elipse) :
    #Deformation elipsis gnuplot script for simus
    for i in range(box_total) :
        if axis_a[i] != 0.0 : 
            x = i % box_per_line_x
            y = i / box_per_line_x
            vid_def.write("set object %i ellipse at %i,%i size %f,1.0 angle %f \n" % (i+1, x, y, axis_b[i] / (axis_a[i]), ang_elipse[i]))
    vid_def.write("plot \'-\' w d notitle\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            x = i % box_per_line_x
            y = i / box_per_line_x
            vid_def.write("%i %i \n" % (x, y))
    vid_def.write("e \n")
    vid_def.write("pause .1 \n")
    vid_def.write("unset for [i=1:%i] object i \n" % (box_total+1))
        
def  zero_borders_and_obstacle(box_per_line_x, box_per_column_y, r_obst, x_obst, y_obst, density_tot, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type) :
    center_x = box_per_line_x/2
    center_y = box_per_column_y/2
    for i in range(box_total):
        if system_type == 'experiment':
            bx = int(i/box_per_column_y)
            by = i%box_per_column_y
        else:
            bx = i%box_per_line_x
            by = int(i/box_per_line_x)        
        if bx == 0 or bx == box_per_line_x-1 or by == 0 or by == box_per_column_y-1 :
            density_tot[i]    = -10
            vx_tot[i]         = 0.
            vy_tot[i]         = 0.
            axis_a_tot[i]     = 0.
            axis_b_tot[i]     = 0.
            ang_elipse_tot[i] = 0.
        if system_type == 'experiment' :
            if math.sqrt((bx-x_obst)**2 + (by-y_obst)**2) <= r_obst :
                density_tot[i]    = -10
                vx_tot[i]         = 0.
                vy_tot[i]         = 0.
                axis_a_tot[i]     = 0.
                axis_b_tot[i]     = 0.
                ang_elipse_tot[i] = 0.
        if system_type == 'superboids' :
            if math.sqrt((bx-center_x)**2 + (by-center_y)**2) < r_obst :
                density_tot[i]    = -10
                vx_tot[i]         = 0.
                vy_tot[i]         = 0.
                axis_a_tot[i]     = 0.
                axis_b_tot[i]     = 0.
                ang_elipse_tot[i] = 0.
        if system_type == 'vicsek-gregoire' :
            if math.sqrt((bx-x_obst)**2 + (by-y_obst)**2) <= r_obst:
                density_tot[i]    = -10
                vx_tot[i]         = 0.
                vy_tot[i]         = 0.
                axis_a_tot[i]     = 0.
                axis_b_tot[i]     = 0.
                ang_elipse_tot[i] = 0.

    return box_per_line_x, box_per_column_y, density_tot, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot

def average_density_velocity_deformation_experiment(box_per_line_x, box_per_column_y, x, y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, image_counter):
    arrow = 1.5
    box_total = box_per_line_x*box_per_column_y
    for i in range(box_total):
        vx_tot[i] /= image_counter
        vy_tot[i] /= image_counter

    vel_win.write("plot [%f:%f] [%f:%f] \'-\' u ($1):($2):(%f*$3):(%f*$4):($5)  with vectors notitle head size  0.3,20,60  filled palette \n" % (0, box_per_line_x, 0, box_per_column_y, arrow, arrow))

    for i in range(box_total):
        dens_win.write("%d %d %f \n" % (x[i],y[i],density_tot[i]))
        module=math.sqrt(vx_tot[i]**2 + vy_tot[i]**2)
        if module >0 :
            vel_win.write("%d %d %f %f %f \n" % (x[i], y[i], vx_tot[i]/module, vy_tot[i]/module, module))
            def_win.write("%d %d %f %f %f \n" % (x[i], y[i], axis_a_tot[i],    axis_b_tot[i],    ang_elipse_tot[i]))

    vel_win.write("e \n")
    vel_win.write("pause -1 \n")

def average_density_velocity_deformation(box_per_line_x, box_per_column_y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, density_tot, vx_win, vy_win, axis_a_win, axis_b_win, density_win, count_events, v0, vel_win_file_name, dens_win_file_name, path, image_counter) :

    box_total         = box_per_column_y*box_per_line_x
    window_size_h     = window_size/2

    if system_type != 'experiment':
        count_box_win = list(0 for i in range(box_total))
        for bx in range(window_size_h+1, box_per_line_x-window_size_h):
            for by in range(window_size_h+1, box_per_column_y-window_size_h):
                for k in range(-window_size_h, window_size_h):
                    for l in range(-window_size_h, window_size_h):
                        if density_tot[(bx+k)+(by+l)*box_per_line_x] > 0 :
                            box = bx + (by*box_per_line_x)
                            density_win[box] += density_tot[(bx+k)+((by+l)*box_per_line_x)]
		            vx_win[box] += vx_tot[(bx+k)+((by+l)*box_per_line_x)]
		            vy_win[box] += vy_tot[(bx+k)+((by+l)*box_per_line_x)]
		            count_box_win[box] += 1
                        else:
                            box                 = bx + (by*box_per_line_x)
                            count_box_win[box] += 1
                            
        #Average win calculus and data print (gnuplot script for velocity)

        module_mean         = 0
        count_busy_box      = 0
        for box in range(box_total):
	    if density_tot[box] > 0 :
	        module_mean       += math.sqrt(vx_tot[box]*vx_tot[box] + vy_tot[box]*vy_tot[box])
	        count_busy_box    += 1
   
        arrow_size = 16*count_busy_box/module_mean
        #create script gnu to plot velocity-density map
        create_gnu_script(arrow_size, box_per_line_x, box_per_column_y, vel_win_file_name, dens_win_file_name, path)
        
        for bx in range(window_size_h+1, box_per_line_x-window_size_h):
            for by in range(window_size_h+1, box_per_column_y-window_size_h):
            	box = bx + (by*box_per_line_x)
                module = math.sqrt((vx_win[box]*vx_win[box]) + (vy_win[box]*vy_win[box]))
	        if density_win[box] > 0.0 and module > 0.0 :
                    normalization = float(image_counter*count_box_win[box])
                    density_win[box]/=normalization
	            dens_win.write("%d %d %f \n" % (bx, by, density_win[box]))
                    vx_win[box]/=normalization
                    vy_win[box]/=normalization
	            vel_win.write("%d %d %f %f %f %f %f \n"% (bx, by, vx_win[box], vy_win[box], module, density_tot[box]/float(count_events), density_win[box]))
	        else :
	            vx_win[box] = 0.0
	            vy_win[box] = 0.0
	            dens_win.write("%d %d %f \n"%(bx, by, 0.0))
                    vel_win.write("%d %d %f %f %f %f %f \n" % (bx, by, 0.0, 0.0, 0.0,0.0, 0.0))
        vel_win.write("e \n")
        vel_win.write("pause -1 \n")
                 
    return vx_win, vy_win, axis_a_win, axis_b_win, density_win




def five_axis(box_total, box_per_line_x, box_per_column_y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type, image_counter):
    caixas_meia_altura    = box_per_column_y/2
    caixas_quarto_altura  = box_per_column_y/4
    caixas_meia_largura   = box_per_line_x/2
    caixas_quarto_largura = box_per_line_x/4
    vx_axis1, vx_axis2, vx_axis3, vx_axis4, vx_axis5, vx_axis6                                                 = [], [], [], [], [], []
    vy_axis1, vy_axis2, vy_axis3, vy_axis4, vy_axis5, vy_axis6                                                 = [], [], [], [], [], []
    axis_a_axis1, axis_a_axis2, axis_a_axis3, axis_a_axis4, axis_a_axis5, axis_a_axis6                         = [], [], [], [], [], []
    axis_b_axis1, axis_b_axis2, axis_b_axis3, axis_b_axis4, axis_b_axis5, axis_b_axis6                         = [], [], [], [], [], []
    ang_elipse_axis1, ang_elipse_axis2, ang_elipse_axis3, ang_elipse_axis4, ang_elipse_axis5, ang_elipse_axis6 = [], [], [], [], [], []

    
    for i in range(box_total):
        if system_type == 'experiment':
            bx = int(i/box_per_column_y)
            by = i%box_per_column_y
        else:
            bx = i%box_per_line_x
            by = int(i/box_per_line_x)
            
        if by == caixas_meia_altura :
            vx_axis1.append(vx_tot[i])
            vy_axis1.append(vy_tot[i])
            axis_a_axis1.append(axis_a_tot[i])
            axis_b_axis1.append(axis_b_tot[i])
            ang_elipse_axis1.append(ang_elipse_tot[i])

        if by == caixas_quarto_altura :
            vx_axis2.append(vx_tot[i]/2.)
            vy_axis2.append(vy_tot[i]/2.)
            axis_a_axis2.append(axis_a_tot[i]/2.)
            axis_b_axis2.append(axis_b_tot[i]/2.)
            ang_elipse_axis2.append(ang_elipse_tot[i]/2.)

        if by == 3*caixas_quarto_altura :
            vx_axis6.append(vx_tot[i]/2.)
            vy_axis6.append(vy_tot[i]/2.)
            axis_a_axis6.append(axis_a_tot[i]/2.)
            axis_b_axis6.append(axis_b_tot[i]/2.)
            ang_elipse_axis6.append(ang_elipse_tot[i]/2.)

        if bx == caixas_meia_largura :
            vx_axis3.append(vx_tot[i]/2.)
            vy_axis3.append(vy_tot[i]/2.)
            axis_a_axis3.append(axis_a_tot[i]/2.)
            axis_b_axis3.append(axis_b_tot[i]/2.)
            ang_elipse_axis3.append(ang_elipse_tot[i]/2.)

        if bx == caixas_meia_largura-caixas_quarto_largura :
            vx_axis4.append(vx_tot[i]/2.)
            vy_axis4.append(vy_tot[i]/2.)
            axis_a_axis4.append(axis_a_tot[i]/2.)
            axis_b_axis4.append(axis_b_tot[i]/2.)
            ang_elipse_axis4.append(ang_elipse_tot[i]/2.)

        if bx == caixas_meia_largura+caixas_quarto_largura :
            vx_axis5.append(vx_tot[i]/2.)
            vy_axis5.append(vy_tot[i]/2.)
            axis_a_axis5.append(axis_a_tot[i]/2.)
            axis_b_axis5.append(axis_b_tot[i]/2.)
            ang_elipse_axis5.append(ang_elipse_tot[i]/2.)


    for i in range(len(vx_axis2)):
            vx_axis2[i] += vx_axis6[i]
            vy_axis2[i] += vy_axis6[i]
            axis_a_axis2[i] += axis_a_axis6[i]
            axis_b_axis2[i] += axis_b_axis6[i]
            ang_elipse_axis2[i] += ang_elipse_axis6[i]

    for i in range(box_per_line_x):
        file_axis1.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis1[i], vy_axis1[i], axis_a_axis1[i], axis_b_axis1[i], ang_elipse_axis1[i]))
        file_axis2.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis2[i], vy_axis2[i], axis_a_axis2[i], axis_b_axis2[i], ang_elipse_axis2[i]))


    for i in range(box_per_column_y):
        file_axis4.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis3[i], vy_axis3[i], axis_a_axis3[i], axis_b_axis3[i], ang_elipse_axis3[i]))
        file_axis3.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis4[i], vy_axis4[i], axis_a_axis4[i], axis_b_axis4[i], ang_elipse_axis4[i]))
        file_axis5.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis5[i], vy_axis5[i], axis_a_axis5[i], axis_b_axis5[i], ang_elipse_axis5[i]))


def imag_count(system_type) :
    counter = 0
    print "Counting images... wait... it may take 5s to count 1000 images\n"
    if system_type == 'superboids' :
        max_number_particles=0
        part_counter=0
        while 1 :
            line          = fn.readline()
            if not line :
                break #EOF
            if line.replace( '\r', '' ) == '\n' : #identifies blank lines
                counter             += 1
                max_number_particles = max(max_number_particles,part_counter)
                part_counter = 0

                while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                    line = fn.readline()
            else:
                part_counter = int(line.split()[0])

    if system_type == 'experiment' :
        while 1 :
            line = file_input_parameter.readline()
            if not line:
                break #EOF
            line_splitted = line.split()
            if line_splitted[0] == 'Time_start:' :
                counter += 1
                
    if system_type == 'szabo-boids' :
        while 1:
            line = fd.readline()
            if not line:
                break
            line_splitted = line.split()
            if line_splitted[0] == 'x' :
                counter += 1

    if system_type == 'vicsek-gregoire' :
        max_number_particles=0
        while 1:
            line = fd.readline()           
            if not line:
                break #EOF
            line_splitted = line.split()
            counter += 1
            n=int(line_splitted[1])
            max_number_particles=max(max_number_particles,n)
            for i in range(n):
                fd.readline()
            
    print "Counted", counter-1, "images.\n"
    print "Type initial and final image number you want to analyse (min=1, max=",counter-1,") - Use spaces to separate the two numbers"
    return max_number_particles

################## Here starts the main program ###############

#Opening input parameter file

file_input_parameter = open("parameter.in")
line_splitted        = file_input_parameter.readline().split()
system_type          = line_splitted[1]
#if system_type == 'superboids' :
#    line_splitted = file_input_parameter.readline().split()
path    = 'output/'+system_type

#Creating the directory structure for output
os.system('mkdir -p %s' % path)
os.system('cp axis12.par %s' % path)
os.system('cp axis345.par %s' % path)
#velocity_density gnuplot file

vid_veloc_dens = open("%s/video_velocity_density.gnu"%path,"w")

#deformation elipse script header

vid_def = open("%s/video_deformation.gnu"%path,"w")
vid_def.write("unset key \n")

#Opening time averages files
dens_win_file_name = "density-win.dat"
dens_win           = open(path+'/'+dens_win_file_name,"w")
vel_win_file_name  = "velocity-win.dat"
vel_win            = open(path+'/'+vel_win_file_name,"w")
def_win            = open("%s/deformation-win.dat"%path,"w")

# Opening five axis analysis files

file_axis1 = open("%s/axis1.dat"%path,"w")
file_axis2 = open("%s/axis2.dat"%path,"w")
file_axis4 = open("%s/axis4.dat"%path,"w")
file_axis3 = open("%s/axis3.dat"%path,"w")
file_axis5 = open("%s/axis5.dat"%path,"w")

if system_type == 'experiment':
    arq_in      = "%s/%s"%(line_splitted[0], line_splitted[1])
    print "You analise an", system_type, "system, reading data from file:\n", arq_in
    window_size = int(line_splitted[2])
    #Opening the data file
    file_input_parameter = open(arq_in)
    imag_count(system_type)
    file_input_parameter.close()
    line_splitted        = sys.stdin.readline().split()
    image_0              = int(line_splitted[0])
    image_f              = int(line_splitted[1])
    file_input_parameter = open(arq_in)
    
    line_splitted        = ['0']
    
    # Reading file head (I have taken some lines of the header, you may want others)
    while(line_splitted[0] != 'X') : #'X' marks the line just before data in experiment data file, that is, the end of the header
        line_splitted = file_input_parameter.readline().split()
        if(line_splitted[0] == 'Box_end:') :
            box_per_column_y, box_per_line_x = int(line_splitted[1]), int(line_splitted[2])
        if(line_splitted[0] == 'Box_size:') :
            box_size = int(line_splitted[1])/4
        if(line_splitted[0] == 'Obstacle_diameter:') :
            R_OBST   = float(line_splitted[1])
            R_OBST   = R_OBST/2.
        if(line_splitted[0] == 'Maximum_observed_velocity:') :
            speed    = float(line_splitted[1])
        if line_splitted[0] == 'Obstacle_position:' :
            X_OBST   = int(line_splitted[1])
            Y_OBST   = int(line_splitted[2])

        x0, xf, y0, yf=0., float(box_per_line_x), 0., float(box_per_column_y)
            
    box_total, ratio, vx_tot, vy_tot, density_tot, axis_a_tot, axis_b_tot, ang_elipse_tot = \
    box_variables_definition_experiment(box_per_column_y, box_per_line_x)


    #Reading x,y,density,vx,vy data on experiment file

    image = 0
    v0    = 1
    while 1 :
        line          = file_input_parameter.readline()
        if not line : break # EOF
        line_splitted = line.split()
        counter       = 0
        image_counter = image_f - image_0
        x, y, axis_a, axis_b, ang_elipse, vx_now, vy_now, density_now=[], [], [], [], [], [], [], []
        while line_splitted[0] != 'Time_start:' :
            if image >= image_0 :
                if line_splitted[7] == 'NaN': 
                    vx_now.append(0.)
                    vy_now.append(0.)
                if line_splitted[7] != 'NaN' : 
                    vx_now.append(float(line_splitted[7]))
                    vy_now.append(float(line_splitted[8]))
                x.append(int(line_splitted[0]) / box_size)
                y.append(int(line_splitted[1]) / box_size)
                axis_a.append(float(line_splitted[4]))
                axis_b.append(float(line_splitted[5]))
                ang_elipse.append(float(line_splitted[6]) / (math.pi) * 180)
                density_now.append(float(line_splitted[9]))
                vx_tot[counter]         += vx_now[counter]
                vy_tot[counter]         += vy_now[counter]
                density_tot[counter]    += density_now[counter]
                axis_a_tot[counter]     += axis_a[counter]
                ang_elipse_tot[counter] += ang_elipse[counter]
                counter                 += 1
            line = file_input_parameter.readline()
            if not line : break # EOF
            line_splitted = line.split()
        image += 1
        line = file_input_parameter.readline()
        line = file_input_parameter.readline()

        if image < image_0 : print "Skipping image ",image
        if image > image_0 and image <= image_f :
            print "Analising image ",image, "..."

            #Function call to write velocity-density gnu script
            velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0, x0, y0, xf, yf)
            #Function call to write deformation elipse gnu script
            deformation_elipsis_script(x, y, axis_b, axis_a, ang_elipse, system_type)
        elif image > image_f:
            break
        

if system_type == "superboids":
    line_splitted = file_input_parameter.readline().split()
    arq_header_in = "%s/%s.dat"%(system_type,line_splitted[1])
    arq_data_in   = "%s/%s_plainprint.dat"%(system_type,line_splitted[1])
    arq_neigh_in  = "%s/%s_neighbors.dat"%(system_type,line_splitted[1])
    print "\nYou analise a", system_type, "system, data is read from files:\n", arq_header_in," (header)\n", arq_data_in," (data)\n", arq_neigh_in," (neighbors)"
    fh            = open(arq_header_in)
    fd            = open(arq_data_in)
    fn            = open(arq_neigh_in)
    line_splitted = file_input_parameter.readline().split()
    window_size   = int(line_splitted[1])
    max_number_particles=imag_count(system_type)
    fn.close()
    fn            = open(arq_neigh_in)
    line_splitted = sys.stdin.readline().split()
    image_0       = int(line_splitted[0])
    image_f       = int(line_splitted[1])
    v0            = 0.007
    part=list(particle(i) for i in range(max_number_particles))

    # Reading superboids parameter file

    while 1 :
        line = fh.readline()
        if not line : break # EOF
        if line.replace( '\r' , '' ) == "\n" :
            continue
        else :
            line_splitted = line.split()
            if line_splitted[0] == '#' and line_splitted[1]=='RECTANGLE:':
                line_splitted = fh.readline().split()
                Lx, Ly        = int(float(line_splitted[1])), int(float(line_splitted[2]))
                Lx            = 240 # corrigindo por enquanto o erro no arquivo de entrada. REVISAR!!!
            if line_splitted[1] == 'Radial' and line_splitted[2]=='R_Eq:':
                line_splitted = fh.readline().split()
                box_size      = int(float(line_splitted[1]))
            if line_splitted[0] == "#radius:":
                R_OBST        = int(line_splitted[1])
                X_OBST        = float(line_splitted[3])
                Y_OBST        = float(line_splitted[4])
                
    box_per_line_x, box_per_column_y = Lx/box_size, Ly/box_size
    x0, y0 = -Lx/2, -Ly/2
    xf, yf = Lx/2, Ly/2
    
    #defining all matrices
    box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, density_tot, vx_win, vy_win, density_win, axis_a_tot, axis_b_tot, \
    ang_elipse_tot, axis_a_win, axis_b_win, ang_elipse_win = box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)

    
    #Reading superboids plainprint data file
    fat_boids_counter = 0
    image             = 0
    count_events      = 0
    points            = []
    index_particle    = []
    
    while 1 :
        line = fd.readline()
        if not line :
            vid_veloc_dens.write("pause -1 \n")
            break #EOF
        if line.replace( '\r', '' ) == '\n' : #identifies blank lines
            if image > image_0 and image <= image_f:
                #Blank lines indicate end of an image so we
                #calculate the average velocity over boxes

                for box in range(box_total):
                    if density_now[box] > 0 :
                        vx_now[box]         = vx_now[box] / density_now[box]
		        vy_now[box]         = vy_now[box] / density_now[box]
                        texture_box[box]    = texture_box[box] / density_now[box]
                        ax_a,ax_b,angle     = axis_angle(texture_box[box])
                        axis_a_win[box]     = ax_a
                        axis_b_win[box]     = ax_b
                        ang_elipse_win[box] = angle
                #Function call to write velocity-density gnu script
                velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0, x0, y0, xf, yf)
                #Function call to write deformation elipsis gnu script
                deformation_elipsis_script_simu(box_per_line_x, box_total, axis_b_win, axis_a_win, ang_elipse_win)
            #Summing each box at different times
            if image > image_0 and image <= image_f:
                for box in range(box_total) :
                    density_tot[box] += density_now[box]
                    vx_tot[box]      += vx_now[box]
                    vy_tot[box]      += vy_now[box]
                    
            image      += 1
            vx_now      = list(0. for i in range(box_total))
            vy_now      = list(0. for i in range(box_total))
            density_now = list(0  for i in range(box_total))
            texture_box = list(np.zeros((2,2)) for i in range(box_total))
            if image <= image_0 :
                print "Skippin image:",image 
            elif image <= image_f :
                print "Image number",image,". Number of super particles =", fat_boids_counter
                count_events     += 1
                # Calculus of textures, B and T
                number_particles = len(points)
                if number_particles > 0:
                    points=np.array(points)
                    list_neighbors=delaunay(points)
                    map_focus_region_to_part(points,list_neighbors,index_particle)
                    map(lambda i:i.texture(), part)
                    for i in index_particle:
                        xx  = int((part[i].r[0]-x0) / box_size)
                        yy  = int((part[i].r[1]-y0) / box_size) * box_per_line_x
                        box = xx+yy
                        texture_box[box]+=part[i].M
                    if  count_events > 1 :
                        
                        map(lambda i:i.UT(), part)
                        
                    map(lambda i:i.copy_to_old(), part)


                points            = []
                index_particle    = []
                

            else:
                break
            fat_boids_counter = 0
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = fd.readline()
        else:
            line_splitted = line.split()
            if image >= image_0 and image <= image_f:
                if float(line_splitted[6]) > 0.5 :
                    fat_boids_counter += 1
                    x, y = float(line_splitted[0]), float(line_splitted[1])
                    if x > x0 and x < xf and y > y0 and y < yf:
                        xx  = int((x-x0) / box_size)
                        yy  = int((y-y0) / box_size) * box_per_line_x
                        box = xx + yy
                        vx_now[box]      += float(line_splitted[2])
                        vy_now[box]      += float(line_splitted[3])
                        density_now[box] += 1.0
                        points.append([x,y])
                        index_particle.append(fat_boids_counter-1)
    image_counter = image_f - image_0

if system_type == "szabo-boids":
    arq_data_in = "%s/%s.dat"% (line_splitted[0], line_splitted[1])
    print "\nYou analise a", line_splitted[0], "system, data is read from files:\n", arq_data_in
    fd            = open(arq_data_in)
#    fn=open(arq_neigh_in)
    window_size   = int(line_splitted[2])
    imag_count(system_type)
    fd.close()
    fd            = open(arq_data_in)
    line_splitted = sys.stdin.readline().split()
    image_0       = int(line_splitted[0])
    image_f       = int(line_splitted[1])
    image_counter = image_f-image_0
    image         = 0
    line_counter  = 0
    count_events  = 0
    v0            = 0.1
    #Reading szabo-boids  data file
    while 1 :
        line = fd.readline()
        if not line : break
        line_counter += 1
        if line.replace( '\r' , '' ) == "\n" :
            continue
        line_splitted = line.split()
        if line_counter < 7 :
            if line_splitted[0] == "Number_of_particles:" :             N = int(line_splitted[1])
            if line_splitted[0] == "Box_size:" :                 box_size = int(line_splitted[1])
            if line_splitted[0] == "Steps_between_images:" : delta_images = int(line_splitted[1])
            if line_splitted[0] == "Dimensions:" :
                Lx = int(line_splitted[1])
                Ly = int(line_splitted[2])
                box_per_line_x, box_per_column_y = Lx / box_size, Ly / box_size
                x0 = -Lx/2
                y0 = -Ly/2
                xf = Lx/2
                yf = Ly/2
                box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, density_tot, vx_win, vy_win, density_win, axis_a_tot, axis_b_tot, \
                ang_elipse_tot, axis_a_win, axis_b_win, ang_elipse_win = box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)


            if line_splitted[0] == "Radius:" : R_OBST = int(line_splitted[1])
            if line_splitted[0] == "Obst_position:" :
                X_OBST = int(line_splitted[1])
                Y_OBST = int(line_splitted[2])
        if line_counter > 6 :
            if image <= image_0 :
                print "Skipping image:",image 
                while line_splitted[0] != 'x' :
                    line = fd.readline()
                    if not line : break
                    line_splitted = line.split()
                image += 1
            elif image <= image_f :
                boids_counter = 0
                vx_now        = list(0. for i in range(box_total))
                vy_now        = list(0. for i in range(box_total))
                density_now = list(0 for i in range(box_total))
                while line_splitted[0] != 'x' :
                    boids_counter += 1
                    x, y = float(line_splitted[0]), float(line_splitted[1])
                    if x > x0 and x < xf and y > y0 and y < yf :
                        xx   = int((x-x0) / box_size)
                        yy   = int((y-y0) / box_size) * box_per_line_x
                        box  = xx + yy
                        vx_now[box]      += float(line_splitted[2])
                        vy_now[box]      += float(line_splitted[3])
                        density_now[box] += 1.0
                        line = fd.readline()
                        if not line : break
                        line_splitted = line.split()
                image        += 1
                count_events += 1
                #Calculate the average velocity over boxes
                for box in range(box_total):
                    if density_now[box] > 0 :
                        vx_now[box] = vx_now[box] / density_now[box]
		        vy_now[box] = vy_now[box] / density_now[box]
                        #Function call to write velocity-density gnu script
                velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0, x0, y0, xf, yf)
                #Summing each box at different times
                for box in range(box_total) :
                    density_tot[box] += density_now[box]
                    vx_tot[box] += vx_now[box]
                    vy_tot[box] += vy_now[box]

            else:
                break

if system_type == "vicsek-gregoire":

    window_size, time_0, time_f, obstacle, voronoi, x0, xf, y0, yf, filename = read_param(file_input_parameter)
    file_input_parameter.close()
    aux           = "%s/include/%s"% (system_type, filename)
    file_par_simu = open(aux)
    Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size = read_param_vic_greg(file_par_simu)
    file_par_simu.close()

    #definindo as caixas e as matrizes
    box_per_line_x, box_per_column_y = int((xf-x0) / box_size), int((yf-y0) / box_size)
    box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, density_tot, vx_win, vy_win, density_win, axis_a_tot, axis_b_tot, ang_elipse_tot, \
    axis_a_win, axis_b_win, ang_elipse_win = box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)
    #arquivo de posicoes e velocidades
    arq_data_in = "%s/data/posicoes.dat"% (system_type)
    print "\nYou analise a", system_type, "system, data is read from files:\n", arq_data_in
    fd            = open(arq_data_in)
    max_number_particles=imag_count(system_type) #conta o numero de imagens
    fd.close()
    fd            = open(arq_data_in)  #reabre o arquivo para leituras das posicoes e vel.
    line_splitted = sys.stdin.readline().split() #le da linha de comando o intervalo de imagens desejado
    image_0       = int(line_splitted[0])
    image_f       = int(line_splitted[1])
    image_counter = image_f - image_0
    image         = 1
    line_counter  = 0
    count_events  = 0
    v0            = 0.05
    part=list(particle(i) for i in range(max_number_particles))
    while 1 :
        line   = fd.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        nlines        = int(line_splitted[1])
        if image <= image_0 :
            print "Skipping image:",image
            for i in range(nlines) :
                fd.readline()
            image += 1
        elif image <= image_f :
            print "Reading image:",image
            vx_now = list(0. for i in range(box_total))
            vy_now = list(0. for i in range(box_total))
            density_now = list(0. for i in range(box_total))
            texture_box = list(np.zeros((2,2)) for i in range(box_total))
            points = []
            index_particle = []
            for i in range(nlines) :
                line_splitted    = fd.readline().split()
                x, y = float(line_splitted[0]), float(line_splitted[1])
                if x > x0 and x < xf and y > y0 and y < yf : 
                    xx  = int((x-x0) / box_size)
                    yy  = int((y-y0) / box_size) * box_per_line_x
                    box = xx+yy
                    vx_now[box]      += float(line_splitted[2])
                    vy_now[box]      += float(line_splitted[3])
                    density_now[box] += 1.0
                    points.append([x,y])
                    index_particle.append(i)
            image        += 1
            count_events += 1
            
            # Calculus of textures, B and T##################
            number_particles = len(points)
            if number_particles > 0:
                points=np.array(points)
                list_neighbors=delaunay(points)
                map_focus_region_to_part(points,list_neighbors,index_particle)
                map(lambda i:i.texture(), part)
                for i in index_particle:
                    xx  = int((part[i].r[0]-x0) / box_size)
                    yy  = int((part[i].r[1]-y0) / box_size) * box_per_line_x
                    box = xx+yy
                    texture_box[box]+=part[i].M
                if  count_events > 1 :
                        
                    map(lambda i:i.UT(), part)
                        
                map(lambda i:i.copy_to_old(), part)

                #Calculate the average velocity over boxes
            for box in range(box_total):
                if density_now[box] > 0 :
                    vx_now[box] = vx_now[box] / density_now[box]
	            vy_now[box] = vy_now[box] / density_now[box]
                    #Function call to write velocity-density gnu script
            velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0, x0, y0, xf, yf)
                #Summing each box at different times
            for box in range(box_total) :
                density_tot[box] += density_now[box]
                vx_tot[box]      += vx_now[box]
                vy_tot[box]      += vy_now[box]

        else:
                break

# Before starting time averages we exclude box at the borders. 
r_obst = R_OBST / box_size
x_obst = (X_OBST-x0) / box_size
y_obst = (Y_OBST-y0) / box_size

zero_borders_and_obstacle(box_per_line_x, box_per_column_y, r_obst, x_obst, y_obst, density_tot, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type)

if system_type == 'experiment':

    # Here we write the time averages of density, velocity and deformation elipse for experiment
    average_density_velocity_deformation_experiment(box_per_line_x, box_per_column_y, x, y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, image_counter)

    
    # Five axis analysis for experiment
    five_axis(box_total, box_per_line_x, box_per_column_y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type, image_counter)
    
else:
    # Here we write the time averages of density, velocity and deformation elipse for simus
    vx_win, vy_win, axis_a_win, axis_b_win, density_win = average_density_velocity_deformation(box_per_line_x, box_per_column_y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, \
                                                                                               density_tot, vx_win, vy_win,axis_a_win, axis_b_win, density_win, count_events, v0, vel_win_file_name, dens_win_file_name, path, image_counter)

    # Five axis analysis for simulations
    five_axis(box_total, box_per_line_x, box_per_column_y, vx_win, vy_win, axis_a_win, axis_b_win, ang_elipse_win, system_type, image_counter)


file_input_parameter.close()
vid_veloc_dens.close()
vid_def.close()
dens_win.close()
vel_win.close()
def_win.close()
file_axis1.close()
file_axis2.close()
file_axis4.close()
file_axis3.close()
file_axis5.close()



        
