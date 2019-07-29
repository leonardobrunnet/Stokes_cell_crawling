import copy
import math as math
import os
import sys
from scipy.spatial import Delaunay
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

#************************************************************
# diagonalization of matrices. Input: Matrix (numpy 2x2)
# output: axis_a, axis_b and angle in degrees
def matrix_three_component(M):
    xxpyy=M.trace()/2.
    xxmyy=(M[0][0]-M[1][1])/2.
    xy = M[0][1]
    return xxpyy,xxmyy,xy

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
                self.M+=self.calc_m(i)
            self.M/=n
            
    def calc_m(self,i):
        l=self.r-part[i].r
        m=np.outer(l,l)
        return m
    
    def copy_to_old(self):
        self.list_neigh_old=copy.deepcopy(self.list_neigh)
        self.r_old=copy.deepcopy(self.r)

    #function averaging l in t and t+dt
    #equations C2 and C6 of graner tools 
    def l_av_dl(self,i):
        l=self.r-part[i].r
        l_old=self.r_old-part[i].r_old
        lav=(l+l_old)/2.
        dl=l-l_old
        return lav,dl
        
    def calc_B_and_T(self):
        list_c=list(set(self.list_neigh).intersection(self.list_neigh_old))
        list_a=list(set(self.list_neigh).difference(self.list_neigh_old))
        list_d=list(set(self.list_neigh_old).difference(self.list_neigh))
        self.zeros()
        Nc=len(list_c) #numero de links conservados
        Na=len(list_a) #numero de links adquiridos
        Nd=len(list_d) #numero de links desaparecidos
        Ntot=Nc+Na+Nd
        if Ntot != 0 :
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
                ma=self.calc_m(i)
                self.T+=ma
            for i in list_d:
                md=self.calc_m(i)
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



def delaunay(points,max_dist):
    tri = Delaunay(points)
    # x,y=[],[]
    # z,zz=[],[]
    # for i,w in enumerate(points):
    #     if i%50 == 0:
    #         x.append(w[0])
    #         y.append(w[1])
    #     else :
    #         z.append(w[0])
    #         zz.append(w[1])
    # fig=plt.scatter(x,y,s=30,c='b')
    # fig=plt.scatter(z,zz,s=30,c='g')
    # plt.show()
    list_neigh = [ [] for i in range(len(points)) ]
    for i in tri.simplices:
        for j in i:
            for l in i:
                if l != j :
                    if np.linalg.norm(points[j]-points[l]) < max_dist : #
                        if l not in list_neigh[j]:
                            list_neigh[j].append(l)
                        if j not in list_neigh[l]:
                            list_neigh[l].append(j)
    # x,y=[],[]
    # for i,w in enumerate(list_neigh) :
    #     if i%50==0 :
    #         for j in w :
    #             x.append(points[j][0])
    #             y.append(points[j][1])
    # fig=plt.scatter(x,y,s=30,c='r')
    # fig=plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
    # plt.show()
    # exit()

    return list_neigh

def create_gnu_script_fluct_vel(arrow_size, box_per_line_x, box_per_column_y, vel_fluct_win_file_name, dens_win_file_name, path):
    proportion_x, proportion_y             = 1.0, 0.7
    grid_x, grid_y, levels                 = 200, 200, 4
    image_resolution_x, image_resolution_y = 1024, 1024
    name_output_map            = "density-velocity-fluct.png"
    file_script_den_vel_fluct  = open(path+"/scriptdenvel_fluct.gnu","w")
   
    file_script_den_vel_fluct.write("set size %1.2f,%1.2f \n"% (proportion_x, proportion_y))
    file_script_den_vel_fluct.write("set palette defined ( 0 '#000000',\\\n")
    file_script_den_vel_fluct.write("                      1 '#0000ff',\\\n")
    file_script_den_vel_fluct.write("                      2 '#00ffff',\\\n")
    file_script_den_vel_fluct.write("                      3 '#00ff00',\\\n")
    file_script_den_vel_fluct.write("                      4 '#ffff00',\\\n")
    file_script_den_vel_fluct.write("                      5 '#ff0000')\n")
    file_script_den_vel_fluct.write("set nokey \n")
    file_script_den_vel_fluct.write("set dgrid3d %d,%d,%d \n"% (grid_x, grid_y, levels))
    file_script_den_vel_fluct.write("set pm3d explicit \n")
    file_script_den_vel_fluct.write("set output \"| head -n -2 > toto.dat\" \n")
    file_script_den_vel_fluct.write("set table \n")
    file_script_den_vel_fluct.write("splot \"%s\" using 1:2:3 \n"% dens_win_file_name)
    file_script_den_vel_fluct.write("unset table \n")
    file_script_den_vel_fluct.write("unset dgrid3d \n")
    file_script_den_vel_fluct.write("mtf = %f \n"% arrow_size)
    file_script_den_vel_fluct.write(" \n")
    file_script_den_vel_fluct.write(" \n")
    file_script_den_vel_fluct.write("set pm3d map \n")
    file_script_den_vel_fluct.write("splot [%d:%d][%d:%d] \"toto.dat\" \n"% (0, box_per_line_x, 0, box_per_column_y))
    file_script_den_vel_fluct.write("replot \"%s\" u($1):($2):(0.0):(mtf*$3):(mtf*$4):(0.0) with vectors head size 1.5,5,60 lt rgb \"black\" \n"% vel_fluct_win_file_name)
#    file_script_den_vel_fluct.write("pause -1 \n")
    file_script_den_vel_fluct.write("set terminal pngcairo  size %d,%d enhanced font 'Verdana, 14' crop\n"% (image_resolution_x, image_resolution_y)) 
    file_script_den_vel_fluct.write("set output \"%s\" \n"% name_output_map)
    file_script_den_vel_fluct.write("replot \n")  
    
def create_gnu_script(arrow_size, box_per_line_x, box_per_column_y, vel_win_file_name, dens_win_file_name, path):
    proportion_x, proportion_y             = 1.0, 0.7
    grid_x, grid_y, levels                 = 200, 200, 4
    image_resolution_x, image_resolution_y = 1024, 1024
    name_output_map            = "density-velocity.png"
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
    file_script_den_vel.write("set output \"| head -n -2 > toto1.dat\" \n")
    file_script_den_vel.write("set table \n")
    file_script_den_vel.write("splot \"%s\" using 1:2:3 \n"% dens_win_file_name)
    file_script_den_vel.write("unset table \n")
    file_script_den_vel.write("unset dgrid3d \n")
    file_script_den_vel.write("mtf = %f \n"% arrow_size)
    file_script_den_vel.write(" \n")
    file_script_den_vel.write(" \n")
    file_script_den_vel.write("set pm3d map \n")
    file_script_den_vel.write("splot [%d:%d][%d:%d] \"toto1.dat\" \n"% (0, box_per_line_x, 0, box_per_column_y))
    file_script_den_vel.write("replot \"%s\" u($1):($2):(0.0):(mtf*$3):(mtf*$4):(0.0) with vectors head size 1.5,5,60 lt rgb \"black\" \n"% vel_win_file_name)
#    file_script_den_vel.write("pause -1 \n")
    file_script_den_vel.write("set terminal pngcairo  size %d,%d enhanced font 'Verdana, 14' crop\n"% (image_resolution_x, image_resolution_y))
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
            max_dist = box_size
        if line_splitted[1] == 'SNAPSHOT' :
            Delta_t = int(line_splitted[2])
        if line_splitted[1] == 'V1' :
            v0 = float(line_splitted[2])
    return Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size, Delta_t, v0, max_dist

def read_param(file_input_parameter) :
    box_mag = 1.0
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
        if line_splitted[0] == 'obstacle' :
            obstacle = line_splitted[1]
        if line_splitted[0] == 'x0' :
            x0 = int(line_splitted[1])
        if line_splitted[0] == 'xf' :
            xf = int(line_splitted[1])
        if line_splitted[0] == 'file' :
            filename = line_splitted[1]
        if line_splitted[0] == 'box_magnification' :
            box_mag = float(line_splitted[1])
    return window_size, time_0, time_f, obstacle, x0, xf, filename, box_mag


def read_param_potts(file_par_simu) :
    while 1 :
        line = file_par_simu.readline()
        if not line:
            break #EOF
        if line == "":
            line = file_par_simu.readline()
            if not line:
                break #EOF
        if line.replace( '\r', '' ) == '\n' : #identifies blank lines
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = file_par_simu.readline()
        if line.replace( '\t', '' ) == '\n' : #identifies blank lines
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = file_par_simu.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        if line_splitted[0] == 'multiplicador' :
            multiplier = float(line_splitted[2])
        if line_splitted[0] == 'size_x' :
            Lx = int(line_splitted[4]) * multiplier
        if line_splitted[0] == 'size_y' :
            Ly = int(line_splitted[4]) * multiplier
        if line_splitted[0] == 'size_obst' :
            R_OBST = float(line_splitted[4]) * multiplier
        if line_splitted[0] == 'center_x' :
            X_OBST = float(line_splitted[4]) * multiplier
        if line_splitted[0] == 'center_y' :
            Y_OBST = float(line_splitted[4]) * multiplier
        if line_splitted[0] == 'target_v' :
            box_size = math.sqrt(float(line_splitted[6]))*multiplier
            box_size = 2*box_size #This seems to be the best box_size for potts
            max_dist = box_size
        if line_splitted[0] == 'CompuCell3DElmnt' :
            break
    Delta_t = 100 # this information is not in the simulation files for compucell :/
    v0 = 1.0 # There is no v0 in Potts model
    return Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size, max_dist, Delta_t, v0

def read_param_voronoi(filename):
    X_OBST = 0.0 
    Y_OBST = 0.0
    v0 = 1.0
    # first file to be read make_initial.py
    aux           = system_type+"/"+filename
    file_par_simu = open(aux)
    while 1 :
        line = file_par_simu.readline()
        # the voronoi has the config file 
        if not line:
            break #EOF
        if line.replace( '\r', '' ) == '\n' : #identifies blank lines
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = file_par_simu.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        if line_splitted[0] == 'size_mul':
            size_mul = float(line_splitted[2])
        if line_splitted[0] == 'wall':
            # make sure there are spaces between names and information
            line = line.replace( '*' , ' * ')
            line = line.replace( ',' , ' , ')
            line_splitted = line.split()
            Lx = float(line_splitted[4]) * size_mul
            Ly = float(line_splitted[8]) * size_mul
        if line_splitted[0] == 'circle' :
            line = line.replace( '*' , ' * ')
            line = line.replace( ',' , ' , ')
            line = line.replace( ')' , ' ) ')
            line_splitted = line.split()
            R_OBST = float(line_splitted[17]) * size_mul
    # the rest of the information comes from another file
    file_par_simu.close()
    filename = "config.conf"
    aux           = system_type+"/"+filename
    file_par_simu = open(aux)
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
        if len(line_splitted) > 1:
            if line_splitted[0] == 'dump' and line_splitted[1] == 'cell' :
                line = line.replace( ';' , ' ; ')
                line = line.replace( '=' , ' = ')
                line_splitted = line.split()
                Delta_t = int(line_splitted[13])

            if line_splitted[0] == 'pair_potential' and line_splitted[1] == 'soft':
                line = line.replace( ';' , ' ; ')
                line = line.replace( '}' , ' } ')
                line_splitted = line.split()
                box_size = float(line_splitted[8])
                max_dist = box_size
    return Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size, Delta_t, v0, max_dist


def box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf) :
    box_total              = box_per_column_y * box_per_line_x
    vx_now                 = list(0. for i in range(box_total))
    vy_now                 = list(0. for i in range(box_total))
    density_now            = list(0  for i in range(box_total))
    vx_tot                 = list(0. for i in range(box_total))
    vy_tot                 = list(0. for i in range(box_total))
    vx2_tot                = list(0. for i in range(box_total))
    vy2_tot                = list(0. for i in range(box_total))
    density_tot            = list(0  for i in range(box_total))
    vx_win                 = list(0. for i in range(box_total))
    vy_win                 = list(0. for i in range(box_total))
    vx2_win                = list(0. for i in range(box_total))
    vy2_win                = list(0. for i in range(box_total))
    density_win            = list(0  for i in range(box_total))
    texture_box = list(np.zeros((2,2)) for i in range(box_total))
    B_box = list(np.zeros((2,2)) for i in range(box_total))
    T_box = list(np.zeros((2,2)) for i in range(box_total))
    texture_tot = list(np.zeros((2,2)) for i in range(box_total))
    B_tot = list(np.zeros((2,2)) for i in range(box_total))
    T_tot = list(np.zeros((2,2)) for i in range(box_total))
    texture_win = list(np.zeros((2,2)) for i in range(box_total))
    B_win = list(np.zeros((2,2)) for i in range(box_total))
    T_win = list(np.zeros((2,2)) for i in range(box_total))
    # ang_elipse_tot         = list(0. for i in range(box_total))
    # axis_a_win_texture     = list(0. for i in range(box_total))
    # axis_b_win_texture     = list(0. for i in range(box_total))
    # ang_elipse_win_texture = list(0. for i in range(box_total))
    # axis_a_win_B           = list(0. for i in range(box_total))
    # axis_b_win_B           = list(0. for i in range(box_total))
    # ang_elipse_win_B       = list(0. for i in range(box_total))
    # axis_a_win_T           = list(0. for i in range(box_total))
    # axis_b_win_T           = list(0. for i in range(box_total))
    # ang_elipse_win_T       = list(0. for i in range(box_total))
    #    ratio=float(box_per_column_y)/box_per_line_x
    ratio = float((yf-y0)) / (xf-x0)
#    vid_def.write("set size 1,%f\n"% (ratio))
    image_resolution_x, image_resolution_y=1300,1300*ratio
    vid_def.write("set terminal png large size %d,%d \n"% (image_resolution_x, image_resolution_y)) 
    vid_def.write("set xrange [0:%f]  \n" % box_per_line_x)
    vid_def.write("set yrange [0:%f]  \n" % box_per_column_y)    
    vid_B.write("set terminal png large size %d,%d \n"% (image_resolution_x, image_resolution_y)) 
    vid_B.write("set xrange [0:%f]  \n" % box_per_line_x)
    vid_B.write("set yrange [0:%f]  \n" % box_per_column_y)    
    vid_T.write("set terminal png large size %d,%d \n"% (image_resolution_x, image_resolution_y)) 
    vid_T.write("set xrange [0:%f]  \n" % box_per_line_x)
    vid_T.write("set yrange [0:%f]  \n" % box_per_column_y)    
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

    return box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, vx2_tot, vy2_tot,\
        density_tot, vx_win, vy_win,vx2_win, vy2_win, density_win, texture_box, B_box, \
T_box, texture_tot, B_tot, T_tot, texture_win, B_win, T_win

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
    vid_def.write("set term png  \n")
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

    return box_total, ratio, vx_tot, vy_tot, density_tot, axis_a_tot,\
        axis_b_tot, ang_elipse_tot


def velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0):
    #Here we write each image to the gnuplot velocity-density movie script
    vid_veloc_dens.write("plot [%f:%f] [%f:%f] \'-\' u ($1):($2):(arrow*$3):(arrow*$4):($5) with vectors head size  0.6,20,60  filled palette title \"%d\"\n" % \
    (0, box_per_line_x, 0, box_per_column_y, image))
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

def deformation_elipsis_script(x, y, axis_a, axis_b, ang_elipse, system_type) :
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

def texture_elipsis_script_simu(box_per_line_x, box_total, axis_a, axis_b, ang_elipse, image, points, x0, y0, box_size) :
    #Texture elipsis gnuplot script for simus
    vid_def.write("set output \"text-%d.png\"\n"%image)
    vid_def.write("set multiplot\n")
    vid_def.write("plot \'-\' using 1:2:3:4:5 with ellipses\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            x = i % box_per_line_x + 0.5
            y = i / box_per_line_x + 0.5
            vid_def.write("%f %f 1.0 %f %f \n" % (x,y,axis_b[i]/(axis_a[i]),ang_elipse[i]))
    vid_def.write("e \n")
    
    vid_def.write("set style line 1 lc rgb 'blue' pt 7\n")
    vid_def.write("plot \'-\' w points ls 1 notitle\n")
    points=(points-np.array([x0,y0]))/box_size
    for i in range(len(points)) :
        vid_def.write("%f %f\n"%(points[i][0],points[i][1]))
    vid_def.write("pause .1 \n")
    vid_def.write("e \n")
    #    vid_def.write("unset for [i=1:%i] object i \n" % (box_total+1))
    vid_def.write("unset multiplot \n")    


        
def B_elipsis_script_simu(box_per_line_x, box_total, axis_a, axis_b, ang_elipse, image, points, x0, y0, box_size) :
    #B elipsis gnuplot script for simus
    vid_B.write("set output \"B-%d.png\"\n"%image)
    vid_B.write("set multiplot\n")
    vid_B.write("plot \'-\' using 1:2:3:4:5 with ellipses lt rgb \"red\"\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            if axis_a[i] > 0 and axis_b[i]  >0:
                x = i % box_per_line_x + 0.5
                y = i / box_per_line_x + 0.5
                vid_B.write("%f %f 1.0 %f %f \n" % (x,y,axis_b[i]/(axis_a[i]),ang_elipse[i]))
    vid_B.write("e \n")
    vid_B.write("plot \'-\' using 1:2:3:4:5 with ellipses lt rgb \"blue\"\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            if axis_a[i] < 0 and axis_b[i]  < 0:
                x = i % box_per_line_x + 0.5
                y = i / box_per_line_x + 0.5
                vid_B.write("%f %f 1.0 %f %f \n" % (x,y,axis_a[i]/(axis_b[i]),ang_elipse[i]))
    vid_B.write("e \n")
    vid_B.write("plot \'-\' using 1:2:3:4:5 with ellipses lt rgb \"black\"\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            if axis_a[i] > 0 and axis_b[i]  < 0:
                x = i % box_per_line_x + 0.5
                y = i / box_per_line_x + 0.5
                if abs(axis_a[i])>abs(axis_b[i]):
                    vid_B.write("%f %f 1.0 %f %f \n" % (x,y,-axis_b[i]/(axis_a[i]),ang_elipse[i]))
                else :
                    vid_B.write("%f %f 1.0 %f %f \n" % (x,y,-axis_a[i]/(axis_b[i]),ang_elipse[i]))
    vid_B.write("e \n")
    #Particles position at the present image
    vid_B.write("set style line 1 lc rgb 'blue' pt 7 ps 1\n")
    vid_B.write("plot \'-\' w points ls 1 notitle\n")
    points=(points-np.array([x0,y0]))/box_size
    for i in range(len(points)) :
        if points[i][0] < box_per_line_x:
            vid_B.write("%f %f\n"%(points[i][0],points[i][1]))
    vid_B.write("e \n")

    #Particles position of the previous image

    # vid_B.write("set style line 2 lc rgb 'green' pt 7 ps 1\n")
    # vid_B.write("plot \'-\' w points ls 2 notitle\n")
    # points_old=[]
    # points_old=map(lambda i:(i.r_old-np.array([x0,y0]))/box_size, part)
    # for i in range(len(points_old)) :
    #     if points_old[i][0] < box_per_line_x:
    #         vid_B.write("%f %f\n"%(points_old[i][0],points_old[i][1]))
    # #     vid_B.write("pause .1 \n")
    # vid_B.write("e \n")
    vid_B.write("unset multiplot \n")    
        
def T_elipsis_script_simu(box_per_line_x, box_total, axis_a, axis_b, ang_elipse, image, points, x0, y0, box_size) :
    #T elipsis gnuplot script for simus, with particle center positions in two consecutive images
    vid_T.write("set output \"T-%d.png\"\n"%image)
    vid_T.write("set multiplot\n")

    vid_T.write("plot \'-\' using 1:2:3:4:5 with ellipses lt rgb \"red\"\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            if axis_a[i] > 0 and axis_b[i]  >0:
                x = i % box_per_line_x + 0.5
                y = i / box_per_line_x + 0.5
                vid_T.write("%f %f 1.0 %f %f \n" % (x,y,axis_b[i]/(axis_a[i]),ang_elipse[i]))
    vid_T.write("e \n")
    vid_T.write("plot \'-\' using 1:2:3:4:5 with ellipses lt rgb \"blue\"\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            if axis_a[i] < 0 and axis_b[i]  < 0:
                x = i % box_per_line_x + 0.5
                y = i / box_per_line_x + 0.5
                vid_T.write("%f %f 1.0 %f %f \n" % (x,y,axis_a[i]/(axis_b[i]),ang_elipse[i]))
    vid_T.write("e \n")
    vid_T.write("plot \'-\' using 1:2:3:4:5 with ellipses lt rgb \"black\"\n")
    for i in range(box_total) :
        if axis_a[i] != 0 :
            if axis_a[i] > 0 and axis_b[i]  < 0:
                x = i % box_per_line_x + 0.5
                y = i / box_per_line_x + 0.5
                if abs(axis_a[i])>abs(axis_b[i]):
                    vid_T.write("%f %f 1.0 %f %f \n" % (x,y,-axis_b[i]/(axis_a[i]),ang_elipse[i]))
                else :
                    vid_T.write("%f %f 1.0 %f %f \n" % (x,y,-axis_a[i]/(axis_b[i]),ang_elipse[i]))
    vid_T.write("e \n")

    #Particles position at the present image
    vid_T.write("set style line 1 lc rgb 'blue' pt 7 ps 1\n")
    vid_T.write("plot \'-\' w points ls 1 notitle\n")
    points=(points-np.array([x0,y0]))/box_size
    for i in range(len(points)) :
        vid_T.write("%f %f\n"%(points[i][0],points[i][1]))
    vid_T.write("pause .1 \n")
    vid_T.write("e \n")

    # #Particles position of the previous image

    # vid_T.write("set style line 2 lc rgb 'green' pt 7 ps 1\n")
    # vid_T.write("plot \'-\' w points ls 2 notitle\n")
    # points_old=map(lambda i:(i.r_old-np.array([x0,y0]))/box_size, part)
    # for i in range(len(points_old)) :
    #     vid_T.write("%f %f\n"%(points_old[i][0],points_old[i][1]))
    # vid_T.write("pause .1 \n")
    # vid_T.write("e \n")
    # vid_T.write("unset multiplot \n")    
    vid_T.write("unset multiplot \n")    

def  zero_borders_and_obstacle_experiment(box_per_line_x, box_per_column_y, r_obst, x_obst, y_obst, density_tot, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type) :
    center_x = box_per_line_x/2
    center_y = box_per_column_y/2
    for i in range(box_total):
        bx = int(i/box_per_column_y)
        by = i%box_per_column_y
        if bx == 0 or bx == box_per_line_x-1 or by == 0 or by == box_per_column_y-1 :
            density_tot[i]    = -10
            vx_tot[i]         = 0.
            vy_tot[i]         = 0.
            axis_a_tot[i]     = 0.
            axis_b_tot[i]     = 0.
            ang_elipse_tot[i] = 0.
        if math.sqrt((bx-x_obst)**2 + (by-y_obst)**2) <= r_obst :
            density_tot[i]    = -10
            vx_tot[i]         = 0.
            vy_tot[i]         = 0.
            axis_a_tot[i]     = 0.
            axis_b_tot[i]     = 0.
            ang_elipse_tot[i] = 0.
    return box_per_line_x, box_per_column_y, density_tot, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot
                              
def  zero_borders_and_obstacle_simu(box_per_line_x, box_per_column_y, r_obst, x_obst, y_obst, density_tot, vx_tot, vy_tot, texture_tot, B_tot, T_tot, system_type) :
    center_x = box_per_line_x/2
    center_y = box_per_column_y/2
    for i in range(box_total):
        bx = i%box_per_line_x
        by = int(i/box_per_line_x)
        caixas_quarto_altura = box_per_column_y/4
        if bx == 0 or bx == box_per_line_x-1 or by == 0 or by == box_per_column_y-1 :
            density_tot[i]    = -10
            vx_tot[i]         = 0.
            vy_tot[i]         = 0.
            texture_tot[i]    = np.zeros((2,2))
            B_tot[i]          = np.zeros((2,2))
            T_tot[i]          = np.zeros((2,2))
        if system_type == 'superboids' :
            if math.sqrt((bx-center_x)**2 + (by-center_y)**2) < int(r_obst) :
                density_tot[i]    = -10
                vx_tot[i]         = 0.
                vy_tot[i]         = 0.
                texture_tot[i]    = np.zeros((2,2))
                B_tot[i]          = np.zeros((2,2))
                T_tot[i]          = np.zeros((2,2))
        if system_type == 'vicsek-gregoire' :
            if math.sqrt((bx-x_obst)**2 + (by-y_obst)**2) <= int(r_obst):
                density_tot[i]    = -10
                vx_tot[i]         = 0.
                vy_tot[i]         = 0.
                texture_tot[i]    = np.zeros((2,2))
                B_tot[i]          = np.zeros((2,2))
                T_tot[i]          = np.zeros((2,2))

    return box_per_line_x, box_per_column_y, density_tot, vx_tot, vy_tot, texture_tot, B_tot, T_tot

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

def average_density_velocity_deformation(box_per_line_x, box_per_column_y, vx_tot, vy_tot,  \
    density_tot, texture_tot, B_tot, T_tot, vx_win, vy_win, vx2_win, vy2_win,  density_win, texture_win, B_win, \
    T_win, count_events, v0, vel_win_file_name, vel_fluct_win_file_name, dens_win_file_name, path, image_counter, window_size) :


    box_total         = box_per_column_y*box_per_line_x
    window_size_h     = window_size/2

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
		        vx2_win[box] += vx_tot[(bx+k)+((by+l)*box_per_line_x)]**2
		        vy2_win[box] += vy_tot[(bx+k)+((by+l)*box_per_line_x)]**2
                        texture_win[box] += texture_tot[(bx+k)+((by+l)*box_per_line_x)]
                        B_win[box] += B_tot[(bx+k)+((by+l)*box_per_line_x)]
                        T_win[box] += T_tot[(bx+k)+((by+l)*box_per_line_x)]
		        count_box_win[box] += 1

    #Average win calculus and data print (gnuplot script for velocity)
    module_mean         = 0
    count_busy_box      = 0
    for box in range(box_total):
	if density_win[box] > 0 :
	    module_mean       += math.sqrt(vx_win[box]*vx_win[box] + vy_win[box]*vy_win[box])
	    count_busy_box    += density_win[box]
   
    arrow_size = count_busy_box/module_mean
    #create script gnu to plot velocity-density
    create_gnu_script(arrow_size, box_per_line_x, box_per_column_y, vel_win_file_name, dens_win_file_name, path)
    create_gnu_script_fluct_vel(arrow_size, box_per_line_x, box_per_column_y, vel_fluct_win_file_name, dens_win_file_name, path)
    ells=[]
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
                vx2_win[box]/=normalization
                vy2_win[box]/=normalization
                vx2_win[box]=vx2_win[box]-vx_win[box]**2
                vy2_win[box]=vy2_win[box]-vy_win[box]**2
                texture_win[box]/=normalization
                B_win[box]/=normalization
                T_win[box]/=normalization
	        vel_win.write("%d %d %f %f %f %f %f \n"% (bx, by, vx_win[box], vy_win[box], module, density_tot[box]/float(count_events), density_win[box]))
	        vel_fluct_win.write("%d %d %f %f %f %f %f \n"% (bx, by, vx2_win[box], vy2_win[box], module, density_tot[box]/float(count_events), density_win[box]))
#                print box, texture_win[box]
                axis_a,axis_b, ang_elipse = axis_angle(texture_win[box])
                ells.append(Ellipse(np.array([bx,by]),1.,axis_b/axis_a,ang_elipse))
	        texture_win_file.write("%d %d %f %f %f \n"% (bx, by, axis_a, axis_b, ang_elipse))
                axis_a,axis_b, ang_elipse = axis_angle(B_win[box])
	        B_win_file.write("%d %d %f %f %f \n"% (bx, by, axis_a, axis_b, ang_elipse))
                axis_a,axis_b, ang_elipse = axis_angle(T_win[box])
	        T_win_file.write("%d %d %f %f %f \n"% (bx, by, axis_a, axis_b, ang_elipse))

	    else :
	        vx_win[box] = 0.0
	        vy_win[box] = 0.0
	        dens_win.write("%d %d %f \n"%(bx, by, 0.0))
                vel_win.write("%d %d %f %f %f %f %f \n" % (bx, by, 0.0, 0.0, 0.0,0.0, 0.0))
                vel_fluct_win.write("%d %d %f %f %f %f %f \n" % (bx, by, 0.0, 0.0, 0.0,0.0, 0.0))
    vel_win.write("e \n")
    vel_win.write("pause -1 \n")

    #Windowed elipsis for texture
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
#        e.set_alpha(np.random.rand())
#        e.set_facecolor(np.random.rand(3))
    ax.set_xlim(0,box_per_line_x)
    ax.set_ylim(0,box_per_column_y)
    plt.savefig(path+"/texture_win.png",bbox_inches="tight")
    return vx_win, vy_win, vx2_win, vy2_win,  density_win, texture_win, B_win, T_win

def five_axis_experiment(box_total, box_per_line_x, box_per_column_y, vx_win, vy_win, axis_a_win, axis_b_win, ang_elipse_win, system_type, image_counter):

    #Attention in experiment no windowing average is done! That's why this function is called using _tot but _win is used.
    
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
        bx = int(i/box_per_column_y)
        by = i%box_per_column_y

        if by == caixas_meia_altura :
            vx_axis1.append(vx_win[i])
            vy_axis1.append(vy_win[i])
            axis_a_axis1.append(axis_a_win[i])
            axis_b_axis1.append(axis_b_win[i])
            ang_elipse_axis1.append(ang_elipse_win[i])

        if by == caixas_quarto_altura :
            vx_axis2.append(vx_win[i])
            vy_axis2.append(vy_win[i])
            axis_a_axis2.append(axis_a_win[i])
            axis_b_axis2.append(axis_b_win[i])
            ang_elipse_axis2.append(ang_elipse_win[i])

        if by == 3*caixas_quarto_altura :
            vx_axis6.append(vx_win[i])
            vy_axis6.append(vy_win[i])
            axis_a_axis6.append(axis_a_win[i])
            axis_b_axis6.append(axis_b_win[i])
            ang_elipse_axis6.append(ang_elipse_win[i])

        if bx == caixas_meia_largura :
            vx_axis3.append(vx_win[i])
            vy_axis3.append(vy_win[i])
            axis_a_axis3.append(axis_a_win[i])
            axis_b_axis3.append(axis_b_win[i])
            ang_elipse_axis3.append(ang_elipse_win[i])

        if bx == caixas_meia_largura-caixas_quarto_largura :
            vx_axis4.append(vx_win[i])
            vy_axis4.append(vy_win[i])
            axis_a_axis4.append(axis_a_win[i])
            axis_b_axis4.append(axis_b_win[i])
            ang_elipse_axis4.append(ang_elipse_win[i])

        if bx == caixas_meia_largura+caixas_quarto_largura :
            vx_axis5.append(vx_win[i])
            vy_axis5.append(vy_win[i])
            axis_a_axis5.append(axis_a_win[i])
            axis_b_axis5.append(axis_b_win[i])
            ang_elipse_axis5.append(ang_elipse_win[i])


#    for i in range(len(vx_axis2)):
#            vx_axis2[i] += vx_axis6[i]
#            vy_axis2[i] += vy_axis6[i]
#            axis_a_axis2[i] += axis_a_axis6[i]
#            axis_b_axis2[i] += axis_b_axis6[i]
#            ang_elipse_axis2[i] += ang_elipse_axis6[i]

    for i in range(box_per_line_x):
        file_axis1.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis1[i], vy_axis1[i], axis_a_axis1[i], axis_b_axis1[i], ang_elipse_axis1[i]))
        file_axis2.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis2[i], vy_axis2[i], axis_a_axis2[i], axis_b_axis2[i], ang_elipse_axis2[i]))
        file_axis6.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis6[i], vy_axis6[i], axis_a_axis6[i], axis_b_axis6[i], ang_elipse_axis6[i]))


    for i in range(box_per_column_y):
        file_axis3.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis3[i], vy_axis3[i], axis_a_axis3[i], axis_b_axis3[i], ang_elipse_axis3[i]))
        file_axis4.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis4[i], vy_axis4[i], axis_a_axis4[i], axis_b_axis4[i], ang_elipse_axis4[i]))
        file_axis5.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis5[i], vy_axis5[i], axis_a_axis5[i], axis_b_axis5[i], ang_elipse_axis5[i]))

        


def five_axis_simu(box_total, box_per_line_x, box_per_column_y, vx_win, vy_win, texture_win, B_win, T_win, system_type, image_counter,path) :

    caixas_meia_altura    = box_per_column_y/2
    caixas_quarto_altura  = box_per_column_y/4
    caixas_meia_largura   = box_per_line_x/2
    caixas_quarto_largura = box_per_line_x/4
    vx_axis1, vx_axis2, vx_axis3, vx_axis4, vx_axis5, vx_axis6                                                 = [], [], [], [], [], []
    vy_axis1, vy_axis2, vy_axis3, vy_axis4, vy_axis5, vy_axis6                                                 = [], [], [], [], [], []
    texture_axis_a_axis1, texture_axis_a_axis2, texture_axis_a_axis3, texture_axis_a_axis4, texture_axis_a_axis5, texture_axis_a_axis6                         = [], [], [], [], [], []
    texture_axis_b_axis1, texture_axis_b_axis2, texture_axis_b_axis3, texture_axis_b_axis4, texture_axis_b_axis5, texture_axis_b_axis6                         = [], [], [], [], [], []
    texture_ang_elipse_axis1, texture_ang_elipse_axis2, texture_ang_elipse_axis3, texture_ang_elipse_axis4, texture_ang_elipse_axis5, texture_ang_elipse_axis6 = [], [], [], [], [], []
    B_axis_a_axis1, B_axis_a_axis2, B_axis_a_axis3, B_axis_a_axis4, B_axis_a_axis5, B_axis_a_axis6                         = [], [], [], [], [], []
    B_axis_b_axis1, B_axis_b_axis2, B_axis_b_axis3, B_axis_b_axis4, B_axis_b_axis5, B_axis_b_axis6                         = [], [], [], [], [], []
    B_ang_elipse_axis1, B_ang_elipse_axis2, B_ang_elipse_axis3, B_ang_elipse_axis4, B_ang_elipse_axis5, B_ang_elipse_axis6 = [], [], [], [], [], []
    T_axis_a_axis1, T_axis_a_axis2, T_axis_a_axis3, T_axis_a_axis4, T_axis_a_axis5, T_axis_a_axis6                         = [], [], [], [], [], []
    T_axis_b_axis1, T_axis_b_axis2, T_axis_b_axis3, T_axis_b_axis4, T_axis_b_axis5, T_axis_b_axis6                         = [], [], [], [], [], []
    T_ang_elipse_axis1, T_ang_elipse_axis2, T_ang_elipse_axis3, T_ang_elipse_axis4, T_ang_elipse_axis5, T_ang_elipse_axis6 = [], [], [], [], [], []

#To avoid rewriting the whole function: axis_a is xxpyy=(xx+yy)/2,
#axis_b is xxmyy=(xx-yy)/2 and ang = xy. These are components of a
# symmetric matrix        xx xy 
#                    M = [     ]
#                         xy yy  

    for i in range(box_total):
        bx = i%box_per_line_x
        by = int(i/box_per_line_x)
           
        if by == caixas_meia_altura :
            vx_axis1.append(vx_win[i])
            vy_axis1.append(vy_win[i])
            xxpyy,xxmyy,xy=matrix_three_component(texture_win[i])
            texture_axis_a_axis1.append(xxpyy)
            texture_axis_b_axis1.append(xxmyy)
            texture_ang_elipse_axis1.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(B_win[i])
            B_axis_a_axis1.append(xxpyy)
            B_axis_b_axis1.append(xxmyy)
            B_ang_elipse_axis1.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(T_win[i])
            T_axis_a_axis1.append(xxpyy)
            T_axis_b_axis1.append(xxmyy)
            T_ang_elipse_axis1.append(xy)

        if by == caixas_quarto_altura :
            vx_axis2.append(vx_win[i])
            vy_axis2.append(vy_win[i])
            xxpyy,xxmyy,xy=matrix_three_component(texture_win[i])
            texture_axis_a_axis2.append(xxpyy)
            texture_axis_b_axis2.append(xxmyy)
            texture_ang_elipse_axis2.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(B_win[i])
            B_axis_a_axis2.append(xxpyy)
            B_axis_b_axis2.append(xxmyy)
            B_ang_elipse_axis2.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(T_win[i])
            T_axis_a_axis2.append(xxpyy)
            T_axis_b_axis2.append(xxmyy)
            T_ang_elipse_axis2.append(xy)
        if by == 3*caixas_quarto_altura :
            vx_axis6.append(vx_win[i])
            vy_axis6.append(vy_win[i])
            xxpyy,xxmyy,xy=matrix_three_component(texture_win[i])
            texture_axis_a_axis6.append(xxpyy)
            texture_axis_b_axis6.append(xxmyy)
            texture_ang_elipse_axis6.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(B_win[i])
            B_axis_a_axis6.append(xxpyy)
            B_axis_b_axis6.append(xxmyy)
            B_ang_elipse_axis6.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(T_win[i])
            T_axis_a_axis6.append(xxpyy)
            T_axis_b_axis6.append(xxmyy)
            T_ang_elipse_axis6.append(xy)
        if bx == caixas_meia_largura :
            vx_axis3.append(vx_win[i])
            vy_axis3.append(vy_win[i])
            xxpyy,xxmyy,xy=matrix_three_component(texture_win[i])
            texture_axis_a_axis3.append(xxpyy)
            texture_axis_b_axis3.append(xxmyy)
            texture_ang_elipse_axis3.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(B_win[i])
            B_axis_a_axis3.append(xxpyy)
            B_axis_b_axis3.append(xxmyy)
            B_ang_elipse_axis3.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(T_win[i])
            T_axis_a_axis3.append(xxpyy)
            T_axis_b_axis3.append(xxmyy)
            T_ang_elipse_axis3.append(xy)

        if bx == caixas_meia_largura-caixas_quarto_largura :
            vx_axis4.append(vx_win[i])
            vy_axis4.append(vy_win[i])
            xxpyy,xxmyy,xy=matrix_three_component(texture_win[i])
            texture_axis_a_axis4.append(xxpyy)
            texture_axis_b_axis4.append(xxmyy)
            texture_ang_elipse_axis4.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(B_win[i])
            B_axis_a_axis4.append(xxpyy)
            B_axis_b_axis4.append(xxmyy)
            B_ang_elipse_axis4.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(T_win[i])
            T_axis_a_axis4.append(xxpyy)
            T_axis_b_axis4.append(xxmyy)
            T_ang_elipse_axis4.append(xy)
        if bx == caixas_meia_largura+caixas_quarto_largura :
            vx_axis5.append(vx_win[i])
            vy_axis5.append(vy_win[i])
            xxpyy,xxmyy,xy=matrix_three_component(texture_win[i])
            texture_axis_a_axis5.append(xxpyy)
            texture_axis_b_axis5.append(xxmyy)
            texture_ang_elipse_axis5.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(B_win[i])
            B_axis_a_axis5.append(xxpyy)
            B_axis_b_axis5.append(xxmyy)
            B_ang_elipse_axis5.append(xy)
            xxpyy,xxmyy,xy=matrix_three_component(T_win[i])
            T_axis_a_axis5.append(xxpyy)
            T_axis_b_axis5.append(xxmyy)
            T_ang_elipse_axis5.append(xy)

    for i in range(box_per_line_x):

        file_axis1.write("%d %f %f %f %f %f %f %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis1[i], vy_axis1[i], texture_axis_a_axis1[i], texture_axis_b_axis1[i], texture_ang_elipse_axis1[i], B_axis_a_axis1[i], B_axis_b_axis1[i], B_ang_elipse_axis1[i], T_axis_a_axis1[i], T_axis_b_axis1[i], T_ang_elipse_axis1[i]))

        file_axis2.write("%d %f %f %f %f %f %f %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis2[i], vy_axis2[i], texture_axis_a_axis2[i], texture_axis_b_axis2[i], texture_ang_elipse_axis2[i], B_axis_a_axis2[i], B_axis_b_axis2[i], B_ang_elipse_axis2[i], T_axis_a_axis2[i], T_axis_b_axis2[i], T_ang_elipse_axis2[i]))

        file_axis6.write("%d %f %f %f %f %f %f %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis6[i], vy_axis6[i], texture_axis_a_axis6[i], texture_axis_b_axis6[i], texture_ang_elipse_axis6[i], B_axis_a_axis6[i], B_axis_b_axis6[i], B_ang_elipse_axis6[i], T_axis_a_axis6[i], T_axis_b_axis6[i], T_ang_elipse_axis6[i]))


    for i in range(box_per_column_y):
        file_axis3.write("%d %f %f %f %f %f %f %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis3[i], vy_axis3[i], texture_axis_a_axis3[i], texture_axis_b_axis3[i], texture_ang_elipse_axis3[i], B_axis_a_axis3[i], B_axis_b_axis3[i], B_ang_elipse_axis3[i], T_axis_a_axis3[i], T_axis_b_axis3[i], T_ang_elipse_axis3[i]))

        file_axis4.write("%d %f %f %f %f %f %f %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis4[i], vy_axis4[i], texture_axis_a_axis4[i], texture_axis_b_axis4[i], texture_ang_elipse_axis4[i], B_axis_a_axis4[i], B_axis_b_axis4[i], B_ang_elipse_axis4[i], T_axis_a_axis4[i], T_axis_b_axis4[i], T_ang_elipse_axis4[i]))
        
        file_axis5.write("%d %f %f  %f %f %f %f %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis5[i], vy_axis5[i], texture_axis_a_axis5[i], texture_axis_b_axis5[i], texture_ang_elipse_axis5[i], B_axis_a_axis5[i], B_axis_b_axis5[i], B_ang_elipse_axis5[i], T_axis_a_axis5[i], T_axis_b_axis5[i], T_ang_elipse_axis5[i]))

    plt.subplot(211)
    plt.ylabel('Vx')
    plt.plot(vx_axis1,'k',vx_axis2,'r',vx_axis6,'g')
    plt.subplot(212)
    plt.ylabel('Vy')
    plt.xlabel('X')
    plt.plot(vy_axis1,'k',vy_axis2,'r',vy_axis6,'g')
    plt.savefig(path+"/six-axis-velocity-field-X-direction.png",bbox_inches="tight")
    plt.close()

    plt.subplot(211)
    plt.ylabel('Vx')
    plt.plot(vx_axis3,'k',vx_axis4,'r',vx_axis5,'g')
    plt.subplot(212)
    plt.ylabel('Vy')
    plt.xlabel('Y')
    plt.plot(vy_axis3,'k',vy_axis4,'r',vy_axis5,'g')
    plt.savefig(path+"/six-axis-velocity-field-Y-direction.png",bbox_inches="tight")
    plt.close()

    

    
    plt.subplot(211)
    plt.title('(XX+YY)/2')
    plt.plot(texture_axis_a_axis1,'k',texture_axis_a_axis2,'r',texture_axis_a_axis6,'g')
    plt.subplot(212)
    plt.plot(texture_axis_a_axis3,'k',texture_axis_a_axis4,'r',texture_axis_a_axis5,'g')
    plt.savefig(path+"/six-axis-texture-xxpyy.png",bbox_inches="tight")
    plt.close()

    plt.subplot(211)
    plt.title('(XX-YY)/2')
    plt.plot(texture_axis_b_axis1,'k',texture_axis_b_axis2,'r',texture_axis_b_axis6,'g')
    plt.subplot(212)
    plt.plot(texture_axis_b_axis3,'k',texture_axis_b_axis4,'r',texture_axis_b_axis5,'g')
    plt.savefig(path+"/six-axis-texture-xxmyy.png",bbox_inches="tight")
    plt.close()

    plt.subplot(211)
    plt.title('XY')
    plt.plot(texture_ang_elipse_axis1,'k',texture_ang_elipse_axis2,'r',texture_ang_elipse_axis6,'g')
    plt.subplot(212)
    plt.plot(texture_ang_elipse_axis3,'k',texture_ang_elipse_axis4,'r',texture_ang_elipse_axis5,'g')
    plt.savefig(path+"/six-axis-texture-xy.png",bbox_inches="tight")
    plt.close()
    #    plt.show()

def imag_count(system_type) :
    counter = 0
    max_number_particles = 0
    part_counter=0
    print "Counting images... wait... it may take 5s to count 1000 images on an I7 \n"
    if system_type == 'superboids' :
        while 1 :
            line          = file_arq_neigh_in.readline()
            if not line :
                break #EOF
            if line.replace( '\r', '' ) == '\n' : #identifies blank lines
                counter             += 1
                max_number_particles = max(max_number_particles,part_counter)
                part_counter = 0

                while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                    line = file_arq_neigh_in.readline()
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
            line = file_arq_data_in.readline()
            if not line:
                break
            line_splitted = line.split()
            if line_splitted[0] == 'x' :
                counter += 1

    if system_type == 'vicsek-gregoire' :
        name_arq_data_in = "%s/data/posicoes.dat"% (system_type)
        file_arq_data_in       = open(name_arq_data_in)

        max_number_particles=0
        while 1:
            line = file_arq_data_in.readline()           
            if not line:
                break #EOF
            line_splitted = line.split()
            counter += 1
            n=int(line_splitted[1])
            max_number_particles=max(max_number_particles,n)
            for i in range(n):
                file_arq_data_in.readline()
        file_arq_data_in.close()
    if system_type == 'voronoi' :
        max_number_particles=0
        os.system("ls voronoi/cell_*.dat > files.dat")
        file_names = open("files.dat",'r')
        for line in file_names:
            counter+=1
            part_counter =0 
            for word in line.split():
                name_arq_data_in = word
                file_arq_data_in       = open(name_arq_data_in)
                while 1:
                    line = file_arq_data_in.readline()           
                    if not line:
                        break #EOF
                    line_splitted = line.split()
                    if line_splitted[1] == '1' :
                        part_counter += 1
                max_number_particles=max(max_number_particles,part_counter)
                file_arq_data_in.close()

    if system_type == 'potts' :
        max_number_particles=0
        while 1:
            line = file_arq_data_in.readline()           
            if not line:
                break #EOF
            line_splitted = line.split()
            counter += 1
            n=int(line_splitted[1])
            max_number_particles=max(max_number_particles,n)
            for i in range(n):
                file_arq_data_in.readline()
            
    print "Counted", counter-1, "images.\n"
#    print "Type initial and final image number you want to analyse (min=1, max=",counter-1,") - Use spaces to separate the two numbers"
    print "Counted", max_number_particles, "as max number of particles."
    return max_number_particles

################## Here starts the main program ###############

#Opening input parameter file

file_input_parameter = open("parameter.in")
line_splitted        = file_input_parameter.readline().split()
system_type          = line_splitted[1]
path    = 'output/'+system_type

#Creating the directory structure for output
os.system('mkdir -p %s' % path)
#os.system('cp axis12.par %s' % path)
#os.system('cp axis345.par %s' % path)

#velocity_density gnuplot file
vid_veloc_dens = open("%s/video_velocity_density.gnu"%path,"w")

#Texture, B and  elipsis script header

vid_def = open("%s/video_deformation.gnu"%path,"w")
vid_def.write("unset key \n")
vid_B = open("%s/video_B.gnu"%path,"w")
vid_B.write("unset key \n")
vid_T = open("%s/video_T.gnu"%path,"w")
vid_T.write("unset key \n")

#Opening time averages files
dens_win_file_name    = "density-win.dat"
dens_win              = open(path+'/'+dens_win_file_name,"w")
vel_win_file_name     = "velocity-win.dat"
vel_win               = open(path+'/'+vel_win_file_name,"w")
vel_fluct_win_file_name     = "velocity-fluct-win.dat"
vel_fluct_win               = open(path+'/'+vel_fluct_win_file_name,"w")
def_win               = open("%s/deformation-win.dat"%path,"w")
texture_win_file_name = "texture-win.dat"
texture_win_file      = open(path+'/'+texture_win_file_name,"w")
B_win_file_name = "B-win.dat"
B_win_file      = open(path+'/'+B_win_file_name,"w")
T_win_file_name = "T-win.dat"
T_win_file      = open(path+'/'+T_win_file_name,"w")

# Opening six axis analysis files

file_axis1 = open("%s/axis1.dat"%path,"w")
file_axis2 = open("%s/axis2.dat"%path,"w")
file_axis4 = open("%s/axis4.dat"%path,"w")
file_axis3 = open("%s/axis3.dat"%path,"w")
file_axis5 = open("%s/axis5.dat"%path,"w")
file_axis6 = open("%s/axis6.dat"%path,"w")

if system_type == 'experiment':
    arq_in      = "%s/%s"%(line_splitted[1], line_splitted[2])
    print "You analise an", system_type, "system, reading data from file:\n", arq_in
    window_size = int(line_splitted[3])
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
                axis_b_tot[counter]     += axis_b[counter]
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
            velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0)
            #Function call to write deformation elipse gnu script
            deformation_elipsis_script(x, y, axis_a, axis_b, ang_elipse, system_type)
        elif image > image_f:
            break
        

if system_type == "superboids":
    line_splitted = file_input_parameter.readline().split()
    name_arq_header_in = "%s/%s.dat"%(system_type,line_splitted[1])
    name_arq_data_in   = "%s/%s_plainprint.dat"%(system_type,line_splitted[1])
    name_arq_neigh_in  = "%s/%s_neighbors.dat"%(system_type,line_splitted[1])
    box_size      = float(file_input_parameter.readline().split()[1])
    window_size   = int(file_input_parameter.readline().split()[1])
    max_dist      = float(file_input_parameter.readline().split()[1])
    
    print "\nYou analise a", system_type, "system, data is read from files:\n", name_arq_header_in," (header)\n", name_arq_data_in," (data)\n", name_arq_neigh_in," (neighbors)"
    file_arq_header_in    = open(name_arq_header_in)
    file_arq_data_in  = open(name_arq_data_in)
    file_arq_neigh_in = open(name_arq_neigh_in)
    line_splitted = file_input_parameter.readline().split()
    x0 = int(line_splitted[1])
    line_splitted = file_input_parameter.readline().split()
    xf = int(line_splitted[1])
    file_arq_neigh_in = open(name_arq_neigh_in)
    max_number_particles=imag_count(system_type)
    file_arq_neigh_in.close()
    line_splitted = sys.stdin.readline().split() #lendo do teclado imagens final e inicial depois da analise do numero de imagens (funcao imag_count())
    image_0       = int(line_splitted[0])
    image_f       = int(line_splitted[1])
    v0            = 0.007
    part=list(particle(i) for i in range(max_number_particles))

    # Reading superboids parameter file

    while 1 :
        line = file_arq_header_in.readline()
        if not line : break # EOF
        if line.replace( '\r' , '' ) == "\n" :
            continue
        else :
            line_splitted = line.split()
            if line_splitted[0] == '#' and line_splitted[1]=='RECTANGLE:':
                line_splitted = file_arq_header_in.readline().split()
                Lx, Ly        = int(float(line_splitted[1])), int(float(line_splitted[2]))
                Lx            = 240 # corrigindo por enquanto o erro no arquivo de entrada. REVISAR!!!
            if line_splitted[1] == 'Radial' and line_splitted[2]=='R_Eq:':
                line_splitted = file_arq_header_in.readline().split()
            if line_splitted[0] == "#radius:":
                R_OBST        = int(line_splitted[1])
                X_OBST        = float(line_splitted[3])
                Y_OBST        = float(line_splitted[4])
    delta_x =int((xf+x0)*R_OBST/box_size)*box_size
    x0 = X_OBST-x0*R_OBST
    xf = x0+delta_x
    y0, yf = box_size*int(-Ly/(2*box_size)), box_size*int(Ly/(2*box_size))
    box_per_line_x, box_per_column_y = int((delta_x)/box_size), int((yf-y0)/box_size)

    if x0 < -Lx/2 or xf > Lx/2 :
        print "Warning: Reseting limits to -Lx/2, Lx/2"
        x0 = -Lx/2
        xf = Lx/2
    #defining all matrices
    box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, vx2_tot, vy2_tot,\
        density_tot, vx_win, vy_win,vx2_win, vy2_win, density_win, texture_box, B_box, \
        T_box, texture_tot, B_tot, T_tot, texture_win, B_win, T_win=\
        box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)

    
    #Reading superboids plainprint data file
    fat_boids_counter = 0
    image             = 0
    count_events      = 0
    points            = []
    index_particle    = []
    
    #local defintions to put diagonalized matrices values
    axis_a_texture     = list(0. for i in range(box_total))
    axis_b_texture     = list(0. for i in range(box_total))
    ang_elipse_texture = list(0. for i in range(box_total))
    axis_a_B           = list(0. for i in range(box_total))
    axis_b_B           = list(0. for i in range(box_total))
    ang_elipse_B       = list(0. for i in range(box_total))
    axis_a_T           = list(0. for i in range(box_total))
    axis_b_T           = list(0. for i in range(box_total))
    ang_elipse_T       = list(0. for i in range(box_total))



    while 1 :
        line = file_arq_data_in.readline()
        if not line :
            vid_veloc_dens.write("pause -1 \n")
            break #EOF
        if not line.replace( '\r', '' ) == '\n' : #identifies not blank lines
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

        else:
            image      += 1
            if image <= image_0 :
                print "Skippin image:",image 
            elif image <= image_f :
                print "Image number",image,". Number of super particles =", fat_boids_counter
                count_events     += 1
                # Calculus of textures, B and T
                number_particles = len(points)
                if number_particles > 0:
                    points=np.array(points)
                    list_neighbors=delaunay(points,max_dist)
                    map_focus_region_to_part(points,list_neighbors,index_particle)
                    map(lambda i:i.texture(), part)
                    if  count_events > 1 :
                        map(lambda i:i.calc_B_and_T(), part)
                    for i in index_particle:
                        xx  = int((part[i].r[0]-x0) / box_size)
                        yy  = int((part[i].r[1]-y0) / box_size) * box_per_line_x
                        box = xx+yy
#                        density_now[box] += 1.0
                        texture_box[box]+=part[i].M
                        if count_events > 1 :
                            B_box[box] += part[i].B
                            T_box[box] += part[i].T

                    map(lambda i:i.copy_to_old(), part)


                for box in range(box_total):
                    if density_now[box] > 0 :

                        vx_now[box]                 = vx_now[box] / density_now[box]
		        vy_now[box]                 = vy_now[box] / density_now[box]
                        texture_box[box]            = texture_box[box] / density_now[box]
                        ax_a,ax_b,angle             = axis_angle(texture_box[box])
                        axis_a_texture[box]     = ax_a
                        axis_b_texture[box]     = ax_b
                        ang_elipse_texture[box] = angle
                        if count_events > 1 :
                            B_box[box]            = B_box[box] / density_now[box]
                            ax_a,ax_b,angle       = axis_angle(B_box[box])
                            axis_a_B[box]     = ax_a
                            axis_b_B[box]     = ax_b
                            ang_elipse_B[box] = angle
                            T_box[box]            = T_box[box] / density_now[box]
                            ax_a,ax_b,angle       = axis_angle(T_box[box])
                            axis_a_T[box]     = ax_a
                            axis_b_T[box]     = ax_b
                            ang_elipse_T[box] = angle

                #Function call to write velocity-density gnu script
                velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0)
                #Function call to write texture elipsis gnu script
                texture_elipsis_script_simu(box_per_line_x, box_total, axis_a_texture, axis_b_texture,\
                                                ang_elipse_texture, image-image_0,points,x0,y0,box_size)
                if count_events > 1 :
                    B_elipsis_script_simu(box_per_line_x, box_total, axis_a_B, axis_b_B,\
                                                ang_elipse_B, image-image_0,points,x0,y0,box_size)
                    T_elipsis_script_simu(box_per_line_x, box_total, axis_a_T, axis_b_T,\
                                                ang_elipse_T, image-image_0,points,x0,y0,box_size)


            #Summing each box at different times
            #if image > image_0 and image <= image_f:
                for box in range(box_total) :
                    density_tot[box] += density_now[box]
                    vx_tot[box]      += vx_now[box]
                    vy_tot[box]      += vy_now[box]
                    texture_tot[box] += texture_box[box]
                    B_tot[box] += B_box[box]
                    T_tot[box] += T_box[box]
                vx_now      = list(0. for i in range(box_total))
                vy_now      = list(0. for i in range(box_total))
                density_now = list(0  for i in range(box_total))
                texture_box = list(np.zeros((2,2)) for i in range(box_total))
                B_box = list(np.zeros((2,2)) for i in range(box_total))
                T_box = list(np.zeros((2,2)) for i in range(box_total))
                points            = []
                index_particle    = []

            else:
                break
            fat_boids_counter = 0
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = file_arq_data_in.readline()

    image_counter = image_f - image_0
    
if system_type == "szabo-boids":
    name_arq_data_in = "%s/%s.dat"% (line_splitted[0], line_splitted[1])
    print "\nYou analise a", line_splitted[0], "system, data is read from files:\n", arq_data_in
    file_arq_data_in           = open(name_arq_data_in)
    window_size   = int(line_splitted[2])
    imag_count(system_type)
    file_arq_data_in.close()
    file_arq_data_in           = open(name_arq_data_in)
    line_splitted = sys.stdin.readline().split()
    image_0       = int(line_splitted[0])
    image_f       = int(line_splitted[1])
    image_counter = image_f-image_0
    image         = 0
    line_counter  = 0
    count_events  = 0
    v0            = 0.1

    #local defintions to put diagonalized matrices values
    axis_a_texture     = list(0. for i in range(box_total))
    axis_b_texture     = list(0. for i in range(box_total))
    ang_elipse_texture = list(0. for i in range(box_total))
    axis_a_B           = list(0. for i in range(box_total))
    axis_b_B           = list(0. for i in range(box_total))
    ang_elipse_B       = list(0. for i in range(box_total))
    axis_a_T           = list(0. for i in range(box_total))
    axis_b_T           = list(0. for i in range(box_total))
    ang_elipse_T       = list(0. for i in range(box_total))

    #Reading szabo-boids  data file
    while 1 :
        line = file_arq_data_in.readline()
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

                box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, vx2_tot, vy2_tot,\
                density_tot, vx_win, vy_win,vx2_win, vy2_win, density_win, texture_box, B_box, \
                T_box, texture_tot, B_tot, T_tot, texture_win, B_win, T_win=\
                box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)


            if line_splitted[0] == "Radius:" : R_OBST = int(line_splitted[1])
            if line_splitted[0] == "Obst_position:" :
                X_OBST = int(line_splitted[1])
                Y_OBST = int(line_splitted[2])
        if line_counter > 6 :
            if image <= image_0 :
                print "Skipping image:",image 
                while line_splitted[0] != 'x' :
                    line = file_arq_data_in.readline()
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
                        line = file_arq_data_in.readline()
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
                velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0)
                #Summing each box at different times
                for box in range(box_total) :
                    density_tot[box] += density_now[box]
                    vx_tot[box] += vx_now[box]
                    vy_tot[box] += vy_now[box]

            else:
                break

if system_type == "vicsek-gregoire":

    window_size, time_0, time_f, obstacle, x0, xf, filename, box_mag = read_param(file_input_parameter)
    file_input_parameter.close()
    aux           = "%s/include/%s"% (system_type, filename)
    file_par_simu = open(aux)
    Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size, Delta_t, v0, max_dist = read_param_vic_greg(file_par_simu)
    file_par_simu.close()
    delta_x =int((xf+x0)*R_OBST/box_size)*box_size
    x0 = X_OBST-x0*R_OBST
    xf = x0+delta_x
    y0 = 0.
    yf =  box_size*int(Ly/box_size)
    box_per_line_x, box_per_column_y = int((delta_x)/box_size), int((yf-y0)/box_size)
    
    if x0 < 0. or xf > Lx :
        print "Warning: Reseting limits to 0, Lx"
        x0 = 0.
        xf = Lx

    #definindo as caixas e as matrizes

    
    box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, vx2_tot, vy2_tot,\
                density_tot, vx_win, vy_win,vx2_win, vy2_win, density_win, texture_box, B_box, \
                T_box, texture_tot, B_tot, T_tot, texture_win, B_win, T_win=\
               box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)

    #arquivo de posicoes e velocidades

    name_arq_data_in = "%s/data/posicoes.dat"% (system_type)
    print "\nYou analise a", system_type, "system, data is read from files:\n", name_arq_data_in
    max_number_particles   = imag_count(system_type) #conta o numero de imagens
    file_arq_data_in       = open(name_arq_data_in)  #reabre o arquivo para leituras das posicoes e vel.
#    line_splitted = sys.stdin.readline().split() #le da linha de comando o intervalo de imagens desejado
    image_0       = int(time_0/Delta_t)
    image_f       = int(time_f/Delta_t)
    image_counter = image_f - image_0
    image         = 1
    line_counter  = 0
    count_events  = 0
    part=list(particle(i) for i in range(max_number_particles))
    #local defintions to put diagonalized matrices values
    axis_a_texture     = list(0. for i in range(box_total))
    axis_b_texture     = list(0. for i in range(box_total))
    ang_elipse_texture = list(0. for i in range(box_total))
    axis_a_B           = list(0. for i in range(box_total))
    axis_b_B           = list(0. for i in range(box_total))
    ang_elipse_B       = list(0. for i in range(box_total))
    axis_a_T           = list(0. for i in range(box_total))
    axis_b_T           = list(0. for i in range(box_total))
    ang_elipse_T       = list(0. for i in range(box_total))

    
    while 1 :
        line   = file_arq_data_in.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        nlines        = int(line_splitted[1])
        if image <= image_0 :
            print "Skipping image:",image
            for i in range(nlines) :
                file_arq_data_in.readline()
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
                line_splitted    = file_arq_data_in.readline().split()
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
                list_neighbors=delaunay(points,max_dist)
                map_focus_region_to_part(points,list_neighbors,index_particle)
                map(lambda i:i.texture(), part)
                for i in index_particle:
                    xx  = int((part[i].r[0]-x0) / box_size)
                    yy  = int((part[i].r[1]-y0) / box_size) * box_per_line_x
                    box = xx+yy
                    texture_box[box]+=part[i].M
                    if  count_events > 1 :
                        B_box[box] += part[i].B
                        T_box[box] += part[i].T
                       
                map(lambda i:i.copy_to_old(), part)

                #Calculate the average velocity over boxes
            for box in range(box_total):
                if density_now[box] > 0 :
                    vx_now[box] = vx_now[box] / density_now[box]
	            vy_now[box] = vy_now[box] / density_now[box]
                    texture_box[box]            = texture_box[box] / density_now[box]
                    ax_a,ax_b,angle             = axis_angle(texture_box[box])
                    axis_a_texture[box]     = ax_a
                    axis_b_texture[box]     = ax_b
                    ang_elipse_texture[box] = angle
                    if count_events > 1 :
                        B_box[box]            = B_box[box] / density_now[box]
                        ax_a,ax_b,angle       = axis_angle(B_box[box])
                        axis_a_B[box]     = ax_a
                        axis_b_B[box]     = ax_b
                        ang_elipse_B[box] = angle
                        T_box[box]            = T_box[box] / density_now[box]
                        ax_a,ax_b,angle       = axis_angle(T_box[box])
                        axis_a_T[box]     = ax_a
                        axis_b_T[box]     = ax_b
                        ang_elipse_T[box] = angle

            #Function call to write velocity-density gnu script
            velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0)


            #Function call to write texture elipsis gnu script
            texture_elipsis_script_simu(box_per_line_x, box_total, axis_a_texture, axis_b_texture,\
                                    ang_elipse_texture, image-image_0,points,x0,y0,box_size)
            if count_events > 1 :
                B_elipsis_script_simu(box_per_line_x, box_total, axis_a_B, axis_b_B,\
                                      ang_elipse_B, image-image_0,points,x0,y0,box_size)
                T_elipsis_script_simu(box_per_line_x, box_total, axis_a_T, axis_b_T,\
                                      ang_elipse_T, image-image_0,points,x0,y0,box_size)
                        

            #Summing each box at different times
            for box in range(box_total) :
                density_tot[box] += density_now[box]
                vx_tot[box]      += vx_now[box]
                vy_tot[box]      += vy_now[box]
                texture_tot[box] += texture_box[box]
                B_tot[box] += B_box[box]
                T_tot[box] += T_box[box]
            #reseting matrices of instaneous measures                
            vx_now      = list(0. for i in range(box_total))
            vy_now      = list(0. for i in range(box_total))
            density_now = list(0  for i in range(box_total))
            texture_box = list(np.zeros((2,2)) for i in range(box_total))
            B_box = list(np.zeros((2,2)) for i in range(box_total))
            T_box = list(np.zeros((2,2)) for i in range(box_total))
            points            = []
            index_particle    = []

        else:
                break

if system_type == "potts":

    window_size, time_0, time_f, obstacle, x0, xf, filename, box_mag = read_param(file_input_parameter)
    file_input_parameter.close()
    aux           = "%s/Simulation/%s"% (system_type, filename)
    file_par_simu = open(aux)
    Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size, max_dist, Delta_t, v0 = read_param_potts(file_par_simu)
    file_par_simu.close()
    delta_x =int((xf+x0)*R_OBST/box_size)*box_size
    x0 = X_OBST-x0*R_OBST
    xf = x0+delta_x
    y0 = 0.
    yf =  box_size*int(Ly/box_size)
    box_per_line_x, box_per_column_y = int((delta_x)/box_size), int((yf-y0)/box_size)

    if x0 < 0. or xf > Lx :
        print "Warning: Reseting limits to 0, Lx"
        x0 = 0.
        xf = Lx

    #definindo as caixas e as matrizes
    box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, vx2_tot, vy2_tot,\
        density_tot, vx_win, vy_win,vx2_win, vy2_win, density_win, texture_box, B_box, \
T_box, texture_tot, B_tot, T_tot, texture_win, B_win, T_win= \
        box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)

    #arquivo de posicoes e velocidades

    name_arq_data_in = "%s/posicoes.dat"% (system_type)
    print "\nYou analise a", system_type, "system, data is read from files:\n", name_arq_data_in
    file_arq_data_in       = open(name_arq_data_in)
    max_number_particles   = imag_count(system_type) #conta o numero de imagens
    print max_number_particles
    file_arq_data_in.close()
    file_arq_data_in       = open(name_arq_data_in)  #reabre o arquivo para leituras das posicoes e vel.
#    line_splitted = sys.stdin.readline().split() #le da linha de comando o intervalo de imagens desejado
    image_0       = int(time_0/Delta_t)
    image_f       = int(time_f/Delta_t)
    image_counter = image_f - image_0
    image         = 1
    line_counter  = 0
    count_events  = 0
    part=list(particle(i) for i in range(max_number_particles))
    #local defintions to put diagonalized matrices values
    axis_a_texture     = list(0. for i in range(box_total))
    axis_b_texture     = list(0. for i in range(box_total))
    ang_elipse_texture = list(0. for i in range(box_total))
    axis_a_B           = list(0. for i in range(box_total))
    axis_b_B           = list(0. for i in range(box_total))
    ang_elipse_B       = list(0. for i in range(box_total))
    axis_a_T           = list(0. for i in range(box_total))
    axis_b_T           = list(0. for i in range(box_total))
    ang_elipse_T       = list(0. for i in range(box_total))

    
    while 1 :
        line   = file_arq_data_in.readline()
        if not line:
            break #EOF
        line_splitted = line.split()
        nlines        = int(line_splitted[1])
        if image <= image_0 :
            print "Skipping image:",image
            for i in range(nlines) :
                file_arq_data_in.readline()
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
                line_splitted    = file_arq_data_in.readline().split()
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
            if points == [] :
                print " "
                print "Your first image has no points in the analysis area!"
                print " "
                exit()
            # Calculus of textures, B and T##################
            number_particles = len(points)
            if number_particles > 0:
                points=np.array(points)
                list_neighbors=delaunay(points,max_dist)
                map_focus_region_to_part(points,list_neighbors,index_particle)

                #calculate textures, B and T
                map(lambda i:i.texture(), part)
                if  count_events > 1 :
                    map(lambda i:i.calc_B_and_T(), part)

                for i in index_particle:
                    xx  = int((part[i].r[0]-x0) / box_size)
                    yy  = int((part[i].r[1]-y0) / box_size) * box_per_line_x
                    box = xx+yy
                    texture_box[box]+=part[i].M
                if  count_events > 1 :
                    B_box[box] += part[i].B
                    T_box[box] += part[i].T

                map(lambda i:i.copy_to_old(), part)

                #Calculate the averages over boxes
            for box in range(box_total):
                if density_now[box] > 0 :
                    vx_now[box] = vx_now[box] / density_now[box]
	            vy_now[box] = vy_now[box] / density_now[box]
                    texture_box[box]            = texture_box[box] / density_now[box]
                    ax_a,ax_b,angle             = axis_angle(texture_box[box])
                    axis_a_texture[box]     = ax_a
                    axis_b_texture[box]     = ax_b
                    ang_elipse_texture[box] = angle
                    if count_events > 1 :
                        B_box[box]            = B_box[box] / density_now[box]
                        ax_a,ax_b,angle       = axis_angle(B_box[box])
                        
                        axis_a_B[box]     = ax_a
                        axis_b_B[box]     = ax_b
                        ang_elipse_B[box] = angle
                        T_box[box]            = T_box[box] / density_now[box]
                        ax_a,ax_b,angle       = axis_angle(T_box[box])
                        axis_a_T[box]     = ax_a
                        axis_b_T[box]     = ax_b
                        ang_elipse_T[box] = angle

            #Function call to write velocity-density gnu script
            velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0)
            
            #Function call to write texture elipsis gnu script
            texture_elipsis_script_simu(box_per_line_x, box_total, axis_a_texture, axis_b_texture,\
                                    ang_elipse_texture, image-image_0,points,x0,y0,box_size)
            if count_events > 1 :
                B_elipsis_script_simu(box_per_line_x, box_total, axis_a_B, axis_b_B,\
                                      ang_elipse_B, image-image_0,points,x0,y0,box_size)
                T_elipsis_script_simu(box_per_line_x, box_total, axis_a_T, axis_b_T,\
                                      ang_elipse_T, image-image_0,points,x0,y0,box_size)
                        

            #Summing each box at different times
            for box in range(box_total) :
                density_tot[box] += density_now[box]
                vx_tot[box]      += vx_now[box]
                vy_tot[box]      += vy_now[box]
                texture_tot[box] += texture_box[box]
                B_tot[box] += B_box[box]
                T_tot[box] += T_box[box]
            vx_now      = list(0. for i in range(box_total))
            vy_now      = list(0. for i in range(box_total))
            density_now = list(0  for i in range(box_total))
            texture_box = list(np.zeros((2,2)) for i in range(box_total))
            B_box = list(np.zeros((2,2)) for i in range(box_total))
            T_box = list(np.zeros((2,2)) for i in range(box_total))
            points            = []
            index_particle    = []

        else:
                break

if system_type == "voronoi":
    # voronoi model works with obstacle at 0.0 , 0.0 and system goes from - Lx/2.0 to Lx/2.0 but the articles can go even further
    window_size, time_0, time_f, obstacle, x0, xf, filename, box_mag = read_param(file_input_parameter)
    file_input_parameter.close()
    Lx, Ly, R_OBST, X_OBST, Y_OBST, box_size, Delta_t, v0, max_dist = read_param_voronoi(filename)
    box_size = box_size*box_mag
    max_dist = max_dist*box_mag
    delta_x =int((xf+x0)*R_OBST/box_size)*box_size
    x0 = X_OBST-x0*R_OBST
    xf = x0+delta_x
    y0 = -box_size*int(Ly/box_size)/2
    yf =  box_size*int(Ly/box_size)/2
    box_per_line_x, box_per_column_y = int((delta_x)/box_size), int((yf-y0)/box_size)

    if x0 < -Lx/2.0 or xf > Lx/2.0 :
        print "Warning: Reseting limits to 0, Lx"
        x0 = -Lx/2.0
        xf =  Lx/2.0

    #definindo as caixas e as matrizes

    box_total, ratio, vx_now, vy_now, density_now, vx_tot, vy_tot, vx2_tot, vy2_tot,\
        density_tot, vx_win, vy_win,vx2_win, vy2_win, density_win, texture_box, B_box, \
T_box, texture_tot, B_tot, T_tot, texture_win, B_win, T_win= \
        box_variables_definition_simu(box_per_column_y, box_per_line_x, x0, y0, xf, yf)
    image_0       = int(time_0/Delta_t)
    image_f       = int(time_f/Delta_t)
    image_counter = image_f - image_0
    image         = 1
    line_counter  = 0
    count_events  = 0
    
    axis_a_texture     = list(0. for i in range(box_total))
    axis_b_texture     = list(0. for i in range(box_total))
    ang_elipse_texture = list(0. for i in range(box_total))
    axis_a_B           = list(0. for i in range(box_total))
    axis_b_B           = list(0. for i in range(box_total))
    ang_elipse_B       = list(0. for i in range(box_total))
    axis_a_T           = list(0. for i in range(box_total))
    axis_b_T           = list(0. for i in range(box_total))
    ang_elipse_T       = list(0. for i in range(box_total))

    print "\nYou analise a", system_type, "system \n"
    max_number_particles = imag_count(system_type) #conta o numero de imagens e o numero maximo de particulas
    nlines        = max_number_particles
    #arquivo de posicoes e velocidades
    os.system("ls voronoi/cell_*.dat > files.dat")
    file_names = open("files.dat",'r')
    # each word in this files.dat is a simulation snapshot datafile, so we need to open file by file and read the information
    for datafile in file_names:
        name_arq_data_in = datafile.split()[0]
        file_arq_data_in = open(name_arq_data_in)  #reabre o arquivo para leituras das posicoes e vel.
        if image <= image_0 :
            print "Skipping image:",image
            image += 1
        elif image <= image_f :
            print "Reading image:",image
            part=list(particle(i) for i in range(max_number_particles))
            vx_now = list(0. for i in range(box_total))
            vy_now = list(0. for i in range(box_total))
            density_now = list(0. for i in range(box_total))
            texture_box = list(np.zeros((2,2)) for i in range(box_total))
            points = []
            index_particle = []
            count_particle = 0
            line   = file_arq_data_in.readline() # first line is not useful for us


            while 1 :  #This scans arq_data_in line by line up to the end
                line   = file_arq_data_in.readline()
                if not line:
                    image        += 1
                    count_events += 1
                    break #EOF
                line_splitted = line.split()
                particle_type, x, y = int(line_splitted[1]), float(line_splitted[2]), float(line_splitted[3])
                if x > x0 and x < xf and y > y0 and y < yf and particle_type == 1: 
                    xx  = int((x-x0) / box_size)
                    yy  = int((y-y0) / box_size) * box_per_line_x
                    box = xx+yy
                    vx_now[box]      += float(line_splitted[5])
                    vy_now[box]      += float(line_splitted[6])
                    density_now[box] += 1.0
                    points.append([x,y])
                    count_particle = count_particle + 1
                    index_particle.append(count_particle)

            # Calculus of textures, B and T##################
            number_particles = len(points)
            if number_particles > 0:
                points         = np.array(points)
                list_neighbors = delaunay(points,max_dist)
                map_focus_region_to_part(points, list_neighbors, index_particle)
                map(lambda i:i.texture(), part)
                for i in index_particle:
                    xx  = int((part[i].r[0]-x0) / box_size)
                    yy  = int((part[i].r[1]-y0) / box_size) * box_per_line_x
                    box = xx + yy
                    texture_box[box]+=part[i].M
                    if  count_events > 1 :
                        B_box[box] += part[i].B
                        T_box[box] += part[i].T
                       
                map(lambda i:i.copy_to_old(), part)

            #Calculate the average velocity over boxes
            for box in range(box_total):
                if density_now[box] > 0 :
                    vx_now[box] = vx_now[box] / density_now[box]
	            vy_now[box] = vy_now[box] / density_now[box]
                    texture_box[box]            = texture_box[box] / density_now[box]
                    ax_a,ax_b,angle             = axis_angle(texture_box[box])
                    axis_a_texture[box]     = ax_a
                    axis_b_texture[box]     = ax_b
                    ang_elipse_texture[box] = angle
                    if count_events > 1 :
                        B_box[box]            = B_box[box] / density_now[box]
                        ax_a,ax_b,angle       = axis_angle(B_box[box])
                        axis_a_B[box]     = ax_a
                        axis_b_B[box]     = ax_b
                        ang_elipse_B[box] = angle
                        T_box[box]            = T_box[box] / density_now[box]
                        ax_a,ax_b,angle       = axis_angle(T_box[box])
                        axis_a_T[box]     = ax_a
                        axis_b_T[box]     = ax_b
                        ang_elipse_T[box] = angle

            #Function call to write velocity-density gnu script
            velocity_density_script(box_per_line_x, box_per_column_y, x, y, vx_now, vy_now, density_now, system_type, image, v0)


            #Function call to write texture elipsis gnu script
            texture_elipsis_script_simu(box_per_line_x, box_total, axis_a_texture, axis_b_texture,\
                ang_elipse_texture, image-image_0,points,x0,y0,box_size)
            if count_events > 1 :
                B_elipsis_script_simu(box_per_line_x, box_total, axis_a_B, axis_b_B,\
                                      ang_elipse_B, image-image_0,points,x0,y0,box_size)
                T_elipsis_script_simu(box_per_line_x, box_total, axis_a_T, axis_b_T,\
                                      ang_elipse_T, image-image_0,points,x0,y0,box_size)
                        

        #Summing each box at different times
            for box in range(box_total) :
                density_tot[box] += density_now[box]
                vx_tot[box]      += vx_now[box]
                vy_tot[box]      += vy_now[box]
                texture_tot[box] += texture_box[box]
                B_tot[box] += B_box[box]
                T_tot[box] += T_box[box]
            #reseting matrices of instaneous measures
            vx_now      = list(0. for i in range(box_total))
            vy_now      = list(0. for i in range(box_total))
            density_now = list(0  for i in range(box_total))
            texture_box = list(np.zeros((2,2)) for i in range(box_total))
            B_box = list(np.zeros((2,2)) for i in range(box_total))
            T_box = list(np.zeros((2,2)) for i in range(box_total))
            points            = []
            index_particle    = []
        
    os.system("rm files.dat");                        


# Before starting time averages we exclude box at the borders. 
r_obst = R_OBST / box_size
x_obst = (X_OBST-x0) / box_size
y_obst = (Y_OBST-y0) / box_size


if system_type == 'experiment':

    zero_borders_and_obstacle_experiment(box_per_line_x, box_per_column_y, r_obst, x_obst, y_obst, density_tot, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type)

    
    # Here we write the time averages of density, velocity and deformation elipse for experiment
    average_density_velocity_deformation_experiment(box_per_line_x, box_per_column_y, x, y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, image_counter)

    
    # Five axis analysis for experiment
    five_axis_experiment(box_total, box_per_line_x, box_per_column_y, vx_tot, vy_tot, axis_a_tot, axis_b_tot, ang_elipse_tot, system_type, image_counter)
    
else:

#    zero_borders_and_obstacle_simu(box_per_line_x, box_per_column_y, r_obst, x_obst, y_obst, density_tot, vx_tot, vy_tot, texture_tot, B_tot, T_tot, system_type)
    
    # Here we write the time averages of density, velocity and deformation elipse for simus
    vx_win, vy_win, vx2_win, vy2_win,  density_win, texture_win, B_win, T_win = average_density_velocity_deformation(box_per_line_x, box_per_column_y, vx_tot, vy_tot,  \
    density_tot, texture_tot, B_tot, T_tot, vx_win, vy_win, vx2_win, vy2_win,  density_win, texture_win, B_win, \
    T_win, count_events, v0, vel_win_file_name, vel_fluct_win_file_name, dens_win_file_name, path, image_counter, window_size)

    # Five axis analysis for simulations
    five_axis_simu(box_total, box_per_line_x, box_per_column_y, vx_win, vy_win, texture_win, B_win, T_win, system_type, image_counter,path)


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



        
