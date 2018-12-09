import math as math
import os
import sys
# To this program work correctly:
# Create the file named 'parameter.in', it will have one line with three words: <system_type>, <data_file>, <window_size>
# <system_type> is: 'experiment' or 'superboids' or 'voronoi' or 'boids' or 'szabo-boids' or ...
# <data_file> is: data file name. It must be in a directory named <system_type>
# <window_size> : this is the space averaging window size
# The outuput files will be in directory named 'output'

def read_param():
    while 1 :
        line = f.readline()
        print(line)
        if not line:
            break #EOF
        a=line.split()
        if a[0] == 'window' :
            window_size = int(a[1])
        if a[0] == 'time_0' :
            time_0 = int(a[1])
        if a[0] == 'time_f' :
            time_f = int(a[1])
        if a[0] == 'voronoi' :
            voronoi = a[1]
        if a[0] == 'x0' :
            x0 = float(a[1])
        if a[0] == 'xf' :
            xf = float(a[1])
        if a[0] == 'y0' :
            y0 = float(a[1])
        if a[0] == 'yf' :
            yf = float(a[1])
        if a[0] == 'file' :
            filename = a[1]
            print(filename)
            exit()
    return window_size,time_0,time_f,obstacle,voronoi,x0,xf,yf,filename



def box_variables_definition_simu(caixas_por_coluna,caixas_por_linha):
    caixas_total = caixas_por_coluna*caixas_por_linha
    vx_now = list(0. for i in range(caixas_total))
    vy_now = list(0. for i in range(caixas_total))
    density_now = list(0 for i in range(caixas_total))
    vx_tot = list(0. for i in range(caixas_total))
    vy_tot = list(0. for i in range(caixas_total))
    density_tot = list(0 for i in range(caixas_total))
    vx_win = list(0. for i in range(caixas_total))
    vy_win = list(0. for i in range(caixas_total))
    density_win = list(0 for i in range(caixas_total))
    eixo_a_tot = list(0. for i in range(caixas_total))
    eixo_b_tot = list(0. for i in range(caixas_total))
    ang_elipse_tot = list(0. for i in range(caixas_total))
    eixo_a_win = list(0. for i in range(caixas_total))
    eixo_b_win = list(0. for i in range(caixas_total))
    ang_elipse_win = list(0. for i in range(caixas_total))
    ratio=float(caixas_por_coluna)/caixas_por_linha
    vid_def.write("set size ratio %f  \n" % ratio)
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

    return caixas_total,ratio,vx_now,vy_now,density_now,vx_tot,vy_tot,density_tot,vx_win,vy_win,density_win,eixo_a_tot,eixo_b_tot,ang_elipse_tot,eixo_a_win,eixo_b_win,ang_elipse_win

def box_variables_definition_experiment(caixas_por_coluna,caixas_por_linha):
    caixas_total=caixas_por_coluna*caixas_por_linha
    ratio=float(caixas_por_coluna)/caixas_por_linha
    vx_tot = list(0. for i in range(caixas_total))
    vy_tot = list(0. for i in range(caixas_total))
    density_tot = list(0 for i in range(caixas_total))
    eixo_a_tot = list(0. for i in range(caixas_total))
    eixo_b_tot = list(0. for i in range(caixas_total))
    ang_elipse_tot = list(0. for i in range(caixas_total))
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

    return caixas_total,ratio,vx_tot,vy_tot,density_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot


def velocity_density_script(caixas_por_linha,caixas_por_coluna,x,y,vx_now,vy_now,density_now,system_type,image,v0):
    #Here we write each image to the gnuplot velocity-density movie script
    vid_veloc_dens.write("plot [%f:%f] [%f:%f] \'-\' u ($1):($2):(arrow*$3):(arrow*$4):($5) with vectors head size  0.6,20,60  filled palette title \"%d\"\n" % (0,caixas_por_linha,0,caixas_por_coluna,image))
    if system_type == "experiment":
        v0=1 #this should be the real velocity, if we can measure...
        density=1 #this should be changed if we measure real density (density_now)
        for i in range(caixas_total):
            module=math.sqrt(vx_now[i]**2+vy_now[i]**2)
            if module > 0. :
                vid_veloc_dens.write('%i %i %f %f %f %f\n' % (x[i],y[i],vx_now[i]/module,vy_now[i]/module,module,density)) #density_now should be used case we have it
    else:
        for box in range(caixas_total):
            if density_now[box] != 0 :
                x,y=box%caixas_por_linha,box/caixas_por_linha
                module=math.sqrt(vx_now[box]**2+vy_now[box]**2)
                if module > 0.:
                    vid_veloc_dens.write('%i %i  %f %f %f %f %d\n' % (x,y,vx_now[box]/module,vy_now[box]/module,module/v0,density_now[box],box))
    vid_veloc_dens.write("e \n")
    vid_veloc_dens.write("pause .1 \n")

def deformation_elipsis_script(x,y,eixo_b,eixo_a,ang_elipse,system_type):
    #Deformation elipsis gnuplot script
    if system_type == 'experiment':
        for i in range(caixas_total):
            vid_def.write("set object %i ellipse at %i,%i size %f,0.5 angle %f \n" % (i+1,x[i],y[i],eixo_b[i]/(2*eixo_a[i]),ang_elipse[i]))
        vid_def.write("plot \'-\' w d notitle\n")
        for i in range(caixas_total):
            vid_def.write("%i %i \n" % (x[i],y[i]))
        vid_def.write("e \n")
        vid_def.write("pause .1 \n")
        vid_def.write("unset for [i=1:%i] object i \n" % (caixas_total+1))

def  zero_borders_and_obstacle(caixas_por_linha,caixas_por_coluna,r_obst,x_obst,y_obst,density_tot,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot,system_type):

    center_x = caixas_por_linha/2
    center_y = caixas_por_coluna/2
    for i in range(caixas_total):
        if system_type == 'experiment':
            bx = int(i/caixas_por_coluna)
            by = i%caixas_por_coluna
        else:
            bx = i%caixas_por_linha
            by = int(i/caixas_por_linha)        
        if bx == 0 or bx == caixas_por_linha-1 or by == 0 or by == caixas_por_coluna-1 :
            density_tot[i] = -10
            vx_tot[i] = 0.
            vy_tot[i] = 0.
            eixo_a_tot[i]=0.
            eixo_b_tot[i]=0.
            ang_elipse_tot[i]=0.
        if system_type == 'experiment' :
            if math.sqrt((bx-x_obst)**2+(by-y_obst)**2) < r_obst :
                density_tot[i] = -10
                vx_tot[i] = 0.
                vy_tot[i] = 0.
                eixo_a_tot[i] = 0.
                eixo_b_tot[i] = 0.
                ang_elipse_tot[i] = 0.
        if system_type == 'superboids' :
            if math.sqrt((bx-center_x)**2+(by-center_y)**2) < r_obst :
                density_tot[i] = -10
                vx_tot[i] = 0.
                vy_tot[i] = 0.
                eixo_a_tot[i]=0.
                eixo_b_tot[i]=0.
                ang_elipse_tot[i]=0.

    return caixas_por_linha,caixas_por_coluna,density_tot,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot

def average_density_velocity_deformation_experiment(caixas_por_linha,caixas_total,x,y,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,image_counter):
    arrow=1.5

    for i in range(caixas_total):
        vx_tot[i]/=image_counter
        vy_tot[i]/=image_counter

    caixas_por_coluna=caixas_total/caixas_por_linha
    vel_win.write("plot [%f:%f] [%f:%f] \'-\' u ($1):($2):(%f*$3):(%f*$4):($5)  with vectors notitle head size  0.3,20,60  filled palette \n" % (0,caixas_por_linha,0,caixas_por_coluna,arrow,arrow))

    for i in range(caixas_total):
        dens_win.write("%d %d %f \n" % (x[i],y[i],density_tot[i]))
        module=math.sqrt(vx_tot[i]**2+vy_tot[i]**2)
        if module >0 :
            vel_win.write("%d %d %f %f %f\n" % (x[i],y[i],vx_tot[i]/module,vy_tot[i]/module,module))
            def_win.write("%d %d %f %f %f \n" % (x[i],y[i], eixo_a_tot[i],eixo_b_tot[i],ang_elipse_tot[i]))


    vel_win.write("e \n")
    vel_win.write("pause -1 \n")

def average_density_velocity_deformation(caixas_por_linha,caixas_total,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,density_tot,vx_win,vy_win,eixo_a_win,eixo_b_win,density_win,count_events,v0):
    caixas_por_coluna=caixas_total/caixas_por_linha
    if system_type != 'experiment':
        count_box_win = list(0 for i in range(caixas_total))
        for bx in range(window_size+1,caixas_por_linha-window_size-1):
            for by in range(window_size+1,caixas_por_coluna-window_size-1):
                for k in range(-window_size,window_size):
                    for l in range(-window_size,window_size):
                        if density_tot[(bx+k)+(by+l)*caixas_por_linha]>0:
                            box = bx + (by*caixas_por_linha)
		            density_win[box] += density_tot[(bx+k)+((by+l)*caixas_por_linha)]
		            vx_win[box]      += vx_tot[(bx+k)+((by+l)*caixas_por_linha)]
		            vy_win[box]      += vy_tot[(bx+k)+((by+l)*caixas_por_linha)]
		            count_box_win[box] += 1
                        else:
                            box = bx + (by*caixas_por_linha)
                            count_box_win[box] += 1
                            
        #Average win calculus and data print (gnuplot script for velocity)

        caixas_por_coluna=caixas_total/caixas_por_linha
        
        vel_win.write("plot [%f:%f] [%f:%f] \'-\' u ($1):($2):(arrow*$3):(arrow*$4):($5)  with vectors notitle  head size  0.3,20,60  filled palette \n" % (0,caixas_por_linha,0,caixas_por_coluna))



        for bx in range(window_size+1,caixas_por_linha-window_size-1):
            for by in range(window_size+1,caixas_por_coluna-window_size-1):
            	box = bx + (by*caixas_por_linha)
                module = math.sqrt((vx_win[box]*vx_win[box])+(vy_win[box]*vy_win[box]))
                if bx%4 == 0 and by%4 == 0 :
	            if density_win[box]>0.0 and module > 0.0 :
                        normalization=float(count_events*count_box_win[box])*v0
                        
	                dens_win.write("%d %d %f \n" % (bx, by,density_win[box]/normalization))
	                vel_win.write("%d %d %f %f %f %f %f \n"% (bx, by, vx_win[box]/module, vy_win[box]/module, module/normalization, density_tot[box]/float(count_events), density_win[box]/normalization))
	            else :
	                vx_win[box] = 0.0
	                vy_win[box] = 0.0
	                dens_win.write("%d %d %f \n"%(bx, by, 0.0))
                        vel_win.write("%d %d %f %f %f %f %f \n" % (bx, by, 0.0, 0.0, 0.0,0.0, 0.0))
        vel_win.write("e \n")
        vel_win.write("pause -1 \n")
                 
    return vx_win,vy_win,eixo_a_win,eixo_b_win,density_win




def five_axis(caixas_total,caixas_por_linha,caixas_por_coluna,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot,system_type,image_counter):
    caixas_meia_altura=caixas_por_coluna/2
    caixas_quarto_altura=caixas_por_coluna/4
    caixas_meia_largura = caixas_por_linha/2
    caixas_quarto_largura = caixas_por_linha/4
    vx_axis1,vx_axis2,vx_axis3,vx_axis4,vx_axis5,vx_axis6=[],[],[],[],[],[]
    vy_axis1,vy_axis2,vy_axis3,vy_axis4,vy_axis5,vy_axis6=[],[],[],[],[],[]
    eixo_a_axis1,eixo_a_axis2,eixo_a_axis3,eixo_a_axis4,eixo_a_axis5,eixo_a_axis6=[],[],[],[],[],[]
    eixo_b_axis1,eixo_b_axis2,eixo_b_axis3,eixo_b_axis4,eixo_b_axis5,eixo_b_axis6=[],[],[],[],[],[]
    ang_elipse_axis1,ang_elipse_axis2,ang_elipse_axis3,ang_elipse_axis4,ang_elipse_axis5,ang_elipse_axis6=[],[],[],[],[],[]

    
    for i in range(caixas_total):
        if system_type == 'experiment':
            bx = int(i/caixas_por_coluna)
            by = i%caixas_por_coluna
        else:
            bx=i%caixas_por_linha
            by=int(i/caixas_por_linha)
        if by == caixas_meia_altura :
            vx_axis1.append(vx_tot[i])
            vy_axis1.append(vy_tot[i])
            eixo_a_axis1.append(eixo_a_tot[i])
            eixo_b_axis1.append(eixo_b_tot[i])
            ang_elipse_axis1.append(ang_elipse_tot[i])

        if by == caixas_quarto_altura:
            vx_axis2.append(vx_tot[i]/2.)
            vy_axis2.append(vy_tot[i]/2.)
            eixo_a_axis2.append(eixo_a_tot[i]/2.)
            eixo_b_axis2.append(eixo_b_tot[i]/2.)
            ang_elipse_axis2.append(ang_elipse_tot[i]/2.)

        if by == 3*caixas_quarto_altura:
            vx_axis6.append(vx_tot[i]/2.)
            vy_axis6.append(vy_tot[i]/2.)
            eixo_a_axis6.append(eixo_a_tot[i]/2.)
            eixo_b_axis6.append(eixo_b_tot[i]/2.)
            ang_elipse_axis6.append(ang_elipse_tot[i]/2.)

        if bx == caixas_meia_largura:
            vx_axis3.append(vx_tot[i]/2.)
            vy_axis3.append(vy_tot[i]/2.)
            eixo_a_axis3.append(eixo_a_tot[i]/2.)
            eixo_b_axis3.append(eixo_b_tot[i]/2.)
            ang_elipse_axis3.append(ang_elipse_tot[i]/2.)

        if bx == caixas_meia_largura-caixas_quarto_largura:
            vx_axis4.append(vx_tot[i]/2.)
            vy_axis4.append(vy_tot[i]/2.)
            eixo_a_axis4.append(eixo_a_tot[i]/2.)
            eixo_b_axis4.append(eixo_b_tot[i]/2.)
            ang_elipse_axis4.append(ang_elipse_tot[i]/2.)

        if bx == caixas_meia_largura+caixas_quarto_largura:
            vx_axis5.append(vx_tot[i]/2.)
            vy_axis5.append(vy_tot[i]/2.)
            eixo_a_axis5.append(eixo_a_tot[i]/2.)
            eixo_b_axis5.append(eixo_b_tot[i]/2.)
            ang_elipse_axis5.append(ang_elipse_tot[i]/2.)


    for i in range(len(vx_axis2)):
            vx_axis2[i] += vx_axis6[i]
            vy_axis2[i] += vy_axis6[i]
            eixo_a_axis2[i] += eixo_a_axis6[i]
            eixo_b_axis2[i] += eixo_b_axis6[i]
            ang_elipse_axis2[i] += ang_elipse_axis6[i]

    for i in range(caixas_por_linha):
        mm.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis1[i], vy_axis1[i], eixo_a_axis1[i], eixo_b_axis1[i], ang_elipse_axis1[i]))
        nn.write("%d %f %f %f %f %f\n" % (i-caixas_meia_largura, vx_axis2[i], vy_axis2[i], eixo_a_axis2[i], eixo_b_axis2[i], ang_elipse_axis2[i]))


    for i in range(caixas_por_coluna):
        oo.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis3[i], vy_axis3[i], eixo_a_axis3[i], eixo_b_axis3[i], ang_elipse_axis3[i]))
        pp.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis4[i], vy_axis4[i], eixo_a_axis4[i], eixo_b_axis4[i], ang_elipse_axis4[i]))
        qq.write("%d %f %f %f %f %f\n" % (i-caixas_meia_altura, vx_axis5[i], vy_axis5[i], eixo_a_axis5[i], eixo_b_axis5[i], ang_elipse_axis5[i]))


def imag_count(system_type):
    counter=0
    print "Counting images... wait... it may take 5s to count 1000 images\n"
    if system_type == 'superboids' :
        while 1 :
            line = fn.readline()
            if not line :
                break #EOF
            if line.replace( '\r', '' ) == '\n' : #identifies blank lines
                counter+=1
                while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                    line = fn.readline()
    if system_type == 'experiment' :
        while 1 :
            line = f.readline()
            if not line:
                break #EOF
            a=line.split()
            if a[0] == 'Time_start:' :
                counter += 1
                
    if system_type == 'szabo-boids' :
        while 1:
            line = fd.readline()
            if not line:
                break
            a=line.split()
            if a[0] == 'x' :
                counter += 1

    if system_type == 'boids-particle' :
        while 1:
            line = fd.readline()           
            if not line:
                break #EOF
            a=line.split()
            counter += 1
            print a
            for i in range(int(a[1])):
                fd.readline()
            
    print "Counted", counter, "images.\n"
    print "Type initial and final image number you want to analyse (min=0, max=",counter,") - Use spaces to separate the two numbers"

################## Here starts the main program ###############

#Opening input parameter file

f=open("parameter.in")
a=f.readline().split()

system_type=a[1]
fileout=a[1]
path='output/'+system_type
print path

#Creating the directory structure for output
os.system('mkdir -p %s' % path)
os.system('cp axis12.par %s' % path)
os.system('cp axis345.par %s' % path)
#velocity_density gnuplot file

vid_veloc_dens=open("%s/video_velocity_density.gnu"%path,"w")

#deformation elipse script header

vid_def=open("%s/video_deformation.gnu"%path,"w")
vid_def.write("unset key \n")

#Opening time averages files

dens_win=open("%s/density-win.dat"%path,"w")
vel_win=open("%s/velocity-win.dat"%path,"w")

def_win=open("%s/deformation-win.dat"%path,"w")

# Opening five axis analysis files

mm = open("%s/axis1.dat"%path,"w")
nn = open("%s/axis2.dat"%path,"w")
oo = open("%s/axis3.dat"%path,"w")
pp = open("%s/axis4.dat"%path,"w")
qq = open("%s/axis5.dat"%path,"w")


if system_type == 'experiment':
    arq_in="%s/%s"%(a[0],a[1])
    print "You analise an", system_type, "system, reading data from file:\n", arq_in
    window_size=int(a[2])
    #Opening the data file
    f=open(arq_in)
    imag_count(system_type)
    f.close()
    a=sys.stdin.readline().split()
    image_0=int(a[0])
    image_f=int(a[1])
    f=open(arq_in)
    
    a=['0']
    
    # Reading file head (I have taken some lines of the header, you may want others)
    while(a[0]!='X'): #'X' marks the line just before data in experiment data file, that is, the end of the header
        a=f.readline().split()
        if(a[0]=='Box_end:'):
            caixas_por_coluna,caixas_por_linha=int(a[1]),int(a[2])
        if(a[0]=='Box_size:'):
            box_size=int(a[1])/4
        if(a[0]=='Obstacle_diameter:'):
            R_OBST=float(a[1])
            R_OBST=R_OBST/2.
        if(a[0]=='Maximum_observed_velocity:'):
            speed=float(a[1])
        if a[0]=='Obstacle_position:':
            X_OBST=int(a[1])
            Y_OBST=int(a[2])


    caixas_total,ratio,vx_tot,vy_tot,density_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot=box_variables_definition_experiment(caixas_por_coluna,caixas_por_linha)


    #Reading x,y,density,vx,vy data on experiment file

    image=0
    v0=1
    while 1 :
        line = f.readline()
        if not line : break # EOF
        a=line.split()
#        print a
        counter=0
        image_counter=image_f-image_0
        x,y,eixo_a,eixo_b,ang_elipse,vx_now,vy_now,density_now=[],[],[],[],[],[],[],[]
        while a[0]!='Time_start:' :
            if image >= image_0 :
                if a[7] == 'NaN': 
                    vx_now.append(0.)
                    vy_now.append(0.)
                if a[7] != 'NaN' : 
                    vx_now.append(float(a[7]))
                    vy_now.append(float(a[8]))
                x.append(int(a[0])/box_size)
                y.append(int(a[1])/box_size)
                eixo_a.append(float(a[4]))
                eixo_b.append(float(a[5]))
                ang_elipse.append(float(a[6])/(math.pi)*180)
                density_now.append(float(a[9]))
                vx_tot[counter]+=vx_now[counter]
                vy_tot[counter]+=vy_now[counter]
                density_tot[counter]+=density_now[counter]
                eixo_a_tot[counter]+=eixo_a[counter]
                ang_elipse_tot[counter]+=ang_elipse[counter]
                counter += 1
            line=f.readline()
            if not line : break # EOF
            a=line.split()
        image += 1
        line = f.readline()
        line = f.readline()

        if image < image_0 : print "Skipping image ",image
        if image > image_0 and image <= image_f :
            print "Analising image ",image, "..."
            #Function call to write velocity-density gnu script
            velocity_density_script(caixas_por_linha,caixas_por_coluna,x,y,vx_now,vy_now,density_now,system_type,image,v0)
            #Function call to write deformation elipse gnu script
            deformation_elipsis_script(x,y,eixo_b,eixo_a,ang_elipse,system_type)
        elif image > image_f:
            break
        

if system_type == "superboids":
    arq_header_in = "%s/%s.dat"%(a[0],a[1])
    arq_data_in = "%s/%s_plainprint.dat"%(a[0],a[1])
    arq_neigh_in = "%s/%s_neighbors.dat"%(a[0],a[1])
    print "\nYou analise a", a[0], "system, data is read from files:\n", arq_header_in," (header)\n", arq_data_in," (data)\n", arq_neigh_in," (neighbors)"
    fh=open(arq_header_in)
    fd=open(arq_data_in)
    fn=open(arq_neigh_in)
    window_size=int(a[2])
    imag_count(system_type)
    fn.close()
    fn=open(arq_neigh_in)
    a=sys.stdin.readline().split()
    image_0=int(a[0])
    image_f=int(a[1])
    v0=0.007
    # Reading superboids parameter file 

    while 1 :
        line = fh.readline()
        if not line : break # EOF
        if line.replace( '\r' , '' ) == "\n" :
            continue
        else :
            a=line.split()
            #            print a
            if a[0] == '#' and a[1]=='RECTANGLE:':
                a=fh.readline().split()
                Lx,Ly = int(float(a[1])),int(float(a[2]))
                Lx = 240 # corrigindo por enquanto o erro no arquivo de entrada. REVISAR!!!
            if a[1] == 'Radial' and a[2]=='R_Eq:':
                a=fh.readline().split()
                box_size  = int(float(a[1]))
                #                print box_size
            if a[0] == "#radius:":
                R_OBST = int(a[1])
                X_OBST = float(a[3])
                Y_OBST = float(a[4])
                
    caixas_por_linha,caixas_por_coluna = Lx/box_size,Ly/box_size

    #defining all matrices
    caixas_total,ratio,vx_now,vy_now,density_now,vx_tot,vy_tot,density_tot,vx_win,vy_win,density_win,eixo_a_tot,eixo_b_tot,ang_elipse_tot,eixo_a_win,eixo_b_win,ang_elipse_win=box_variables_definition_simu(caixas_por_coluna,caixas_por_linha)
    x0,y0 = -Lx/2,-Ly/2
    xf,yf = Lx/2,Ly/2

    
    #Reading superboids plainprint data file
    fat_boids_counter = 0
    image = 0
#    image_0=1
#    image_f=100000
    count_events=0
    #    mx=[0,0,0]
    while 1 :
        line = fd.readline()
        if not line :
            vid_veloc_dens.write("pause -1 \n")
            break #EOF
        if line.replace( '\r', '' ) == '\n' : #identifies blank lines
            if image > image_0 and image <= image_f:
                #Calculate the average velocity over boxes
                for box in range(caixas_total):
                    if density_now[box] > 0 :
                        vx_now[box] = vx_now[box] / density_now[box]
		        vy_now[box] = vy_now[box] / density_now[box]
                #Function call to write velocity-density gnu script
                velocity_density_script(caixas_por_linha,caixas_por_coluna,x,y,vx_now,vy_now,density_now,system_type,image,v0)
            #Summing each box at different times
            if image > image_0 and image <= image_f:
                for box in range(caixas_total) :
                    density_tot[box] += density_now[box]
                    vx_tot[box] += vx_now[box]
                    vy_tot[box] += vy_now[box]

#            print image_counter
            image+=1
            if image <= image_0 :
                print "Skippin image:",image 
            elif image <= image_f :
                print "Image number",image,". Number of super particles =", fat_boids_counter
                count_events+=1

            else:
                break
            fat_boids_counter = 0
            while line.replace( '\r' , '' ) == '\n' : # skip blank lines
                line = fd.readline()
        else:
            a = line.split()
            if image >= image_0 and image <= image_f:
                if float(a[6]) > 0.5 :
                    fat_boids_counter+=1
                    x,y=float(a[0]),float(a[1])
                    if x>x0 and x<xf and y>y0 and y<yf:
                        xx=int((x-x0)/box_size)
                        yy=int((y-y0)/box_size) * caixas_por_linha
                        box = xx+yy
                        #                    if mx[0]<box : mx=[box,xx,yy]
                        vx_now[box]+=float(a[2])
                        vy_now[box]+=float(a[3])
                        density_now[box] += 1.0
    image_counter=image_f-image_0


if system_type == "szabo-boids":
    arq_data_in = "%s/%s.dat"%(a[0],a[1])
    print "\nYou analise a", a[0], "system, data is read from files:\n", arq_data_in
    fd=open(arq_data_in)
#    fn=open(arq_neigh_in)
    window_size=int(a[2])
    imag_count(system_type)
    fd.close()
    fd=open(arq_data_in)
    a=sys.stdin.readline().split()
    image_0=int(a[0])
    image_f=int(a[1])
    image_counter=image_f-image_0
    image =0
    line_counter=0
    count_events=0
    v0=0.1
    #Reading szabo-boids  data file
    while 1 :
        line = fd.readline()
        if not line : break
        line_counter +=1
        if line.replace( '\r' , '' ) == "\n" :
            continue
        a=line.split()
        if line_counter < 7 :
            if a[0] == "Number_of_particles:" : N = int(a[1])
            if a[0] == "Box_size:" : box_size = int(a[1])
            if a[0] == "Steps_between_images:" : delta_images = int(a[1])
            if a[0] == "Dimensions:" :

                Lx = int(a[1])
                Ly = int(a[2])
                caixas_por_linha,caixas_por_coluna = Lx/box_size,Ly/box_size
                x0=-Lx/2
                y0=-Ly/2
                xf=Lx/2
                yf=Ly/2
                caixas_total,ratio,vx_now,vy_now,density_now,vx_tot,vy_tot,density_tot,vx_win,vy_win,density_win,eixo_a_tot,eixo_b_tot,ang_elipse_tot,eixo_a_win,eixo_b_win,ang_elipse_win=box_variables_definition_simu(caixas_por_coluna,caixas_por_linha)


            if a[0] == "Radius:" : R_OBST = int(a[1])
            if a[0] == "Obst_position:" :
                X_OBST = int(a[1])
                Y_OBST = int(a[2])
        if line_counter > 6 :
            if image <= image_0 :
                print "Skipping image:",image 
                while a[0] != 'x' :
                    line = fd.readline()
                    if not line : break
                    a =line.split()
                image+=1
            elif image <= image_f :
                boids_counter = 0
                while a[0] != 'x' :
                    boids_counter += 1
                    x,y=float(a[0]),float(a[1])
                    if x>x0 and x<xf and y>y0 and y<yf:
                        xx=int((x-x0)/box_size)
                        yy=int((y-y0)/box_size) * caixas_por_linha
                        box = xx+yy
                        vx_now[box]+=float(a[2])
                        vy_now[box]+=float(a[3])
                        density_now[box] += 1.0
                        line = fd.readline()
                        if not line : break
                        a =line.split()
                image+=1
                count_events+=1
                print image, boids_counter
                #Calculate the average velocity over boxes
                for box in range(caixas_total):
                    if density_now[box] > 0 :
                        vx_now[box] = vx_now[box] / density_now[box]
		        vy_now[box] = vy_now[box] / density_now[box]
                        #Function call to write velocity-density gnu script
                velocity_density_script(caixas_por_linha,caixas_por_coluna,x,y,vx_now,vy_now,density_now,system_type,image,v0)
                #Summing each box at different times
                for box in range(caixas_total) :
                    density_tot[box] += density_now[box]
                    vx_tot[box] += vx_now[box]
                    vy_tot[box] += vy_now[box]

            else:
                break

if system_type == "vicsek":

    window_size,time_0,time_f,obstacle,voronoi,x0,xf,yf,filename=read_param()
    arq_data_in = "%s/%s.dat"%(a[0],a[1])
    print "\nYou analise a", a[0], "system, data is read from files:\n", arq_data_in
    fd=open(arq_data_in)
#    fn=open(arq_neigh_in)
    imag_count(system_type)
    fd.close()
    fd=open(arq_data_in)
    a=sys.stdin.readline().split()
    image_0=int(a[0])
    image_f=int(a[1])
    image_counter=image_f-image_0
    image =0
    line_counter=0
    count_events=0
    v0=0.1    

# Before starting time averages we exclude box at the borders. 
r_obst=R_OBST/box_size
x_obst=X_OBST/box_size
y_obst=Y_OBST/box_size
zero_borders_and_obstacle(caixas_por_linha,caixas_por_coluna,r_obst,x_obst,y_obst,density_tot,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot,system_type)


if system_type == 'experiment':

    # Here we write the time averages of density, velocity and deformation elipse for experiment
    average_density_velocity_deformation_experiment(caixas_por_linha,caixas_total,x,y,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,image_counter)

#    for i in range(caixas_total):
#        print i%caixas_por_linha,i/caixas_por_linha,x[i],y[i],vx_tot[i]
    
    # Five axis analysis for experiment
    five_axis(caixas_total,caixas_por_linha,caixas_por_coluna,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,ang_elipse_tot,system_type,image_counter)

    
else:
    # Here we write the time averages of density, velocity and deformation elipse for simus
    vx_win,vy_win,eixo_a_win,eixo_b_win,density_win=average_density_velocity_deformation(caixas_por_linha,caixas_total,vx_tot,vy_tot,eixo_a_tot,eixo_b_tot,density_tot,vx_win,vy_win,eixo_a_win,eixo_b_win,density_win,count_events,v0)

    # Five axis analysis for simulations
    five_axis(caixas_total,caixas_por_linha,caixas_por_coluna,vx_win,vy_win,eixo_a_win,eixo_b_win,ang_elipse_win,system_type,image_counter)





f.close()
vid_veloc_dens.close()
vid_def.close()
dens_win.close()
vel_win.close()
def_win.close()
mm.close()
nn.close()
oo.close()
pp.close()
qq.close()



        
