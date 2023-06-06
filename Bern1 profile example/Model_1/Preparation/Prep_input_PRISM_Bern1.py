# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 18:32:22 2023

@author: Dimitri
"""
import numpy as np
import datetime
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as patches
import itertools
from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.path import Path
from matplotlib import gridspec

############ INPUT ##############

#Path to the folder that will contain the model data. If this folder does not exist, it will be created
run_folder_path="../Data_files"

#Name of the model this input is written for
model_name="Bern1_model1"

initial_input_file='initial_input_bern1.txt'

#Indicate which prisms to plot on the LCS 2D map view figure
#Position in the file expected example: [0,2] will show the first prism and the last prism
show_prism=[0]

#Position of the maximum anomaly station, if this option is left on auto, the maximum anomaly station will be automatically plotted
#Use this option to plot more than one station
max_ano_mp="auto"

#Flag to check if a station or drilling is far from the gravity profile
#Provides the distance perpendicular to the profile. The flag_distance variable will mark points further than the set threshold
check_dist="yes"
flag_distance=250

xlimits_2D=[-2.5,2.5]
ylimits_2D=[-2.5,2.5]


########## DATA IMPORT ####################

#Fixed length of the header before the MPs input
length_header=34

with open(initial_input_file) as file_input:
   nb_file_input=file_input.readlines()
   
#Angles
#For the projection onto the gravity profile
   alpha=float(nb_file_input[6])
#For LCS rotation and back rotation to CH coor
   theta=float(nb_file_input[9])
#Gravity profile length - This value is extracted from the DEM profile file
   gp_length=float(nb_file_input[12])
#Load the number of stations from the file, this position is the line number in the file and is fixed
#Used to extract the input data and to check the user's input
   nb_stations=int(nb_file_input[24])
   
#Load the number of drillings from the file, this position is the line number in the file and is fixed
#Only needed to extract the input data. Drillings values are not used in this code
   nb_drillings=int(nb_file_input[27])
   
#Load the number of prisms from the file, this position is the line number in the file and is fixed
#Used to extract the input data and to check the user's input
   nb_prisms=int(nb_file_input[30])


#### CH COORDINATES IMPORT

#IF UNPACK=TRUE ALL X/Y/Z ARE STORED IN THE SAME ARRAY E.G. ARRAY1= [X1,X2,X3,....], ARRAY2= [Y1,Y2,Y3,....], ETC...
#IF UNPACK=FALSE EACH ARRAY STORES COUPLES OF COORDINATES (X,Y,Z)

#Gravity profile and LCS origin data import
#gp_LCS_data=np.loadtxt(str(gp_file_path), unpack=True)
 
#Gravity profile and LCS origin coordinates import
#Skipheader is the fixed number of lines at the start of the file
gp_LCS_data_header=np.array([15,18,21])
#Skipfooter is the number of drillings, the number of prisms, the number of intersections prism/gravity profiles and the max depth coordinates
gp_start_data_footer=5 + nb_stations + nb_drillings + 4 + nb_prisms
gp_end_data_footer=4 + nb_stations + nb_drillings + 4 + nb_prisms
LCS_data_footer=3 + nb_stations + nb_drillings + 4 + nb_prisms
#Gravity profile start coordinates
gp_start=np.genfromtxt(initial_input_file, skip_header=gp_LCS_data_header[0], skip_footer=gp_start_data_footer)
#Gravity profile end coordinates
gp_end=np.genfromtxt(initial_input_file, skip_header=gp_LCS_data_header[1], skip_footer=gp_end_data_footer)
#LCS origin
LCS_ori=np.genfromtxt(initial_input_file, skip_header=gp_LCS_data_header[2], skip_footer=LCS_data_footer)

#Getting the MPs data
#Skipheader is the fixed number of lines at the start of the file
MP_data_header=length_header
#Skipfooter is the number of drillings, the number of prisms, the number of intersections prism/gravity profiles and the max depth coordinates
MP_data_footer=nb_drillings + 4 + nb_prisms 
#Load the stations coordinates
MP_all_data =  np.genfromtxt(initial_input_file, skip_header=MP_data_header, skip_footer=MP_data_footer, unpack=True)

#Getting the drillings data
#Skipheader is the fixed number of lines at the start of the file
drillings_header= length_header + nb_stations + 2
#Skipfooter is the top prism corner coordinates and the number of prisms
drillings_footer=4 + nb_prisms 
#Load the drillings coordinates
drillings = np.genfromtxt(initial_input_file, skip_header=drillings_header, skip_footer=drillings_footer, unpack=True)

#Getting the prism data
#Skipheader is the fixed number of lines at the start of the file and the MPs,the drillings and their header and the prisms coordinates header
prism_ch_data_header= drillings_header + nb_drillings +  3
#Skipfooter is the number of prisms
prism_ch_data_footer=nb_prisms
#Load the top prism corner CH coordinates
pri_ch_data = np.genfromtxt(initial_input_file, skip_header=prism_ch_data_header, skip_footer=prism_ch_data_footer, unpack=True)

#Getting the prism data
#Skipheader is the fixed number of lines at the start of the file and the MPs,the drillings and their header and the prism corner coordinates and the prisms header 
prism_data_header= prism_ch_data_header + 7
#No data to skip at the footer
prism_data_footer=0
#Load the prisms LCS coordinates and density values
pri_lcs_data = np.genfromtxt(initial_input_file, skip_header=prism_data_header, skip_footer=prism_data_footer, unpack=True)


pri_ch_E, pri_ch_N, pri_ch_elev= pri_ch_data

prism_x2, prism_x1, prism_y2, prism_y1, prism_z2, prism_z1, rho_pri=pri_lcs_data

#### ADD to BA_plot ####
#Regional gravity stations
#regio_pts=np.loadtxt(str(first_dir)+"/"+str(regio_MP_file_path), unpack=True)

#### PROJECTION DATA (distance from profile start, CH coordinates and elevation data)

#DEM and bedrock profile import
#DEM_cross_data=distance, DEM_E, DEM_N, elev=np.loadtxt(str(first_dir)+"/"+str(sec_dir)+"/"+str(DEM_file_path), unpack=True)#Elevation profile data
#bedrock_cross_data=distance_b, bedrock_E, bedrock_N, elev_b=np.loadtxt(str(first_dir)+"/"+str(sec_dir)+"/"+str(bedrock_file_path), unpack=True)#Bedrock model profile data

#### ADD OPTION IN DISPLAY_PRISMA TO READ ANOTHER RESIDUAL ANOMALY DATA IF NEEDED
#### BREMGARTEN SPECIFIC - Some (5) stations (4031, 40310, 4039, 4040, 40400) are excluded from the PRISM modelling but included in the residual anomaly plot
#Since they stations are excluded from the PRISM modelling, and not in the PRISM MP list, they have to be added with another input used only for the residual anomaly plot.
#MP_residual_coor=np.loadtxt(str(first_dir)+"/MP_bremgarten.txt", unpack=True)#Position of the stations


############ FUNCTIONS ########################

#Convert the angle into radians
def rad_conv(angle):
    angle_rad=angle*np.pi/180
    return(angle_rad)

#Function to calculate the intermediate coordinates (LCS1) the new origin
def coor_shift(x,xstart):
    x_lcs1=x-xstart
    return(x_lcs1)

#Rotation equations, calculate x and y coordinates in a rotated coordinate system called local coordinate system (LCS)
#For the x coordinate calculations
def cos_x_rot(x,alpha):
    cos_x=x*np.cos(alpha)
    return(cos_x)

def sin_x_rot(y,alpha):
    sin_x=y*np.sin(alpha)
    return(sin_x)

def x_rotation(x,y,alpha):
    x_rotated=cos_x_rot(x,alpha)+sin_x_rot(y,alpha)
    return(x_rotated)

#For the y coordinate calculations
def cos_y_rot(y,alpha):
    cos_y=y*np.cos(alpha)
    return(cos_y)

def sin_y_rot(x,alpha):
    sin_y=-x*np.sin(alpha)
    return(sin_y)

def y_rotation(x,y,alpha):
    y_rotated=sin_y_rot(x,alpha)+cos_y_rot(y,alpha)
    return(y_rotated)

#Functions to flip the x and z coordinates
def invert_coor(x):
    x_invert=x*-1
    return(x_invert)

#Functions to calculate the a and b constants of y=ax+b, returns a and b 
#Takes two numpy arrays as input x=[x1,x2] and y=[y1,y2] with x1<x2 and y1<y2
def get_a_line_equation(x,y):
    y_diff = y[1] - y[0]
    x_diff = x[1] - x[0]
    if x_diff!=0:
        a=y_diff/x_diff
    return(a)

def get_b_line_equation(a,x,y):
    b=y-a*x
    return(b)

#function to calculate the x coordinate of the intersection between two lines
#Takes two numpy arrays as input a=[a1,a2] and b=[b1,b2], the constants from the linear equation of the two lines
def get_x_intersection(a,b):
    x_inter=-(b[1]-b[0])/(a[1]-a[0])
    return(x_inter)

#function to calculate the y coordinate of the intersection between two lines
#Takes three single variables as input a, b and x_inter, a and b being from one of the two intersecting lines
#No arrays
def get_y_intersection(a,b,x_inter):
    y_inter=a*x_inter+b
    return(y_inter)

def proj_gp(x,y,xgp,ygp,alpha):
    
    xlcs1=coor_shift(x,xgp)
    ylcs1=coor_shift(y,ygp)
    dist_ygp=y_rotation(xlcs1, ylcs1, alpha)
    if check_dist=="yes" or "YES" or "Yes":
        dist_xgp=x_rotation(xlcs1, ylcs1, alpha)
        return(dist_ygp,dist_xgp)
    else:
        return(dist_ygp)


def ch_to_lcs(x,y,z,xori,yori,zori,theta):

    xlcs1=coor_shift(x,xori)
    ylcs1=coor_shift(y,yori)
    zlcs1=coor_shift(z,zori)
    
    xlcs2=x_rotation(xlcs1,ylcs1,theta)
    ylcs2=y_rotation(xlcs1,ylcs1,theta)
    
    xlcs=invert_coor(xlcs2)
    ylcs=ylcs2
    zlcs=invert_coor(zlcs1)
    
    return(xlcs,ylcs,zlcs)

def lcs_to_ch(xlcs,ylcs,zlcs,xori,yori,zmin,theta):
    
    b_theta=invert_coor(theta)
    b_xori=invert_coor(xori)
    b_yori=invert_coor(yori)
    b_zmin=invert_coor(zmin)

    xlcs1=invert_coor(xlcs)
    ylcs1=ylcs
    zlcs1=invert_coor(zlcs)
            
    xlcs2=x_rotation(xlcs1,ylcs1,b_theta)
    ylcs2=y_rotation(xlcs1,ylcs1,b_theta)

    xb_ch=coor_shift(xlcs2,b_xori)
    yb_ch=coor_shift(ylcs2,b_yori)
    zb_ch=coor_shift(zlcs1,b_zmin)

    
    return(xb_ch,yb_ch,zb_ch)

def pri_corners(x1, x2, y1, y2):
    cornerbl=np.asarray([x2,y1])
    cornertl=np.asarray([x1,y1])
    cornertr=np.asarray([x1,y2])
    cornerbr=np.asarray([x2,y2])
    return(cornerbl,cornertl,cornertr,cornerbr)

def prism_verts(prism_x, prism_y):
    verts=[(prism_x[0], prism_y[0]),
           (prism_x[1], prism_y[1]),
           (prism_x[2], prism_y[2]),
           (prism_x[3], prism_y[3]),
           (prism_x[0], prism_y[0])]
    
    return verts

############ INPUT PROCESSING ########################
#Store the angles in degree
alpha_deg=alpha
theta_deg=theta

#convert angles in degree to radian
alpha=rad_conv(alpha)
theta=rad_conv(theta)

#Assign a variable to the MP numbers
MP_num=MP_all_data[0]

#Assign variables to the MPs x, y and z CH-coordinates
mp_xch=MP_all_data[1]
mp_ych=MP_all_data[2]
mp_zch=MP_all_data[3]

#Assign variables to the Bouguer anomaly, regional anomaly and residual anomaly
bouguer=MP_all_data[4]
regional=MP_all_data[5]
residual=MP_all_data[6]

#Assign variables to the drillings x, y and z CH-coordinates
drilling_number=drillings[0]
drillings_xch=drillings[1]
drillings_ych=drillings[2]
drillings_zch=drillings[3]

#ADD TO BA_PLOT
#Assign variables to the regional gravity stations x, y and z CH-coordinates
#regio_mp_xch=regio_pts[1]
#regio_mp_ych=regio_pts[2]


#Assign variables to the gravity profile and LCS origin CH-coordinates
gp_xstart=gp_start[0]
gp_ystart=gp_start[1]

gp_x_ch=np.array([gp_start[0],gp_end[0]])
gp_y_ch=np.array([gp_start[1],gp_end[1]])
gp_z_ch=np.array([gp_start[2],gp_end[2]])

LCS_xori=LCS_ori[0]
LCS_yori=LCS_ori[1]
zmin=LCS_ori[2]

user_angle_check=0

if user_angle_check == 1:
    user_alpha=float(input('Please confirm the value of the projection angle alpha: '))

    if user_alpha != alpha_deg:
        sys.exit("The projection angle you provided is not the same as given in the input, please check your data!")
    else:
        print("The projection angle values are consistent.")
    
    user_theta=float(input('Please confirm the value of the rotation angle tetha: '))
    if user_theta != theta_deg:
        sys.exit("The rotation angle you provided is not the same as given in the input, please check your data!")
    else:
        print("All angles values are consistent, proceeding with the input file creation.")
        
################### CH COORDINATES PROJECTION ONTO GRAVITY PROFILE ################### 
#Calculating projected coordinates on the gravity profile for the gravity stations
mp_gp=proj_gp(mp_xch , mp_ych, gp_xstart, gp_ystart, alpha)

#Calculating projected coordinates on the gravity profile for the drillings
drillings_gp=proj_gp(drillings_xch , drillings_ych, gp_xstart, gp_ystart, alpha)

#ADD TO BA_PLOT
#Calculating projected coordinates on the gravity profile for the regional gravity stations
#regio_mp_gp=proj_gp(regio_mp_xch , regio_mp_ych, gp_xstart, gp_ystart, alpha)



################### CH COORDINATES ROTATION TO THE LCS ################### 

#Gravity stations
mp_lcs=ch_to_lcs(mp_xch, mp_ych, mp_zch, LCS_xori, LCS_yori, zmin, theta)

#Drillings
drillings_lcs=ch_to_lcs(drillings_xch, drillings_ych, drillings_zch, LCS_xori, LCS_yori, zmin, theta)

#Gravity profile
gp_lcs=ch_to_lcs(gp_x_ch, gp_y_ch, gp_z_ch, LCS_xori, LCS_yori, zmin, theta)

#Place x and y in arrays for the intersection calculations
gp_x_lcs=gp_lcs[0]
gp_y_lcs=gp_lcs[1]

mp_xlcs=mp_lcs[0]
mp_ylcs=mp_lcs[1]
mp_zlcs=mp_lcs[2]

drillings_xlcs=drillings_lcs[0]
drillings_ylcs=drillings_lcs[1]
drillings_zlcs=drillings_lcs[2]

################### PLOTTING THE LCS DATA ON A 2D MAP ####################
matplotlib.rcParams.update({'font.size': 20})
hfont = {'fontname':'DejaVu Sans'}
plt.rc('font', weight='bold') 

convertkm=1/1000

prism_x1_2D=prism_x1*convertkm
prism_x2_2D=prism_x2*convertkm
prism_y1_2D=prism_y1*convertkm
prism_y2_2D=prism_y2*convertkm

mp_xlcs_2D=mp_xlcs*convertkm
mp_ylcs_2D=mp_ylcs*convertkm

drillings_xlcs_2D=drillings_xlcs*convertkm
drillings_ylcs_2D=drillings_ylcs*convertkm

gp_lcs=np.asarray(gp_lcs)

gp_lcs_2D=gp_lcs*convertkm

botleftcorners=np.zeros((nb_prisms, 2))
topleftcorners=np.zeros((nb_prisms, 2))
toprightcorners=np.zeros((nb_prisms, 2))
botrightcorners=np.zeros((nb_prisms, 2))

#Prepraring the prisms LCS corners data
all_pri_corners_lcs=np.zeros((nb_prisms,4, 2))

for i in range(nb_prisms):
    botleftcorners[i],topleftcorners[i],toprightcorners[i],botrightcorners[i]=pri_corners(prism_x1[i], prism_x2[i], prism_y1[i], prism_y2[i])

for i in range(nb_prisms):
    all_pri_corners_lcs[i,0][:2]=botleftcorners[i]
    all_pri_corners_lcs[i,1][:2]=topleftcorners[i]
    all_pri_corners_lcs[i,2][:2]=toprightcorners[i]
    all_pri_corners_lcs[i,3][:2]=botrightcorners[i]

all_pri_corners_lcs_2D=all_pri_corners_lcs*convertkm
vertices_LCS=np.zeros((nb_prisms, 5, 2))

for i in show_prism:
    vertices_LCS[i]=prism_verts(all_pri_corners_lcs_2D[i,:,0],all_pri_corners_lcs_2D[i,:,1])

for i in show_prism:
    vertices_top_pri = [
        (prism_x1_2D[0], prism_y1_2D[0]), # left, bottom
        (prism_x1_2D[0], prism_y2_2D[0]), # left, top
        (prism_x2_2D[0], prism_y2_2D[0]), # right, top
        (prism_x2_2D[0], prism_y1_2D[0]), # right, bottom
        (prism_x1_2D[0], prism_y1_2D[0]), # ignored
        ]

codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]


fig1,ax1=plt.subplots()

ax1.scatter(mp_xlcs_2D, mp_ylcs_2D, c='#0080ff',s=100, alpha=1, label='Gravity station')

ax1.scatter(0,0,marker='P', s=200, zorder=1, linewidth=2, edgecolor="k" ,c="w", label='LCS origin')

ax1.scatter(drillings_xlcs_2D, drillings_ylcs_2D,marker='d', s=50, color='#db1e2a', label='Drilling')

ax1.plot(gp_lcs_2D[0,:],gp_lcs_2D[1,:],'-', linewidth=5, color='#fbd409', zorder=0, label='Gravity profile')

ax1.set_xlabel('Distance along the $x_{LCS}$ axis (km)',labelpad=15,fontsize=18)
ax1.set_ylabel('Distance along the $y_{LCS}$ axis (km)',labelpad=15,fontsize=18)

ax1.scatter(gp_lcs_2D[0,0],gp_lcs_2D[1,0], s=300, marker='*', color='#fbd40c', edgecolors='k', label='Profile start')
ax1.scatter(gp_lcs_2D[0,1],gp_lcs_2D[1,1], s=150, marker='o', color='#fbd40c', edgecolors='k', label='Profile end')

#path = Path(vertices_top_pri, codes)
#patch = patches.PathPatch(path, facecolor='none', lw=2, label='Top prism')

for i in show_prism:
    paths= Path(vertices_LCS[i], codes)
    patch_i = patches.PathPatch(paths, edgecolor='k', facecolor='none', lw=1, label='Prism '+str(i+1))
    ax1.add_patch(patch_i)

ax1.set_xlim(xlimits_2D[0],xlimits_2D[1])
ax1.set_ylim(ylimits_2D[0],ylimits_2D[1])
ax1.set_aspect('equal')
ax1.legend(title=str(model_name), fontsize=14, title_fontsize=16, loc="lower right")
ax1.invert_xaxis()

plt.show()


#### Zoom of the center of the LCS, to check the angle between the gravity profile and the x axis of the LCS


fig2,ax2=plt.subplots()

#Plot the maximum anomaly
if max_ano_mp == 'auto':
    mp_max_ano=np.argmin(residual)
else:
    mp_max_ano=max_ano_mp

ax2.scatter(mp_xlcs_2D[mp_max_ano], mp_ylcs_2D[mp_max_ano], c='#0080ff',s=200, alpha=1, zorder=1, label='Max. anomaly station')

if max_ano_mp == "auto":
    ax2.annotate(str('%.f' % MP_num[mp_max_ano]),xy=[mp_xlcs_2D[mp_max_ano], mp_ylcs_2D[mp_max_ano]], size='15')
else:
    for i in max_ano_mp:
        ax2.annotate(str('%.f' % MP_num[i]),xy=[mp_xlcs_2D[i], mp_ylcs_2D[i]], size='15')
    

#Plot the origin of the LCS
ax2.scatter(0,0,marker='P', s=400, linewidth=2, edgecolor="k" ,c="w", zorder=1, label='LCS origin')

#Plot the gravity profile
ax2.plot(gp_lcs_2D[0,:],gp_lcs_2D[1,:],'-', linewidth=5, color='#fbd409', zorder=0, label='Gravity profile')

#plot the x and y axes
#ax2.plot([0,0],[-2,2], color='#db1e2a', linewidth=4, zorder=0, label='LCS axes')
#ax2.plot([-2,2],[0,0], color='#db1e2a', linewidth=4, zorder=0)

#adds arrows to the axes
ax2.arrow(-0.2,0.0, 0.38, 0.0, overhang=0.3, head_width=0.01, color='#db1e2a', linewidth=4, zorder=0, label='LCS axes')
ax2.arrow(0.0,-0.2, 0.0,0.38, overhang=0.3, head_width=0.01, color='#db1e2a', linewidth=4, zorder=0)

ax2.set_xlabel('Distance along the $x_{LCS}$ axis (km)',labelpad=15,fontsize=18)
ax2.set_ylabel('Distance along the $y_{LCS}$ axis (km)',labelpad=15,fontsize=18)
ax2.set_xlim(-0.2,0.2)
ax2.set_ylim(-0.2,0.2)
ax2.set_aspect('equal')
ax2.legend(title=str(model_name)+' LCS origin location' , fontsize=14, title_fontsize=16, loc="lower right")
ax2.invert_xaxis()

plt.show()


########### WRITING THE INPUT FILE ################### 

#Create or open file to store the back rotated data
#Name is given by the user's input for the variable "model_folder"
file_prisma=open(str(run_folder_path)+'/'+str(model_name)+'_PRISM_input.txt','w')
np.set_printoptions(formatter={'float': '{: .4f}'.format})
run_day_time=datetime.datetime.now()
#Write the gravity stations coordinates
print("#Input for "+str(model_name), file=file_prisma)
print("", file=file_prisma)

print("#Input creation date and time:", file=file_prisma)
print(str("%s"%run_day_time), file=file_prisma)
print("", file=file_prisma)

print("#Initial profile data input file", file=file_prisma)
print(str(initial_input_file), file=file_prisma)
print("", file=file_prisma)

print("#Angles", file=file_prisma)
print("#For the projection onto the gravity profile (°), alpha", file=file_prisma)
print(alpha_deg, file=file_prisma)
print("", file=file_prisma)

print("#For LCS rotation and back rotation to CH coor (°), theta", file=file_prisma)
print(theta_deg, file=file_prisma)
print("", file=file_prisma)

print("#Gravity profile length (m)", file=file_prisma)
print(gp_length, file=file_prisma)
print("", file=file_prisma)

print("#Gravity profile start CH-coordinates and LCS coordinates", file=file_prisma)
print(str(gp_x_ch[0]), str(gp_y_ch[0]), str(gp_z_ch[0]),file=file_prisma)
print(str("{0:<10.3f}".format(gp_x_lcs[0])), str("{0:<10.3f}".format(gp_y_lcs[0])), str(gp_lcs[2][0]),file=file_prisma)
print("", file=file_prisma)

print("#Gravity profile end CH-coordinates and LCS coordinates", file=file_prisma)
print(str(gp_x_ch[1]), str(gp_y_ch[1]), str(gp_z_ch[1]),file=file_prisma)
print(str("{0:<10.3f}".format(gp_x_lcs[1])), str("{0:<10.3f}".format(gp_y_lcs[1])), str(gp_lcs[2][1]),file=file_prisma)
print("", file=file_prisma)

print("#LCS origin CH-coordinates", file=file_prisma)
print(str(LCS_xori), str(LCS_yori), str(zmin),file=file_prisma)
print("", file=file_prisma)

print("#Number of stations", file=file_prisma)
print(nb_stations, file=file_prisma)
print("", file=file_prisma)

print("#Number of drillings", file=file_prisma)
print(nb_drillings, file=file_prisma)
print("", file=file_prisma)

#DISPLAY BA ROUTINE
#print("#Regional gravity stations input file", file=file_prisma)
#print("#"+str(regio_MP_file_path), file=file_prisma)
#print("", file=file_prisma)
#DISPLAY BA ROUTINE
#print("#Number of stations", file=file_prisma)
#print(nb_regio_stations, file=file_prisma)
#print("", file=file_prisma)

#### NEED TO ADD BOUGUER ANOMALY Bouguer anomaly (mGal)\t| #####
#### NEED TO ADD SINGLE MP/DRILLINGS/PRISM case ####
if check_dist=="yes" or check_dist=="YES" or check_dist=="Yes":
    
    MP_dist_check=[]
    drill_dist_check=[]
    
#Gravity stations data    
    print("#MP Station \t| CH coordinates(m) E,\tN,\tElev\t| Distance(m) along\taway from the profile\t| Distance check \t|" +
          "LCS coordinates(m) x,\ty,\tz\t| Bouguer anomaly(mGal)\t| Regional anomaly(mGal)\t| Residual anomaly(mGal)|", file=file_prisma)
    for i in range(nb_stations):
        
        #Check if stations projected along the gravity profile are farther than the distance given by the user, from the profile
        if np.abs(mp_gp[1][i]) > flag_distance:
            #If the drilling is over this distance, flag it as far
            MP_dist_check.append("far")
        else:
            #Otherwise mark it as ok
            MP_dist_check.append("ok")

              #Original input coordinates
        print(str("{0:<10.0f}".format(MP_num[i])),"\t",str("{0:<10.3f}".format(mp_xch[i])),"\t",str("{0:<10.3f}".format(mp_ych[i])),"\t",str("{0:<10.3f}".format(mp_zch[i])),"\t",
              #Projected input coordinates on the gravity profile and distance check results
              str("{0:<10.2f}".format(mp_gp[0][i])),"\t",str("{0:<10.2f}".format(mp_gp[1][i])),"\t", MP_dist_check[i],"\t",
              #LCS coordinates
              str("{0:<10.2f}".format(mp_xlcs[i])),"\t",str("{0:<10.2f}".format(mp_ylcs[i])),"\t",str("{0:<10.2f}".format(mp_zlcs[i])),"\t",
              #Bouguer (TBA) and residual anomaly values
              str("{0:<10.2f}".format(bouguer[i])),"\t",str("{0:<10.2f}".format(regional[i])),"\t",str("{0:<10.2f}".format(residual[i])),file=file_prisma)
    print(" ", file=file_prisma)

#Drillings data    
    print("#Drilling number\t| CH coordinates(m) E,\tN,\tElev\t| Distance(m) along\taway from the profile\t| Distance check \t| LCS coordinates(m) x,\ty,\tz", file=file_prisma)
    for i in range(len(drillings[0])):
        
        #Check if drillings projected along the gravity profile are farther than the distance given by the user, from the profile
        if np.abs(drillings_gp[1][i]) > flag_distance:
            #If the drilling is over this distance, flag it as far
            drill_dist_check.append("far")
        else:
            #Otherwise mark it as ok
            drill_dist_check.append("ok")
            
              #Drillings number
        print(str("{0:<10.0f}".format(drillings[0,i])),"\t",
              #Original input coordinates
              str("{0:<10.2f}".format(drillings_xch[i])),"\t",str("{0:<10.2f}".format(drillings_ych[i])),"\t",str("{0:<10.2f}".format(drillings_zch[i])),"\t",
              #Projected input coordinates on the gravity profile and distance check results
              str("{0:<10.2f}".format(drillings_gp[0][i])),"\t",str("{0:<10.2f}".format(drillings_gp[1][i])),"\t", drill_dist_check[i],"\t",
              #LCS coordinates
              str("{0:<10.2f}".format(drillings_xlcs[i])),"\t",str("{0:<10.2f}".format(drillings_ylcs[i])),"\t",str("{0:<10.2f}".format(drillings_zlcs[i])),file=file_prisma)
    print(" ", file=file_prisma)
    
#Prisms data   
    print("#Prism LCS coordinates, input for PRISMA\t" ,file=file_prisma)
    print("#Number of prisms", file=file_prisma)
    print(nb_prisms, file=file_prisma)
    print("", file=file_prisma)
    #Top prism corner CH coordinates
    #Used to check the LCS angle input in the Display routine
    print("#Top prism CH coordinates", file=file_prisma)
    print("#E\tN\t\tElevation", file=file_prisma)
    for i in range(4):
        print(str("{0:<10.3f}".format(pri_ch_E[i])),"\t",str("{0:<10.3f}".format(pri_ch_N[i])),"\t",str("{0:<10.3f}".format(pri_ch_elev[i])), file=file_prisma)
    print("", file=file_prisma)
    #Input LCS coordinates for the prism
    #These coordinates will be used by PRISMA for the calculations 
    print("#xmax(m)\t xmin(m)\t ymax(m)\t ymin(m)\t zmax(m)\t zmin(m)\t density contrast(kg/m3)",file=file_prisma)
    for i in range(nb_prisms):
        print(str("{0:<10.0f}".format(prism_x2[i])),"\t",str("{0:<10.0f}".format(prism_x1[i])),"\t",
              str("{0:<10.0f}".format(prism_y2[i])),"\t",str("{0:<10.0f}".format(prism_y1[i])),"\t",
              str("{0:<10.2f}".format(prism_z2[i])),"\t",str("{0:<10.2f}".format(prism_z1[i])),"\t",
              str("{0:<10.0f}".format(rho_pri[i])),"\t",file=file_prisma)
    print(" ", file=file_prisma)

#### THIS SECTION NEEDS TO BE UPDATED ####
else:
    print("#MP Station \t| Back rotated CH coordinates E,\tN,\tz\t| Distance along gravity profile\t| Residual anomaly (mGal)\t| Modelled residual anomaly (mGal)", file=file_prisma)
    for i in range(nb_stations):
#Gravity stations data
              #Original input coordinates
        print(str("{0:<10.0f}".format(MP_num[i])),"\t",str("{0:<10.3f}".format(mp_xch[i])),"\t",str("{0:<10.3f}".format(mp_ych[i])),"\t",str("{0:<10.3f}".format(mp_zch[i])),"\t",
              #Projected input coordinates on the gravity profile
              str("{0:<10.2f}".format(mp_gp[i])),"\t",
              #LCS coordinates
              str("{0:<10.2f}".format(mp_xlcs[i])),"\t",str("{0:<10.2f}".format(mp_ylcs[i])),"\t",str("{0:<10.2f}".format(mp_zlcs[i])),"\t",
              #Bouguer (TBA) and residual anomaly values
              str("{0:<10.2f}".format(residual[1,i])),file=file_prisma) #,"\t",str("{0:<10.2f}".format(prism_results[4,i]))
    print(" ", file=file_prisma)    

#Write the drillings coordinates
    print("#Drilling number\t| Back rotated CH coordinates x,\ty,\tz\t| Distance along gravity profile", file=file_prisma)
    for i in range(len(drillings[0])):
              #Drillings number
        print(str("{0:<10.0f}".format(drillings[0,i])),"\t",
              #Original input coordinates
              str("{0:<10.3f}".format(drillings_xch[i])),"\t",str("{0:<10.3f}".format(drillings_ych[i])),"\t",str("{0:<10.3f}".format(drillings_zch[i])),"\t",
              #Projected input coordinates on the gravity profile
              str("{0:<10.2f}".format(drillings_gp[i])),"\t",
              #LCS coordinates
              str("{0:<10.2f}".format(drillings_xlcs[i])),"\t",str("{0:<10.2f}".format(drillings_ylcs[i])),"\t",str("{0:<10.2f}".format(drillings_zlcs[i])),file=file_prisma)
    print(" ", file=file_prisma)

    
file_prisma.close()












