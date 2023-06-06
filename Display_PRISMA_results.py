# -*- coding: utf-8 -*-
"""
Created on 13/04/2021

@author: Dimitri
"""

import numpy as np
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as patches
import itertools
import shutil
import random
from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.path import Path
from matplotlib import gridspec


############ INPUT PROCESSING ##########################

#DEM and bedrock profile import
DEM_cross_data=distance, DEM_E, DEM_N, DEM_elev=np.loadtxt(str(DEM_file), unpack=True)#Elevation profile data
bedrock_cross_data=distance_b, bedrock_E, bedrock_N, bedrock_elev=np.loadtxt(str(bedrock_file), unpack=True)#Bedrock model profile data

PRISMA_output_file="Display/"+str(model_name)+"_PRISM_output.txt"

#Fixed length of the header before the MPs input
length_header=48

with open(PRISMA_output_file) as file_input:
   nb_file_input=file_input.readlines()
#Load the number of stations from the file, this position is the line number in the file and is fixed
#Used to extract the input data and to check the user's input
   nb_stations=int(nb_file_input[42])
   
#Load the number of drillings from the file, this position is the line number in the file and is fixed
#Only needed to extract the input data. Drillings values are not used in this code
   nb_drillings=int(nb_file_input[45])
   
#Load the number of prisms from the file, this position is the line number in the file and is fixed
#Used to extract the input data and to check the user's input
   nb_prisms=int(nb_file_input[length_header + nb_stations + 2 + nb_drillings + 3])

   alpha=float(nb_file_input[22])
   
   theta=float(nb_file_input[25])
   
#Load the profile length value from the file, this position is the line number in the file and is fixed
   profile_length=float(nb_file_input[28])
   
   gp_chstart=np.fromstring(nb_file_input[31], dtype='float', sep=' ')
   gp_lcsstart=np.fromstring(nb_file_input[32], dtype='float', sep=' ')
   gp_chend=np.fromstring(nb_file_input[35], dtype='float', sep=' ')
   gp_lcsend=np.fromstring(nb_file_input[36], dtype='float', sep=' ')
   LCS_ori=np.fromstring(nb_file_input[39], dtype='float', sep=' ')

#Assign variables to coordinates
gp_xstart=gp_chstart[0]
gp_ystart=gp_chstart[1]
LCS_xori=LCS_ori[0]
LCS_yori=LCS_ori[1]
zmin=LCS_ori[2]



#Getting the MPs data
#Skipheader is the fixed number of lines at the start of the file
MPs_data_header=length_header
#Skipfooter is the number of drillings, the top prism corner coordinates and the number of prisms
MPs_data_footer=nb_drillings + 4 + nb_prisms
#Load the stations coordinates
MPs_all_data =  np.genfromtxt(PRISMA_output_file, skip_header=MPs_data_header, skip_footer=MPs_data_footer, unpack=True)
#Load the distance check results separately, as they are string data and are loaded above as "NaN"
MP_distance_check = np.genfromtxt(PRISMA_output_file, delimiter='\t', skip_header=MPs_data_header, skip_footer=MPs_data_footer, usecols=6, dtype=str)

#Getting the drillings data
#Skipheader is the fixed number of lines at the start of the file
drillings_data_header= length_header + nb_stations + 2
#Skipfooter is the top prism corner coordinates and the number of prisms
drillings_data_footer=4 + nb_prisms 
#Load the drillings coordinates
drillings_data = np.genfromtxt(PRISMA_output_file, skip_header=drillings_data_header, skip_footer=drillings_data_footer, unpack=True)
#Load the distance check results separately, as they are string data and are loaded above as "NaN"
drillings_distance_check = np.genfromtxt(PRISMA_output_file, delimiter='\t', skip_header=drillings_data_header, skip_footer=drillings_data_footer, usecols=6, dtype=str)

#Getting the prism data
#Skipheader is the fixed number of lines at the start of the file and the MPs,the drillings and their header and the prisms coordinates header and prism number
prism_ch_data_header= drillings_data_header + nb_drillings +  7
#Skipfooter is the number of prisms
prism_ch_data_footer=nb_prisms
#Load the top prism corner CH coordinates
pri_ch_data=np.genfromtxt(PRISMA_output_file, skip_header=prism_ch_data_header, skip_footer=prism_ch_data_footer, unpack=True)

#Getting the prism data
#Skipheader is the fixed number of lines at the start of the file and the MPs,the drillings and their header and the prism corner coordinates and the prisms header and prisms number
prism_data_header= prism_ch_data_header + 6 
#No data to skip at the footer
prism_data_footer=0
#Load the prisms LCS coordinates and density values
pri_lcs_data = np.genfromtxt(PRISMA_output_file, skip_header=prism_data_header, skip_footer=prism_data_footer, unpack=True)


np.set_printoptions(suppress=True)
check_dist="yes"

#Assign data to array variables

MPs_plot_data=MPs_all_data[:14]

MP_num, mp_xch, mp_ych, mp_elev, Dymp, Dxmp, mp_nan, mp_xlcs, mp_ylcs, mp_zlcs, bouguer, regional, residual, model_results = MPs_plot_data
drill_num, drillxch, drillych, drillelev, Dydrill, Dxdrill, drill_nan, drillings_xlcs, drillings_ylcs, drillings_zlcs = drillings_data
pri_ch_E, pri_ch_N, pri_ch_elev= pri_ch_data
prism_x2, prism_x1, prism_y2,  prism_y1, prism_z2, prism_z1, rho_pri= pri_lcs_data


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
    if check_dist=="yes" or check_dist=="YES" or check_dist=="Yes":
        dist_xgp=x_rotation(xlcs1, ylcs1, alpha)
        return(dist_ygp,dist_xgp)
    elif check_dist=="yes" or check_dist=="YES" or check_dist=="Yes":
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

def backro_pri_co(cornerbl, cornertl, cornertr, cornerbr, prism_z1, prism_z2, LCS_xori, LCS_yori, zmin, theta):
    cornerbl_b_ch=lcs_to_ch(cornerbl[0], cornerbl[1], prism_z1, LCS_xori, LCS_yori, zmin, theta)
    cornertl_b_ch=lcs_to_ch(cornertl[0], cornertl[1], prism_z2, LCS_xori, LCS_yori, zmin, theta) 
    cornertr_b_ch=lcs_to_ch(cornertr[0], cornertr[1], prism_z1, LCS_xori, LCS_yori, zmin, theta)
    cornerbr_b_ch=lcs_to_ch(cornerbr[0], cornerbr[1], prism_z2, LCS_xori, LCS_yori, zmin, theta) 
    return(cornerbl_b_ch, cornertl_b_ch, cornertr_b_ch, cornerbr_b_ch)

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

########### PROJECTION OF PROFILE DATA - CHECKING PROJECTION ANGLE VALUE ###################
#Projecting the bedrock profile data to check the projection angle given by the user
#The back projected distance value will then be compared with the existing value

#Projecting the data long the gravity profile
bedrock_b_gp=proj_gp(bedrock_E , bedrock_N, gp_xstart, gp_ystart, alpha)

#Comparing the distance values
if check_dist=="yes" or "YES" or "Yes":
    check_proj=distance_b-bedrock_b_gp[0]
else:
    check_proj=distance_b-bedrock_b_gp

#Check if any value has a difference over 0.5 m or 50 cm, which is the uncertainty of the Swiss 2m DEM
if np.abs(check_proj).all() < 0.5 :
    print(" ")
    print("The projection angle alpha is ok")
else:
    print(" ")
    print("WARNING:The projection angle alpha is not correct and the data will not be projected correctly, please check your alpha value")

########### CALCULATION PRISM/GRAVITY PROFILE INTERSECTION ################### 
#Place x and y in arrays for the intersection calculations
gp_lcs=np.zeros([2,2])
gp_lcs[0,0]=gp_lcsstart[0]
gp_lcs[0,1]=gp_lcsend[0]
gp_lcs[1,0]=gp_lcsstart[1]
gp_lcs[1,1]=gp_lcsend[1]

gp_x_lcs=gp_lcs[0]
gp_y_lcs=gp_lcs[1]

#Get linear equations for the intersection points
a_gp=get_a_line_equation(gp_x_lcs,gp_y_lcs)
b_gp=get_b_line_equation(a_gp,gp_x_lcs[0],gp_y_lcs[0])

#Create arrays to store the intersection points
prism_y1_inter=np.zeros(nb_prisms)
prism_y2_inter=np.zeros(nb_prisms)

##### Prisms calculation of the intersection prisms and gravity profile #####


#Getting intersection coordinates prism/gp and getting their LCS coordinates rotated for the projection on the gravity profile
for i in range(nb_prisms):
    prism_y1_inter[i]=get_y_intersection(a_gp,b_gp,prism_x1[i])
    prism_y2_inter[i]=get_y_intersection(a_gp,b_gp,prism_x2[i])

#Search for prisms that do not cross the gravity profile, and exclude them from the cross-section plots
#based on the user input
if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    print(' ')
    print('Only the prisms crossing the gravity profile will be shown on the cross-section plots')
    #Create arrays to store prism x coordinates that cross the gravity profile
    prism_x1_cross=prism_x1
    prism_x2_cross=prism_x2
    prism_z1_cross=prism_z1
    prism_z2_cross=prism_z2
    #Number of prisms that cross the gravity profile
    nb_prisms_cross=nb_prisms
    prisms_exclu=np.array([])
    
    #Check if the y coordinates of a prism crosses at least the lowest position of the gravity profile (gravity profile starting point)
    for i in range(nb_prisms):
        if prism_y1[i] < gp_y_lcs[0]  and prism_y2[i] > gp_y_lcs[0] :
            prisms_exclu=prisms_exclu
        else:
            prisms_exclu=np.append(prisms_exclu, i)
    if len(prisms_exclu) !=0:
        print(' ')
        print("Prism number "+str(prisms_exclu)+" isn''t crossing the gravity profile, removing it from the plot")
        prism_x1_cross=np.delete(prism_x1_cross, prisms_exclu)
        prism_x2_cross=np.delete(prism_x2_cross, prisms_exclu)
        prism_z1_cross=np.delete(prism_z1_cross, prisms_exclu)
        prism_z2_cross=np.delete(prism_z2_cross, prisms_exclu)
    else:
        prism_x1_cross=prism_x1
        prism_x2_cross=prism_x2
        prism_z1_cross=prism_z1
        prism_z2_cross=prism_z2
        print(' ')
        print("All the prisms are crossing the gravity profile")
    
    #Calculating the intersection with the prisms x coordinates left
    nb_prisms_cross=nb_prisms_cross-len(prisms_exclu)
    prism_y1_inter_cross=np.zeros(nb_prisms_cross)
    prism_y2_inter_cross=np.zeros(nb_prisms_cross)
    print(' ')
    print('Removed '+str(nb_prisms-nb_prisms_cross)+' prisms for the cross-section plot')
    for i in range(nb_prisms_cross):
        prism_y1_inter_cross[i]=get_y_intersection(a_gp,b_gp,prism_x1_cross[i])
        prism_y2_inter_cross[i]=get_y_intersection(a_gp,b_gp,prism_x2_cross[i])
else:
    print('All the prisms will be shown on the cross-section plots')
    
########### CALCULATIONS PROJECTION PERPENDICULAR LCS Y AXIS ###############
#In case the angle between the LCS x-axis and the gravity profile is over 5°, project the MPS perpendicularly to the LCS y axis, onto the gravity profile
#Then project the gravity profile coordinates along the x axis
#Check plot switch
if perpendicular_plot=="yes" or perpendicular_plot=="yes" or perpendicular_plot=="yes" :
#Create arrays to store data crossing gravity profile, perpendicular to LCS y axis
    #Gravity stations
    mp_ylcs_perp=np.zeros(nb_stations)
    #Drillings
    dril_ylcs_perp=np.zeros(nb_drillings)
#Calculate the intersection y coordinates
    for i in range(nb_stations):
        mp_ylcs_perp[i]=get_y_intersection(a_gp,b_gp,mp_xlcs[i])
    for i in range(nb_drillings):
        dril_ylcs_perp[i]=get_y_intersection(a_gp,b_gp,drillings_xlcs[i])
    
#Calculate the LCS coordinates of the DEM and bedrock profiles
    DEM_lcs=ch_to_lcs(DEM_E , DEM_N , DEM_elev, LCS_xori, LCS_yori, zmin, theta)
    DEM_xlcs=DEM_lcs[0]
    DEM_ylcs=DEM_lcs[1]
    
    bedrock_lcs=ch_to_lcs(bedrock_E , bedrock_N , bedrock_elev, LCS_xori, LCS_yori, zmin, theta)
    bedrock_xlcs=bedrock_lcs[0]
    bedrock_ylcs=bedrock_lcs[1]
    
#Project data onto the LCS x axis
    #Calculate the angle between the LCS x axis and the gravity profile
    angle_diff=np.abs(theta_deg)+np.abs(alpha_deg)-90
    #Assign the proper sign to the angle, by checking the x lcs coordinate of the gravity profile start
        #Following the convention of the start of the profile always being the most toward the South
        #If the gravity profile start x is negative, the angle will always be clockwise, therefore negative
        #If the start has a positive x, the angle will always be counterclockwise, therefore positive 
    if gp_lcsstart[0] < 0:
        gamma=angle_diff *-1
    else:
        gamma=angle_diff
        
    gamma_deg=gamma
    gamma=rad_conv(gamma)
    
    #Gravity stations
    mp_perp=proj_gp(mp_xlcs,mp_ylcs_perp, gp_x_lcs[0], gp_y_lcs[0], gamma)
    mp_x_perp=mp_perp[1]
    mp_y_perp=mp_perp[0]
    
    #Intermediate step of the projection to check the results
    mp_x_perp_check=x_rotation(mp_xlcs,mp_ylcs_perp, gamma)
    mp_y_perp_check=y_rotation(mp_xlcs,mp_ylcs_perp, gamma)
    
    #Drillings
    drill_perp=proj_gp(drillings_xlcs,dril_ylcs_perp, gp_x_lcs[0], gp_y_lcs[0], gamma)
    drill_x_perp=drill_perp[1]
    drill_y_perp=drill_perp[0]

    #DEM
    DEM_perp=proj_gp(DEM_xlcs,DEM_ylcs, gp_x_lcs[0], gp_y_lcs[0], gamma)
    DEM_x_perp=DEM_perp[1]
    DEM_y_perp=DEM_perp[0]
    #bedrock
    bedrock_perp=proj_gp(bedrock_xlcs,bedrock_ylcs, gp_x_lcs[0], gp_y_lcs[0], gamma)
    bedrock_x_perp=bedrock_perp[1]
    bedrock_y_perp=bedrock_perp[0]
    
    
################### LCS COORDINATES BACK ROTATION TO THE CH COORDINATE SYSTEM ################### 

#Gravity profile
gp_b_ch=lcs_to_ch(gp_lcs[0], gp_lcs[1], zmin, LCS_xori, LCS_yori, zmin, theta)

#Gravity stations
mp_b_ch=lcs_to_ch(mp_xlcs, mp_ylcs, mp_zlcs, LCS_xori, LCS_yori, zmin, theta)    

#Drillings
drillings_b_ch=lcs_to_ch(drillings_xlcs, drillings_ylcs, drillings_zlcs, LCS_xori, LCS_yori, zmin, theta)

#Top prism corners
tpc1_b_ch=lcs_to_ch(prism_x2[0], prism_y2[0], prism_z1, LCS_xori, LCS_yori, zmin, theta)
tpc2_b_ch=lcs_to_ch(prism_x1[0], prism_y2[0], prism_z2, LCS_xori, LCS_yori, zmin, theta) 
tpc3_b_ch=lcs_to_ch(prism_x2[0], prism_y1[0], prism_z1, LCS_xori, LCS_yori, zmin, theta)
tpc4_b_ch=lcs_to_ch(prism_x1[0], prism_y1[0], prism_z2, LCS_xori, LCS_yori, zmin, theta) 

#Prisms
all_pri_corners_ch=np.zeros((nb_prisms,4, 3))

botleftcorners=np.zeros((nb_prisms, 2))
topleftcorners=np.zeros((nb_prisms, 2))
toprightcorners=np.zeros((nb_prisms, 2))
botrightcorners=np.zeros((nb_prisms, 2))

for i in range(nb_prisms):
    botleftcorners[i],topleftcorners[i],toprightcorners[i],botrightcorners[i]=pri_corners(prism_x1[i], prism_x2[i], prism_y1[i], prism_y2[i])
    all_pri_corners_ch[i]=backro_pri_co(botleftcorners[i],topleftcorners[i],toprightcorners[i],botrightcorners[i],
                                            prism_z1[i], prism_z2[i], LCS_xori, LCS_yori, zmin, theta)    

#Intersection prisms/gravity profile
prism1_int_b_ch=lcs_to_ch(prism_x1, prism_y1_inter, prism_z1, LCS_xori, LCS_yori, zmin, theta) 
prism2_int_b_ch=lcs_to_ch(prism_x2, prism_y2_inter, prism_z2, LCS_xori, LCS_yori, zmin, theta) 

if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    prism1_int_b_ch_cross=lcs_to_ch(prism_x1_cross, prism_y1_inter_cross, prism_z1_cross, LCS_xori, LCS_yori, zmin, theta) 
    prism2_int_b_ch_cross=lcs_to_ch(prism_x2_cross, prism_y2_inter_cross, prism_z2_cross, LCS_xori, LCS_yori, zmin, theta) 

#Assigning results to variables
mp_b_xch=mp_b_ch[0]
mp_b_ych=mp_b_ch[1]
mp_b_zch=mp_b_ch[2]

drillings_b_xch=drillings_b_ch[0]
drillings_b_ych=drillings_b_ch[1]
drillings_b_zch=drillings_b_ch[2]

tpc_b_x1ch=tpc1_b_ch[0]
tpc_b_x2ch=tpc2_b_ch[0]
tpc_b_x3ch=tpc3_b_ch[0]
tpc_b_x4ch=tpc4_b_ch[0]
tpc_b_y1ch=tpc1_b_ch[1]
tpc_b_y2ch=tpc2_b_ch[1]
tpc_b_y3ch=tpc3_b_ch[1]
tpc_b_y4ch=tpc4_b_ch[1]

prism_int_b_x1ch=prism1_int_b_ch[0]
prism_int_b_x2ch=prism2_int_b_ch[0]
prism_int_b_y1ch=prism1_int_b_ch[1]
prism_int_b_y2ch=prism2_int_b_ch[1]
prism_int_b_z1ch=prism1_int_b_ch[2]
prism_int_b_z2ch=prism2_int_b_ch[2]

if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    prism_int_b_x1ch_cross=prism1_int_b_ch_cross[0]
    prism_int_b_x2ch_cross=prism2_int_b_ch_cross[0]
    prism_int_b_y1ch_cross=prism1_int_b_ch_cross[1]
    prism_int_b_y2ch_cross=prism2_int_b_ch_cross[1]
    prism_int_b_z1ch_cross=prism1_int_b_ch_cross[2]
    prism_int_b_z2ch_cross=prism2_int_b_ch_cross[2]

################### BACK ROTATED CH COORDINATES PROJECTION ONTO GRAVITY PROFILE ################### 

#Gravity stations
mp_b_gp=proj_gp(mp_b_xch, mp_b_ych, gp_xstart, gp_ystart, alpha)

#Drillings
drillings_b_gp=proj_gp(drillings_b_xch, drillings_b_ych, gp_xstart, gp_ystart, alpha)
 
#Intersection prisms/gravity profile                       
prism_int_b_1gp=proj_gp(prism_int_b_x1ch, prism_int_b_y1ch, gp_xstart, gp_ystart, alpha) 
prism_int_b_2gp=proj_gp(prism_int_b_x2ch, prism_int_b_y2ch, gp_xstart, gp_ystart, alpha) 

#LCS origin position along the profile - Not back rotated coordinates
LCS_proj=proj_gp(LCS_ori[0], LCS_ori[1], gp_xstart, gp_ystart, alpha)

if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    prism_int_b_1gp_cross=proj_gp(prism_int_b_x1ch_cross, prism_int_b_y1ch_cross, gp_xstart, gp_ystart, alpha) 
    prism_int_b_2gp_cross=proj_gp(prism_int_b_x2ch_cross, prism_int_b_y2ch_cross, gp_xstart, gp_ystart, alpha)
    
#Assigning results to variables
prism_int_b_x1gp=prism_int_b_1gp[0]
prism_int_b_x2gp=prism_int_b_2gp[0]
prism_int_b_y1gp=prism_int_b_1gp[1]
prism_int_b_y2gp=prism_int_b_2gp[1]
prism_b_z1gp=prism_int_b_z1ch
prism_b_z2gp=prism_int_b_z2ch

if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    prism_int_b_x1gp_cross=prism_int_b_1gp_cross[0]
    prism_int_b_x2gp_cross=prism_int_b_2gp_cross[0]
    prism_int_b_y1gp_cross=prism_int_b_1gp_cross[1]
    prism_int_b_y2gp_cross=prism_int_b_2gp_cross[1]
    prism_b_z1gp_cross=prism_int_b_z1ch_cross
    prism_b_z2gp_cross=prism_int_b_z2ch_cross
    
################### CHECKING COORDINATES ################
#This section checks the difference in value between back rotated coordinates and the original coordinates
#The difference due to rounding should not be higher than 0.01 m

#Gravity stations
MP_data_check= np.abs(mp_b_ch-MPs_all_data[1:4])

#for each x, y, z coordinates
for j in range(3):
    #for each station
    for i in range(nb_stations):
        #check that for values above 1 cm
        if np.any(MP_data_check[j,i] > 0.01)==True:
            print("Please check your angle for the back rotation, station "+str(MP_num[i])+" back rotated coordinates are incorrect")
            
#Drillings
drill_data_check= np.abs(drillings_b_ch-drillings_data[1:4])

#for each x, y, z coordinates
for j in range(3):
    #for each station
    for i in range(nb_drillings):
        #check that for values above 1 cm
        if np.any(drill_data_check[j,i] > 0.01)==True:
            print("Please check your angle for the back rotation, drilling "+str(drill_num[i])+" back rotated coordinates are incorrect")

            
################### PLOTTING ROUTINE ################### 
#Set font size
matplotlib.rcParams.update({'font.size': 14})
#Set the font for the annotations
hfont = {'fontname':'Arial'}
#Set font weight to bold
plt.rc('font', weight='bold')

##### 
fig3 = plt.figure(constrained_layout=True)
widths = [3, 2]
heights = [3, 1.5]
gs = fig3.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
                          height_ratios=heights)
f3_ax1 = fig3.add_subplot(gs[0, :])
f3_ax2 = fig3.add_subplot(gs[1, :])

############## RESIDUAL ANOMALY/MODELLED ANOMALY PLOT ##################

#############Top plot###############

#Observed residual anomaly
f3_ax1.scatter(Dymp,residual, color='#0081f1', s=300,zorder=1, label="Observed residual anomaly")
f3_ax1.errorbar(Dymp,residual, yerr=0.13, color='k', capsize=5,linestyle='None',zorder=0)

#Modelled residual anomaly
f3_ax1.scatter(Dymp,model_results, color='#df7116', s=300,zorder=1, label="Modelled residual anomaly")


f3_ax1.plot(xlimits_plot, [0,0], color='k', linestyle="--", linewidth=2, dashes=[10,3,10, 3])
    
#Text for orientation
f3_ax1.text(0,left_side, '%s' %str(left_label),horizontalalignment='center',verticalalignment='top', size= 30, fontweight='bold',zorder=1)
f3_ax1.text(profile_length,right_side, '%s' %str(right_label),horizontalalignment='center',verticalalignment='top', size= 30, fontweight='bold',zorder=1)

#Write the different density contrast values used for this model

rho=np.unique(rho_pri)


if len(rho)==1:
    f3_ax1.text(xrho,yrho, '%s' %'Δρ(kg/m3)='+str("{0:<10.0f}".format(rho[0])),horizontalalignment='left',verticalalignment='top',fontweight='bold', **hfont, size= 24,zorder=1)
else:
    f3_ax1.text(xrho,yrho, '%s' %'Δρ(kg/m3)= '+','.join(list(map("{0:^10.0f}".format,rho))),horizontalalignment='left',verticalalignment='top',fontweight='bold', **hfont, size= 24,zorder=1)
        

#Write the axes titles
f3_ax1.set_xlabel("Distance along the gravity profile (m)", size=20)
f3_ax1.set_ylabel("Residual Anomaly (mGal)", size=20)
f3_ax1.legend(loc='lower left')



################## Bottom plot ##################

convertkm=1/1000

#Get position of the maximum anomaly value
posmin=np.where(residual==residual.min())

pos_max_ano=LCS_proj[0]

#Subtract coordinate for to have the distance from the maximum anomaly stations
#DEM and bedrock profiles
bot_dist=distance-pos_max_ano
bot_dist_b=distance_b-pos_max_ano
#Prism coordinates
if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    bot_pri_x2_cross=prism_int_b_x2gp_cross-pos_max_ano
    bot_pri_x1_cross=prism_int_b_x1gp_cross-pos_max_ano

else:
    bot_pri_x2=prism_int_b_x2gp-pos_max_ano
    bot_pri_x1=prism_int_b_x1gp-pos_max_ano
#Gravity stations
bot_mp=Dymp-pos_max_ano
#Drillings
bot_drill=Dydrill-pos_max_ano
#x axis limits for the bottom plot
xlimits_plot_bot=xlimits_plot-pos_max_ano

#DEM and bedrock plots
f3_ax2.plot(bot_dist,DEM_elev, color='#004b88',label='Elevation')
f3_ax2.plot(bot_dist_b,bedrock_elev, '--',color='#800000',label='Bedrock model')

#Plot the prisms with an exception for a single prism plot
if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
    if nb_prisms_cross==1:
        f3_ax2.plot([bot_pri_x2_cross,bot_pri_x1_cross],[prism_b_z1gp_cross,prism_b_z1gp_cross], linewidth=2,color='k') #top
        f3_ax2.plot([bot_pri_x1_cross,bot_pri_x1_cross],[prism_b_z1gp_cross,prism_b_z2gp_cross], linewidth=2,color='k') #right
        f3_ax2.plot([bot_pri_x2_cross,bot_pri_x1_cross],[prism_b_z2gp_cross,prism_b_z2gp_cross], linewidth=2,color='k') #bottom
        f3_ax2.plot([bot_pri_x2_cross,bot_pri_x2_cross],[prism_b_z2gp_cross,prism_b_z1gp_cross], linewidth=2,color='k') #left
    else:
        for i in range(nb_prisms_cross):
            f3_ax2.plot([bot_pri_x2_cross[i],bot_pri_x1_cross[i]],[prism_b_z1gp_cross[i],prism_b_z1gp_cross[i]], linewidth=2,color='k') #top
            f3_ax2.plot([bot_pri_x1_cross[i],bot_pri_x1_cross[i]],[prism_b_z1gp_cross[i],prism_b_z2gp_cross[i]], linewidth=2,color='k') #right
            f3_ax2.plot([bot_pri_x2_cross[i],bot_pri_x1_cross[i]],[prism_b_z2gp_cross[i],prism_b_z2gp_cross[i]], linewidth=2,color='k') #bottom
            f3_ax2.plot([bot_pri_x2_cross[i],bot_pri_x2_cross[i]],[prism_b_z2gp_cross[i],prism_b_z1gp_cross[i]], linewidth=2,color='k') #left
else:
    if nb_prisms==1:
        f3_ax2.plot([bot_pri_x2,bot_pri_x1],[prism_b_z1gp,prism_b_z1gp], linewidth=2,color='k') #top
        f3_ax2.plot([bot_pri_x1,bot_pri_x1],[prism_b_z1gp,prism_b_z2gp], linewidth=2,color='k') #right
        f3_ax2.plot([bot_pri_x2,bot_pri_x1],[prism_b_z2gp,prism_b_z2gp], linewidth=2,color='k') #bottom
        f3_ax2.plot([bot_pri_x2,bot_pri_x2],[prism_b_z2gp,prism_b_z1gp], linewidth=2,color='k') #left
    else:
        for i in range(nb_prisms):
            f3_ax2.plot([bot_pri_x2[i],bot_pri_x1[i]],[prism_b_z1gp[i],prism_b_z1gp[i]], linewidth=2,color='k') #top
            f3_ax2.plot([bot_pri_x1[i],bot_pri_x1[i]],[prism_b_z1gp[i],prism_b_z2gp[i]], linewidth=2,color='k') #right
            f3_ax2.plot([bot_pri_x2[i],bot_pri_x1[i]],[prism_b_z2gp[i],prism_b_z2gp[i]], linewidth=2,color='k') #bottom
            f3_ax2.plot([bot_pri_x2[i],bot_pri_x2[i]],[prism_b_z2gp[i],prism_b_z1gp[i]], linewidth=2,color='k') #left
            
#Plot the gravity stations
f3_ax2.scatter(bot_mp, mp_elev, color='#0081f1')
#Plot the drillings
f3_ax2.scatter(bot_drill,drillelev,marker='d', s=50, color='#db1e2a')
#Plot the gravity profile at the bottom
f3_ax2.plot()
f3_ax2.plot([bot_dist[0], bot_dist[-1]],[ylimits_plot_bot[0]+10, ylimits_plot_bot[0]+10],'-', linewidth=5, color='#fbd409', zorder=0)
f3_ax2.scatter(bot_dist[0],ylimits_plot_bot[0]+10, s=300, marker='*', color='#fbd40c', edgecolors='k')
f3_ax2.scatter(bot_dist[-1],ylimits_plot_bot[0]+10, s=150, marker='o', color='#fbd40c', edgecolors='k')

#Write the model name
f3_ax2.text(0,name_height, '%s' %str(model_name), color= '#df7116ff', horizontalalignment='center',verticalalignment='top', size= 30, fontweight='bold', **hfont, zorder=1)

#Write the axes labels
f3_ax2.set_ylabel("Elevation (m)", size=20)
f3_ax2.set_xlabel("Distance from the maximum anomaly, along the gravity profile (m)", size=20)


f3_ax1.set_ylim(ylimits_plot_top)
f3_ax2.set_ylim(ylimits_plot_bot)
f3_ax1.set_xlim(xlimits_plot)
f3_ax2.set_xlim(xlimits_plot_bot)
#f3_ax2.invert_xaxis()
f3_ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

#Set figure as fullscreen
fig3.set_size_inches(25.7, 13.1)

#make plot full screen
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
#Show figure
plt.show()
#plt.close()

#Save the plot as a .png and .svg figures
plt.savefig("Data_files/"+str(model_name)+'.svg', bbox_inches='tight')
plt.savefig("Data_files/"+str(model_name)+'.png', bbox_inches='tight')


##################### 2D map view plot ################

plt.rc('font', weight='bold') 

#Gravity profile
gp_chstart_2D=gp_chstart[:2]*convertkm
gp_chend_2D=gp_chend[:2]*convertkm

gp_b_ch_2D=np.asarray(gp_b_ch)*convertkm

#Top prism
pri_ch_E_2D=pri_ch_E*convertkm
pri_ch_N_2D=pri_ch_N*convertkm

prism_b_ch_2D=all_pri_corners_ch[:,:,:2]*convertkm

tpc_b_x1ch_2D=tpc_b_x1ch*convertkm
tpc_b_x2ch_2D=tpc_b_x2ch*convertkm
tpc_b_x3ch_2D=tpc_b_x3ch*convertkm
tpc_b_x4ch_2D=tpc_b_x4ch*convertkm
tpc_b_y1ch_2D=tpc_b_y1ch*convertkm
tpc_b_y2ch_2D=tpc_b_y2ch*convertkm
tpc_b_y3ch_2D=tpc_b_y3ch*convertkm
tpc_b_y4ch_2D=tpc_b_y4ch*convertkm

#Gravity stations
mp_xch_2D=mp_xch*convertkm
mp_ych_2D=mp_ych*convertkm

mp_b_xch_2D=mp_b_xch*convertkm
mp_b_ych_2D=mp_b_ych*convertkm

#Drillings
drillxch_2D=drillxch*convertkm
drillych_2D=drillych*convertkm

drillings_b_xch_2D=drillings_b_xch*convertkm
drillings_b_ych_2D=drillings_b_ych*convertkm

#Plotting the prism
vertices=np.zeros((nb_prisms, 5, 2))

for i in range(nb_prisms):
    vertices[i]=prism_verts(prism_b_ch_2D[i,:,0],prism_b_ch_2D[i,:,1])

verts_sqa_1 = [
    (pri_ch_E_2D[2], pri_ch_N_2D[2]), # left, bottom
    (pri_ch_E_2D[3], pri_ch_N_2D[3]), # left, top
    (pri_ch_E_2D[1], pri_ch_N_2D[1]), # right, top
    (pri_ch_E_2D[0], pri_ch_N_2D[0]), # right, bottom
    (pri_ch_E_2D[2], pri_ch_N_2D[2]), # ignored
    ]

#Initial top prism 
verts_sqa_2 = [
    (tpc_b_x3ch_2D, tpc_b_y3ch_2D), # left, bottom
    (tpc_b_x4ch_2D, tpc_b_y4ch_2D), # left, top
    (tpc_b_x2ch_2D, tpc_b_y2ch_2D), # right, top
    (tpc_b_x1ch_2D, tpc_b_y1ch_2D), # right, bottom
    (tpc_b_x3ch_2D, tpc_b_y3ch_2D), # ignored
    ]

codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]



    
path_1 = Path(verts_sqa_1, codes)
path_2 = Path(verts_sqa_2, codes)

#Intersection 

prism_int_b_x1ch_2D=prism_int_b_x1ch*convertkm
prism_int_b_x2ch_2D=prism_int_b_x2ch*convertkm
prism_int_b_y1ch_2D=prism_int_b_y1ch*convertkm
prism_int_b_y2ch_2D=prism_int_b_y2ch*convertkm

fig1,ax=plt.subplots()

ax.scatter(mp_xch_2D, mp_ych_2D, c='#0080ff',s=150, alpha=1, label='MP original')
ax.scatter(mp_b_xch_2D, mp_b_ych_2D, c='#df7116ff',s=50, alpha=1, label='Back rotated results')

ax.scatter(LCS_xori*convertkm,LCS_yori*convertkm,marker='P', s=200, zorder=1, linewidth=2, edgecolor="k" ,c="w" ,label='LCS origin')

ax.scatter(drillxch_2D, drillych_2D,marker='d', s=100, color='#db1e2a', label='Drillings original')
ax.scatter(drillings_b_xch_2D, drillings_b_ych_2D,marker='d', s=50, color='k', label='Back rotated results')

ax.plot([gp_chstart_2D[0], gp_chend_2D[0]],[gp_chstart_2D[1], gp_chend_2D[1]],'-', linewidth=5, color='#fbd409', zorder=0, label='Profile original')
ax.plot(gp_b_ch_2D[0],gp_b_ch_2D[1],'-', linewidth=2, color='c', zorder=0, label='Back rotated results')
ax.scatter(gp_chstart[0]*convertkm,gp_chstart[1]*convertkm, s=300, marker='*', color='#fbd40c', edgecolors='k', label='Profile start')
ax.scatter(gp_chend[0]*convertkm,gp_chend[1]*convertkm, s=150, marker='o', color='#fbd40c', edgecolors='k', label='Profile end')

ax.scatter(prism_int_b_x1ch_2D,prism_int_b_y1ch_2D, color='k')
ax.scatter(prism_int_b_x2ch_2D,prism_int_b_y2ch_2D, color='k')

ax.set_xlabel('Distance along the $x_{CH}$ axis (km)',labelpad=15,fontsize=18)
ax.set_ylabel('Distance along the $y_{CH}$ axis (km)',labelpad=15,fontsize=18)


patch = patches.PathPatch(path_1,edgecolor='#bcbcbc', facecolor='#ededed', lw=5, zorder=0, label='Prism original')
patch2 = patches.PathPatch(path_2, edgecolor='k',facecolor='none', lw=1, zorder=0, label='Back rotated results')

#Create a path for the prisms to be plotted
paths=np.zeros((nb_prisms,2,5))

#Generate a list of RGB colors randomly depending on the number of prisms
color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(nb_prisms)]

#Plot top prism and initial prism
ax.add_patch(patch)
ax.add_patch(patch2)

#Plot all the prisms
for i in range(1,nb_prisms):
    paths= Path(vertices[i], codes)
    patch_i = patches.PathPatch(paths, edgecolor=color[i], facecolor='none', lw=1, label='Prism '+str(i+1))
    ax.add_patch(patch_i)
    


ax.set_xlim(xlimits_2D[0],xlimits_2D[1])
ax.set_ylim(ylimits_2D[0],ylimits_2D[1])
ax.set_aspect('equal')
ax.legend(title=str(model_name), loc='lower right')

plt.show()

####################### Same plot as above but with stations numbers ###########################

fig2,ax=plt.subplots()

ax.scatter(mp_xch_2D, mp_ych_2D, c='#0080ff',s=150, alpha=1, label='MP original')
ax.scatter(mp_b_xch_2D, mp_b_ych_2D, c='#df7116ff',s=50, alpha=1, label='Back rotated results')

for i in range(nb_stations):
    ax.annotate(str('%.f' % MP_num[i]),xy=[mp_xch_2D[i], mp_ych_2D[i]], size='10')

ax.scatter(LCS_xori*convertkm,LCS_yori*convertkm,marker='P', s=200, zorder=1, linewidth=2, edgecolor="k" ,c="w" ,label='LCS origin')

ax.scatter(drillxch_2D, drillych_2D,marker='d', s=100, color='#db1e2a', label='Drillings original')
ax.scatter(drillings_b_xch_2D, drillings_b_ych_2D,marker='d', s=50, color='k', label='Back rotated results')

ax.plot([gp_chstart_2D[0], gp_chend_2D[0]],[gp_chstart_2D[1], gp_chend_2D[1]],'-', linewidth=5, color='#fbd409', zorder=0, label='Profile original')
ax.plot(gp_b_ch_2D[0],gp_b_ch_2D[1],'-', linewidth=2, color='c', zorder=0, label='Back rotated results')
ax.scatter(gp_chstart[0]*convertkm,gp_chstart[1]*convertkm, s=300, marker='*', color='#fbd40c', edgecolors='k', label='Profile start')
ax.scatter(gp_chend[0]*convertkm,gp_chend[1]*convertkm, s=150, marker='o', color='#fbd40c', edgecolors='k', label='Profile end')

ax.scatter(prism_int_b_x1ch_2D,prism_int_b_y1ch_2D, color='k')
ax.scatter(prism_int_b_x2ch_2D,prism_int_b_y2ch_2D, color='k')

ax.set_xlabel('Distance along the $x_{CH}$ axis (km)',labelpad=15,fontsize=18)
ax.set_ylabel('Distance along the $y_{CH}$ axis (km)',labelpad=15,fontsize=18)


patch = patches.PathPatch(path_1,edgecolor='#bcbcbc', facecolor='#ededed', lw=5, zorder=0, label='Prism original')
patch2 = patches.PathPatch(path_2, edgecolor='k',facecolor='none', lw=1, zorder=0, label='Back rotated results')

#Create a path for the prisms to be plotted
paths=np.zeros((nb_prisms,2,5))

#Plot top prism
ax.add_patch(patch)
ax.add_patch(patch2)

#Plot all the prisms
for i in range(1,nb_prisms):
    paths= Path(vertices[i], codes)
    patch_i = patches.PathPatch(paths, edgecolor=color[i], facecolor='none', lw=1, label='Prism '+str(i+1))
    ax.add_patch(patch_i)
    



ax.set_xlim(xlimits_2D[0],xlimits_2D[1])
ax.set_ylim(ylimits_2D[0],ylimits_2D[1])
ax.set_aspect('equal')
ax.legend(title=str(model_name), loc='lower right')

plt.show()


##################################### 2D map view of the LCS ########################################

convertkm=1/1000

prism_x1_2D=prism_x1*convertkm
prism_x2_2D=prism_x2*convertkm
prism_y1_inter_2D=prism_y1_inter*convertkm
prism_y2_inter_2D=prism_y2_inter*convertkm

mp_xlcs_2D=mp_xlcs*convertkm
mp_ylcs_2D=mp_ylcs*convertkm

drillings_xlcs_2D=drillings_xlcs*convertkm
drillings_ylcs_2D=drillings_ylcs*convertkm

gp_lcs=np.asarray(gp_lcs)

gp_lcs_2D=gp_lcs*convertkm

#Prepraring the prisms LCS corners data
all_pri_corners_lcs=np.zeros((nb_prisms,4, 2))

for i in range(nb_prisms):
    all_pri_corners_lcs[i,0][:2]=botleftcorners[i]
    all_pri_corners_lcs[i,1][:2]=topleftcorners[i]
    all_pri_corners_lcs[i,2][:2]=toprightcorners[i]
    all_pri_corners_lcs[i,3][:2]=botrightcorners[i]

all_pri_corners_lcs_2D=all_pri_corners_lcs*convertkm
vertices_LCS=np.zeros((nb_prisms, 5, 2))

for i in range(nb_prisms):
    vertices_LCS[i]=prism_verts(all_pri_corners_lcs_2D[i,:,0],all_pri_corners_lcs_2D[i,:,1])

codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]

#Setting up the 2D map view LCS figure
fig3,ax=plt.subplots()

ax.scatter(mp_xlcs_2D, mp_ylcs_2D, c='#0080ff',s=100, alpha=1, zorder=0)

ax.scatter(drillings_xlcs_2D, drillings_ylcs_2D,marker='d', s=50, color='#db1e2a', zorder=0)

#Plot the LCS origin 
ax.scatter(0,0,marker='P', s=200, zorder=1, linewidth=2, edgecolor="k" ,c="w", label='LCS origin')

ax.plot(gp_lcs_2D[0,:],gp_lcs_2D[1,:],'o-', linewidth=5, color='#fbd409', zorder=0)
ax.scatter(gp_lcs_2D[0,0],gp_lcs_2D[1,0], s=300, marker='*', color='#fbd40c', edgecolors='k', label='Profile start')
ax.scatter(gp_lcs_2D[0,1],gp_lcs_2D[1,1], s=150, marker='o', color='#fbd40c', edgecolors='k', label='Profile end')

#Plotting the intersection prisms and gravity profile
ax.scatter(prism_x1_2D,prism_y1_inter_2D, color='k')
ax.scatter(prism_x2_2D,prism_y2_inter_2D, color='k')

#Plotting the prisms
#Top prism
path_top_pri=Path(vertices_LCS[0], codes)
patch_top_pri=patches.PathPatch(path_top_pri, edgecolor='k', facecolor='none', lw=1, label='Top prism ')
ax.add_patch(patch_top_pri)
#The remaining prisms
for i in range(1,nb_prisms):
    paths= Path(vertices_LCS[i], codes)
    patch_i = patches.PathPatch(paths, edgecolor=color[i], facecolor='none', lw=1, label='Prism '+str(i+1))
    ax.add_patch(patch_i)

ax.set_xlabel('Distance along the $x_{LCS}$ axis (km)',labelpad=15,fontsize=18)
ax.set_ylabel('Distance along the $y_{LCS}$ axis (km)',labelpad=15,fontsize=18)
ax.set_xlim(xlimits_lcs_2D[0],xlimits_lcs_2D[1])
ax.set_ylim(ylimits_lcs_2D[0],ylimits_lcs_2D[1])
ax.legend(title=str(model_name), loc='lower right')
ax.set_aspect('equal')
ax.invert_xaxis()

plt.show()

#################### PLOT PERPENDICULAR TO THE VALLEY DIRECTION ####################
################ LCS 2D MAP VIEW ################
#Check plot switch
if perpendicular_plot=="yes" or perpendicular_plot=="yes" or perpendicular_plot=="yes" :
    mp_x_perp_test=mp_xlcs*convertkm
    mp_y_perp_test=mp_y_perp*convertkm
    mp_y_perp_2D=mp_ylcs_perp*convertkm
    drill_x_perp_2D=drill_x_perp*convertkm
    drill_y_perp_2D=drill_y_perp*convertkm
    DEM_x_perp_2D=DEM_x_perp*convertkm
    DEM_y_perp_2D=DEM_y_perp*convertkm
    bedrock_x_perp_2D=bedrock_x_perp*convertkm
    bedrock_y_perp_2D=bedrock_y_perp*convertkm
    
    #Setting up the 2D map view LCS figure
    fig4,ax=plt.subplots()
    
    #MPs with the LCS coordinates
    ax.scatter(mp_xlcs_2D, mp_ylcs_2D, c='#0080ff',s=200, alpha=1, zorder=0)
    #MPs projected onto the gravity profile, parallel to the LCS y-axis...
    ax.scatter(mp_xlcs_2D, mp_y_perp_2D, c='#df7116ff',s=100, alpha=1, zorder=0)
    #... then projected onto the LCS x axis. These are used for the plot parallel to the valley direction
    ax.scatter(mp_x_perp_test, mp_y_perp_test, c='g',s=100, alpha=1, zorder=0)
    
    ax.scatter(drillings_xlcs_2D, drillings_ylcs_2D,marker='d', s=100, color='#db1e2a', zorder=0)
    ax.scatter(drillings_xlcs_2D, drill_y_perp_2D,marker='d', s=50, color='k', zorder=0)
    
    #Plot the LCS origin 
    ax.scatter(0,0,marker='P', s=200, zorder=1, linewidth=2, edgecolor="k" ,c="w", label='LCS origin')

    ax.plot(gp_lcs_2D[0,:],gp_lcs_2D[1,:],'o-', linewidth=5, color='#fbd409', zorder=0)
    ax.scatter(gp_lcs_2D[0,0],gp_lcs_2D[1,0], s=300, marker='*', color='#fbd40c', edgecolors='k', label='Profile start')
    ax.scatter(gp_lcs_2D[0,1],gp_lcs_2D[1,1], s=150, marker='o', color='#fbd40c', edgecolors='k', label='Profile end')

    #Plotting the intersection prisms and gravity profile
    ax.scatter(prism_x1_2D,prism_y1_inter_2D, color='k')
    ax.scatter(prism_x2_2D,prism_y2_inter_2D, color='k')

#Plotting the prisms
#Top prism
    path_top_pri=Path(vertices_LCS[0], codes)
    patch_top_pri=patches.PathPatch(path_top_pri, edgecolor='k', facecolor='none', lw=1, label='Top prism ')
    ax.add_patch(patch_top_pri)
#The remaining prisms
    for i in range(1,nb_prisms):
        paths= Path(vertices_LCS[i], codes)
        patch_i = patches.PathPatch(paths, edgecolor=color[i], facecolor='none', lw=1, label='Prism '+str(i+1))
        ax.add_patch(patch_i)

#plot the lcs axes
    ax.plot([0,0],[-2,2], color='#db1e2a', linewidth=4, zorder=0, label='LCS axes')
    ax.plot([-2,2],[0,0], color='#db1e2a', linewidth=4, zorder=0)


    ax.set_xlabel('Distance along the $x_{LCS}$ axis (km)',labelpad=15,fontsize=18)
    ax.set_ylabel('Distance along the $y_{LCS}$ axis (km)',labelpad=15,fontsize=18)
    ax.set_xlim(xlimits_lcs_2D[0],xlimits_lcs_2D[1])
    ax.set_ylim(ylimits_lcs_2D[0],ylimits_lcs_2D[1])
    ax.legend(title=str(model_name), loc='lower right')
    ax.set_aspect('equal')
    ax.invert_xaxis()

    plt.show()
    
################## Residual anomaly plot for parallel ylcs ##################
    fig5 = plt.figure(constrained_layout=True)
    widths = [3, 2]
    heights = [3, 1.5]
    gs = fig5.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
                          height_ratios=heights)
    ax1 = fig5.add_subplot(gs[0, :])
    ax2 = fig5.add_subplot(gs[1, :])
################## Top plot ##################  
#Observed residual anomaly
    ax1.scatter(mp_x_perp,residual, color='#0081f1', s=300,zorder=1, label="Observed residual anomaly")
    ax1.errorbar(mp_x_perp,residual, yerr=0.13, color='k', capsize=5,linestyle='None',zorder=0)

#Modelled residual anomaly
    ax1.scatter(mp_x_perp,model_results, color='#df7116', s=300,zorder=1, label="Modelled residual anomaly")


    ax1.plot(xlimits_plot, [0,0], color='k', linestyle="--", linewidth=2, dashes=[10,3,10, 3])
    
#Text for orientation
    ax1.text(0,left_side, '%s' %str(left_label),horizontalalignment='center',verticalalignment='top', size= 30, fontweight='bold',zorder=1)
    ax1.text(profile_length,right_side, '%s' %str(right_label),horizontalalignment='center',verticalalignment='top', size= 30, fontweight='bold',zorder=1)

#Write the different density contrast values used for this model


    if len(rho)==1:
        ax1.text(xrho,yrho, '%s' %'Δρ(kg/m3)='+str("{0:<10.0f}".format(rho[0])),horizontalalignment='left',verticalalignment='top',fontweight='bold', **hfont, size= 24,zorder=1)
    else:
        ax1.text(xrho,yrho, '%s' %'Δρ(kg/m3)= '+','.join(list(map("{0:^10.0f}".format,rho))),horizontalalignment='left',verticalalignment='top',fontweight='bold', **hfont, size= 24,zorder=1)
        

#Write the axes titles
    ax1.set_xlabel("Distance along the gravity profile (m)", size=20)
    ax1.set_ylabel("Residual Anomaly (mGal)", size=20)
    ax1.legend(loc='lower left')

################## Bottom plot ##################

#Subtract coordinate for to have the distance from the maximum anomaly stations
#DEM and bedrock profiles
    bot_dist=DEM_x_perp-pos_max_ano
    bot_dist_b=bedrock_x_perp-pos_max_ano
#Prism coordinates
    if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
        bot_pri_x2_cross=prism_x2_cross
        bot_pri_x1_cross=prism_x1_cross

    else:
        bot_pri_x2=prism_x2
        bot_pri_x1=prism_x1
        
#Gravity stations
    bot_mp=mp_x_perp-pos_max_ano
#Drillings
    bot_drill=drill_x_perp-pos_max_ano
#x axis limits for the bottom plot
    xlimits_plot_bot=xlimits_plot-pos_max_ano

#DEM and bedrock plots
    ax2.plot(bot_dist,DEM_elev, color='#004b88',label='Elevation')
    ax2.plot(bot_dist_b,bedrock_elev, '--',color='#800000',label='Bedrock model')

#Plot the prisms with an exception for a single prism plot
    if all_prisms!='yes' and all_prisms!='Yes' and all_prisms!='YES':
        if nb_prisms_cross==1:
            ax2.plot([bot_pri_x2_cross,bot_pri_x1_cross],[prism_b_z1gp_cross,prism_b_z1gp_cross], linewidth=2,color='k') #top
            ax2.plot([bot_pri_x1_cross,bot_pri_x1_cross],[prism_b_z1gp_cross,prism_b_z2gp_cross], linewidth=2,color='k') #right
            ax2.plot([bot_pri_x2_cross,bot_pri_x1_cross],[prism_b_z2gp_cross,prism_b_z2gp_cross], linewidth=2,color='k') #bottom
            ax2.plot([bot_pri_x2_cross,bot_pri_x2_cross],[prism_b_z2gp_cross,prism_b_z1gp_cross], linewidth=2,color='k') #left
        else:
            for i in range(nb_prisms_cross):
                ax2.plot([bot_pri_x2_cross[i],bot_pri_x1_cross[i]],[prism_b_z1gp_cross[i],prism_b_z1gp_cross[i]], linewidth=2,color='k') #top
                ax2.plot([bot_pri_x1_cross[i],bot_pri_x1_cross[i]],[prism_b_z1gp_cross[i],prism_b_z2gp_cross[i]], linewidth=2,color='k') #right
                ax2.plot([bot_pri_x2_cross[i],bot_pri_x1_cross[i]],[prism_b_z2gp_cross[i],prism_b_z2gp_cross[i]], linewidth=2,color='k') #bottom
                ax2.plot([bot_pri_x2_cross[i],bot_pri_x2_cross[i]],[prism_b_z2gp_cross[i],prism_b_z1gp_cross[i]], linewidth=2,color='k') #left
    else:        
        if nb_prisms==1:
            ax2.plot([bot_pri_x2,bot_pri_x1],[prism_b_z1gp,prism_b_z1gp], linewidth=2,color='k') #top
            ax2.plot([bot_pri_x1,bot_pri_x1],[prism_b_z1gp,prism_b_z2gp], linewidth=2,color='k') #right
            ax2.plot([bot_pri_x2,bot_pri_x1],[prism_b_z2gp,prism_b_z2gp], linewidth=2,color='k') #bottom
            ax2.plot([bot_pri_x2,bot_pri_x2],[prism_b_z2gp,prism_b_z1gp], linewidth=2,color='k') #left
        else:
            for i in range(nb_prisms):
                ax2.plot([bot_pri_x2[i],bot_pri_x1[i]],[prism_b_z1gp[i],prism_b_z1gp[i]], linewidth=2,color='k') #top
                ax2.plot([bot_pri_x1[i],bot_pri_x1[i]],[prism_b_z1gp[i],prism_b_z2gp[i]], linewidth=2,color='k') #right
                ax2.plot([bot_pri_x2[i],bot_pri_x1[i]],[prism_b_z2gp[i],prism_b_z2gp[i]], linewidth=2,color='k') #bottom
                ax2.plot([bot_pri_x2[i],bot_pri_x2[i]],[prism_b_z2gp[i],prism_b_z1gp[i]], linewidth=2,color='k') #left
    
#Plot the gravity stations
    ax2.scatter(bot_mp, mp_elev, color='#0081f1')
#Plot the drillings
    ax2.scatter(bot_drill,drillelev,marker='d', s=50, color='#db1e2a')
#Plot the gravity profile at the bottom
    ax2.plot()
    ax2.plot([bot_dist[0], bot_dist[-1]],[ylimits_plot_bot[0]+10, ylimits_plot_bot[0]+10],'-', linewidth=5, color='#fbd409', zorder=0)
    ax2.scatter(bot_dist[0],ylimits_plot_bot[0]+10, s=300, marker='*', color='#fbd40c', edgecolors='k')
    ax2.scatter(bot_dist[-1],ylimits_plot_bot[0]+10, s=150, marker='o', color='#fbd40c', edgecolors='k')

#Write the model name
    ax2.text(0,name_height, '%s' %str(model_name), color= '#df7116ff', horizontalalignment='center',verticalalignment='top', size= 30, fontweight='bold', **hfont, zorder=1)
#Label for parallel to valley
    ax2.text(xlimits_plot_bot[0]+100,300, 'Perpendicular to the valley direction', color= 'k', horizontalalignment='left',verticalalignment='top', size= 24, fontweight='bold', **hfont, zorder=1)
#Write the axes labels
    ax2.set_ylabel("Elevation (m)", size=20)
    ax2.set_xlabel("Distance from the maximum anomaly, along the gravity profile (m)", size=20)

    ax1.set_ylim(ylimits_plot_top)
    ax2.set_ylim(ylimits_plot_bot)
    ax1.set_xlim(xlimits_plot)
    ax2.set_xlim(xlimits_plot_bot)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

#Set figure as fullscreen
    fig5.set_size_inches(25.7, 13.1)

#make plot full screen
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
#Show figure
    plt.show()

#Save the plot as a .png and .svg figures
    plt.savefig("Data_files/"+str(model_name)+'_perpendicular.svg', bbox_inches='tight')
    plt.savefig("Data_files/"+str(model_name)+'_perpendicular.png', bbox_inches='tight')


#### NEED TO ADD CASE FOR SINGLE MPS OR PRISM ######
################### WRITING RESULTS FILE ################### 

results_file_name=str(model_name)+'_PRISM_results.txt'

file_results=open('Display/'+str(results_file_name),'w')
np.set_printoptions(formatter={'float': '{: .4f}'.format})
run_day_time=datetime.datetime.now()

print("#Display PRISMA results run date and start time:", file=file_results)
print("#"+str("%s"%run_day_time), file=file_results)
print(" ", file=file_results)

print("#PRISMA model data", file=file_results)
print(str(PRISMA_output_file), file=file_results)
print(" ", file=file_results)

print("#DEM and bedrock profile files", file=file_results)
print(str(DEM_file), file=file_results)
print(str(bedrock_file), file=file_results)
print(" ", file=file_results)

with open(PRISMA_output_file, 'r') as file_input:
    #Make a list with each of the lines of the input file
    lines = file_input.read().splitlines()
    #Write the header lines, except the gravity stations header, in the new output file
    for i in range(length_header-1):
        print(str(lines[i]), file=file_results)
        
#Gravity stations  
if nb_prisms==1:
    print("#MP Station \t| CH coordinates(m) E,\tN,\tElev\t| Distance(m) along\taway from the profile\t| Distance check \t|" +
          " LCS coordinates(m) x,\ty,\tz\t| Bouguer anomaly(mGal)\t| Regional anomaly(mGal)\t| Residual anomaly(mGal)\t|"+
          " Prism gravity effect (mGal)\t| Back rotated CH-coordinates E,\tN\t| Back rotated projected coordinates x,\ty", file=file_results)
else:
    print("#MP Station \t| CH coordinates(m) E,\tN,\tElev\t| Distance(m) along\taway from the profile\t| Distance check \t|" +
          " LCS coordinates(m) x,\ty,\tz\t| Bouguer anomaly(mGal)\t| Regional anomaly(mGal)\t| Residual anomaly(mGal)\t|"+
          " All prisms gravity effect (mGal)\t| Invididual prisms gravity effect (mGal)\t| Back rotated CH-coordinates E,\tN\t| Back rotated projected coordinates x,\ty", file=file_results)
for i in range(nb_stations):
    #Write existing data
    print(str("{0:<10.0f}".format(MP_num[i])),"\t","\t".join(str("{0:<10.3f}".format(mpchcoor)) for mpchcoor in MPs_all_data[1:6,i]),"\t", MP_distance_check[i],"\t",
          "\t".join(str("{0:<10.3f}".format(mpchcoor)) for mpchcoor in MPs_all_data[7:,i]),"\t",
          #Back rotated LCS coordinates to initial coordinate system
          str("{0:<10.3f}".format(mp_b_xch[i])),"\t",str("{0:<10.3f}".format(mp_b_ych[i])),"\t",str("{0:<10.3f}".format(mp_b_zch[i])),"\t",
          #Projected back rotated coordinates on the gravity profile
          str("{0:<10.3f}".format(mp_b_gp[0][i])),"\t",str("{0:<10.3f}".format(mp_b_gp[1][i])),"\t", file=file_results)
print(" ", file=file_results)

#Drillings
print("#Drilling number\t| CH coordinates E,\tN,\tElev\t| Distance along \t away from the profile\t| LCS coordinates x,\ty,\tz\t| Back rotated CH-coordinates E,\tN\t| Back rotated projected coordinates x,\ty", file=file_results)
for i in range(nb_drillings):
    #Write existing data
    print(str("{0:<10.0f}".format(drill_num[i])),"\t","\t".join(str("{0:<10.3f}".format(dril_all)) for dril_all in drillings_data[1:6,i]),"\t",drillings_distance_check[i],"\t",
          "\t".join(str("{0:<10.3f}".format(dril_all)) for dril_all in drillings_data[7:,i]),"\t",
          #Back rotated LCS coordinates to initial coordinate system
          str("{0:<10.3f}".format(drillings_b_xch[i])),"\t",str("{0:<10.3f}".format(drillings_b_ych[i])),"\t",str("{0:<10.3f}".format(drillings_b_zch[i])),"\t",
          #Projected back rotated coordinates on the gravity profile
          str("{0:<10.3f}".format(drillings_b_gp[0][i])),"\t",str("{0:<10.3f}".format(drillings_b_gp[1][i])),file=file_results)

with open(PRISMA_output_file, 'r') as file_input:
    #Make a list with each of the lines of the input file
    lines = file_input.read().splitlines()
    #Write the header lines, except the gravity stations header, in the new output file
    for i in range(length_header+nb_stations+2+nb_drillings, len(lines)):
        print(str(lines[i]), file=file_results)
print(" ", file=file_results)

#Intersection coordinates, back rotated CH coordinates and LCS
print("#Intersection prism/gravity profile coordinates\t" ,file=file_results)
print("#Prism number\t| Intersection CH coordinates E,\tN \t| LCS coordinates x, \ty \t| ztop,\tzmean,\tzbot\t| Distance along gravity profile",file=file_results)
for i in range(nb_prisms):
        #This part writes the "maximum" coordinates xmax, ymax
        print(str("{0:<10.0f}".format(i+1)),"\t",
              #LCS coordinates
              str("{0:<10.2f}".format(prism_x2[i])),"\t",str("{0:<10.2f}".format(prism_y2_inter[i])),"\t",
              #Back rotated CH coordinates
              str("{0:<10.2f}".format(prism_int_b_x2ch[i])),"\t",str("{0:<10.2f}".format(prism_int_b_y2ch[i])),"\t",
              #elevation highest value, mean of the prism and lowest value
              str("{0:<10.2f}".format(prism_b_z1gp[i])),"\t",
              str("{0:<10.2f}".format((prism_b_z1gp[i]+prism_b_z2gp[i])/2)),"\t",str("{0:<10.2f}".format(prism_b_z2gp[i])),"\t",
              #Coordinates projected along the profile
              str("{0:<10.2f}".format(prism_int_b_x2gp[i])),"\t",str("{0:<10.2f}".format(prism_int_b_y2gp[i])),file=file_results)

        #This part writes the "minimum" coordinates xmin, ymin
        print(str("{0:<10.0f}".format(i+1)),"\t",
              #LCS coordinates
              str("{0:<10.2f}".format(prism_x1[i])),"\t",str("{0:<10.2f}".format(prism_y1_inter[i])),"\t",
              #Back rotated CH coordinates
              str("{0:<10.2f}".format(prism_int_b_x1ch[i])),"\t",str("{0:<10.2f}".format(prism_int_b_y1ch[i])),"\t",
              #elevation highest value, mean of the prism and lowest value
              str("{0:<10.2f}".format(prism_b_z1gp[i])),"\t",
              str("{0:<10.2f}".format((prism_b_z1gp[i]+prism_b_z2gp[i])/2)),"\t",str("{0:<10.2f}".format(prism_b_z2gp[i])),"\t",
              #Coordinates projected along the profile
              str("{0:<10.2f}".format(prism_int_b_x1gp[i])),"\t",str("{0:<10.2f}".format(prism_int_b_y1gp[i])),file=file_results)
        
#Write the coordinates of the maximum anomaly coordinate with the maximum depth value
print(str("{0:<10.0f}".format(i+1)),"\t",str("{0:<10.2f}".format(LCS_xori)),"\t",str("{0:<10.2f}".format(LCS_yori)),"\t",
      str("{0:<10.2f}".format(prism_b_z2gp[nb_prisms-1])),"\t#Maximum anomaly coordinates and max depth",file=file_results)

#Close the output file after all the results are written
file_results.close()

#Make a copy of the output file in the Data_files folder
shutil.copy2('Display/'+str(results_file_name), "Data_files")
