# -*- coding: utf-8 -*-
"""
Created on 10/03/2021

@author: Dimitri
"""

import numpy as np
import math
import itertools
import sys
import os
import time
import datetime
import shutil
import fileinput


###########################
# LOADING USER INPUTS FILES
###########################
#Number of spaces and headers in the input file. This input is needed to skip the appropriate number of lines when reading the input
#DO NOT CHANGE THIS VALUE UNLESS THE INPUT FORMAT IS MODIFIED TOO
#skipheader skips all the lines comments AND spaces while skipfooter skip only takes into account non commented lines
#The number of spaces and header is fixed while the number of data varies. The added numbers in skipheader corresponds to the comments and spaces to skip on top of the number of points

#Fixed length of the header before the MPs input
length_header=36

#Directly read the number of stations, drillings and prisms to properly extract the MPs, and prisms data arrays
with open(PRISMA_input_file) as file_input:
   nb_file_input=file_input.readlines()
#Load the number of stations from the file, this position is the line number in the file and is fixed
#Used to extract the input data and to check the user's input
   number_stations=int(nb_file_input[30])
   
#Load the number of drillings from the file, this position is the line number in the file and is fixed
#Only needed to extract the input data. Drillings values are not used in this code
   number_drillings=int(nb_file_input[33])
   
#Load the number of prisms from the file, this position is the line number in the file and is fixed
#Used to extract the input data and to check the user's input
   prisms_number=int(nb_file_input[length_header + number_stations + 2 + number_drillings + 3])

#This check is complementary to the check done below which directly compares the number of stations loaded and the number of stations given as input
### ADD CHECK USER INPUT AND FILE INPUT NUMBER HERE ###

#Getting the MPs data
#Skipheader is the fixed number of lines at the start of the file
MPs_data_header=length_header
#Skipfooter is the number of drillings, the number of prisms, the number of intersections prism/gravity profiles and the max depth coordinates
MPs_data_footer=number_drillings + 4 + prisms_number 
#Load the stations coordinates
MPs_all_data =  np.genfromtxt(PRISMA_input_file, skip_header=MPs_data_header, skip_footer=MPs_data_footer)
#Load the distance check results separately, as they are string data and are loaded above as "NaN"
MP_distance_check = np.genfromtxt(PRISMA_input_file, delimiter='\t', skip_header=MPs_data_header, skip_footer=MPs_data_footer, usecols=6, dtype=str)

if routine=='PRISMA':
#Getting the prism data
#Skipheader is the fixed number of lines at the start of the file and the MPs,the drillings and their header, the prisms coordinates header and the prism corner coordinates
    prism_data_header= length_header + number_stations + 2 + number_drillings + 13
#Skipfooter is the number of intersections prism/gravity profiles and the max depth coordinates
    prism_data_footer= 0
#Load the prisms coordinates and density values for PRISMA
    pri_data = x2s, x1s, y2s, y1s, z2s, z1s, rho_pri =np.genfromtxt(PRISMA_input_file, skip_header=prism_data_header, skip_footer=prism_data_footer, unpack=True)

#Load the polygons/isolines/height lines coordinates for BGPOLY
if routine=='BGPOLY':
    poly_data=x_poly, y_poly, z_poly=np.loadtxt(BGPOLY_input_file, unpack=True)
#Load the polygons/isolines/height lines density values for BGPOLY
    rho_poly=np.loadtxt(density_poly_file_name)


#FLAG FOR DEBUGGING MODE 0=OFF, 1=ON 
#CHANGE THIS FLAG ONLY IF YOU WANT ALL INTERMEDIATE CALCULATIONS RESULTS
#IT WILL CREATE INDIVIDUAL FILES FOR EACH STATIONS/MP 
DEBUG=0

#########################################
#FUNCTIONS DEFINITION - Part 1 BGPoly
#########################################
#These are the definitions of the functions used in the routines below

#Check if coordinates have n as their height and store them in isoline if it's the case
def poly_line(poly_list,n,isoline):
    #Its purpose is to create an array of (x,y,z) points for a specific height
    #For all the points available search if z is equal to n
    if poly_list[2]==n:
    #in which case store it in isoline
        isoline.append(poly_list)
    return isoline,poly_list

#To avoid dividing by 0 when calculating the body gravity effect     
def zero_line(x_iso,y_iso):               
    #check if x_iso = 0 and replace it by 0.1  
    if x_iso==0 :
        x_iso=0.1 
    #check if y_iso = 0 and replace it by 0.1        
    if y_iso[i]==0 :
        y_iso[i]=0.1        
    return x_iso, y_iso

#For each MPs, create a local coordinate system where the MP is considered the origin.
#New height lines coordinates (x,y,z) are centered around the MPs
def coor_transform(mp,isoline,iso_norm):
    x_mp=mp[0]
    y_mp=mp[1]
    z_mp=mp[2]
    iso_norm[0] = isoline[0]-x_mp
    iso_norm[1] = isoline[1]-y_mp
    iso_norm[2] = isoline[2]-z_mp
    return iso_norm    

#function to calculate pi,qi,fi, etc.... for each isoline according to Talwani and Ewing (1960)
#Avoid diving by 0
# and calculate i+1 for every value
# Or do a While loop and specify the last value when i=len(x_iso) so it takes value[i+1]=value[0] as it's a polygon 
def r_def(x_iso,y_iso):
    r=(x_iso*x_iso + y_iso*y_iso)**(1/2)
    return r

def sum_ri_def(x_iso1,x_iso2,y_iso1,y_iso2):
    sum_ri=((x_iso1-x_iso2)**2 + (y_iso1-y_iso2)**2)**(1/2)
    return sum_ri

def p_def(x_iso1,x_iso2,y_iso1,y_iso2,sum_ri):
    p=(y_iso1-y_iso2)/sum_ri*x_iso1 - (x_iso1-x_iso2)/sum_ri*y_iso1
    return p

def m_def(x_iso1,x_iso2,y_iso1,y_iso2,r,r_1):
    m=(y_iso1/r)*(x_iso2/r_1) - (y_iso2/r_1)*(x_iso1/r)
    return m

def q_def(x_iso1,x_iso2,y_iso1,y_iso2,sum_ri,r):
    q=(x_iso1-x_iso2)/sum_ri*x_iso1/r + (y_iso1-y_iso2)/sum_ri*y_iso1/r
    return q

def f_def(x_iso1,x_iso2,y_iso1,y_iso2,sum_ri,r_1):
    f=(x_iso1-x_iso2)/sum_ri*x_iso2/r_1 + (y_iso1-y_iso2)/sum_ri*y_iso2/r_1
    return f

#Function to give the sign of W and S depending on the value of m and p respectively
def S_sign(p):
    
    if p < 0:
        S=-1
    else :
        S=1
    return S

def W_sign(m):
    
    if m < 0:
        W=-1
    else :
        W=1
    return W

def bgpoly_arcsin(z_iso, q, p, f, S):
    #Avoiding nan for np.arcsin(x) because of float format where arcsin_exp = 1.000000002
    arcsin_exp1=z_iso*q*S/(p*p + z_iso*z_iso)**(1/2)
    arcsin_exp2=z_iso*f*S/(p*p + z_iso*z_iso)**(1/2)
    diff_arcsin= np.abs(arcsin_exp1-arcsin_exp2)
    
    if diff_arcsin <= 0.000001:
        G_E_arcsin= - np.arcsin(1) + np.arcsin(1)
        if DEBUG==1:
            print('replaced arcsin value', file=file)
            print('iso eff arcsin1 = '+str(arcsin_exp1),file=file)    
            print('iso eff arcsin2 = '+str(arcsin_exp2),file=file)
            print('diff arcsin '+str(G_E_arcsin), file= file)    
    else :
        G_E_arcsin= - np.arcsin(arcsin_exp1) + np.arcsin(arcsin_exp2)
        

    return G_E_arcsin

def bgpoly_arccos(x_iso1, x_iso2, y_iso1, y_iso2, r, r_1):
    #Avoiding nan for np.arccos(x) because of float format where arccos_exp = 1.000000002
    arccos_exp=(x_iso1/r)*(x_iso2/r_1) + (y_iso1/r)*(y_iso2/r_1)
    
    if 1 < arccos_exp <= 1.000001 :
        G_E_arccos = 1
        if DEBUG==1:
            print('replaced arcsin value', file=file)
            print('old arccos expression '+str(arccos_exp), file= file) 
            print('new arccos expression '+str(G_E_arccos), file= file)    
    elif -1.000001 < arccos_exp < -1:
        G_E_arccos = -1
        if DEBUG==1:
            print('replaced arcsin value', file=file)
            print('old arccos expression '+str(arccos_exp), file= file) 
            print('new arccos expression '+str(G_E_arccos), file= file)    
    else :
        G_E_arccos = arccos_exp

    return G_E_arccos

#Function to calculate the gravity effect of height lines
def iso_eff(z_iso, x_iso1, x_iso2, y_iso1, y_iso2, r, r_1, q, p, f, W, S):
    G_E_arcsin = bgpoly_arcsin(z_iso, q, p, f, S)
    G_E_arccos = bgpoly_arccos(x_iso1, x_iso2, y_iso1, y_iso2, r, r_1)
#Original expression
#W*np.arccos((x_iso1/r)*(x_iso2/r_1) + (y_iso1/r)*(y_iso2/r_1)) - np.arcsin(z_iso*q*S/(p*p + z_iso*z_iso)**(1/2)) + np.arcsin(z_iso*f*S/(p*p + z_iso*z_iso)**(1/2))
    G_E=W*np.arccos(G_E_arccos) + G_E_arcsin


    return G_E
#Function to do the interpolation between 3 isolines and return the gravity effect
def body_func(V1,V2,V3,z1,z2,z3):
    Body_eff=1/6*(V1*(z1-z3)/(z1-z2)*(3*z2-z3-2*z1) + V2*(z1-z3)**3/((z2-z3)*(z2-z1)) \
            + V3*(z1-z3)/(z3-z2)*(3*z2-z1-2*z3))
    return Body_eff

#########################################
#FUNCTIONS DEFINITION  part 2 PRISMA
#########################################

#Function to calculate the gravity effect of one prism, from Nagy, 1966 and Banerjee & Das Gulpa (1977).
#It takes the folloying entry parameters : x1,x2,y1,y2,z1,z2
#They correspond to the coordinates of the vertices of the prism.
#Depending on the signs of x1,x2,y1,y2 this function will use a different equation.

#These terms are defined to simplify the final expressions
#They are from the article Nagy, 1966
    
#+(x2_prism, y2_prism, z1_prism)
#-(x1_prism, y2_prism, z1_prism)
#-(x2_prism, y1_prism, z1_prism)
#+(x1_prism, y1_prism, z1_prism)
#-(x2_prism, y2_prism, z2_prism)
#+(x1_prism, y2_prism, z2_prism)
#+(x2_prism, y1_prism, z2_prism)
#-(x1_prism, y1_prism, z2_prism)
    
def prism_func2(x,y,z):

# Full expression is: x*np.log(y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) \
#                     + y*np.log(x + np.sqrt(np.abs(x)*np.abs(x)+np.abs(y)*np.abs(y)+z*z)) \
#                     - z*math.atan2(y*x, z*ra)
#
# The full expression is separated into the following parts
    radius=np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y) + np.abs(z)*np.abs(z))
    
    if y==0:
       expr_prism1= x*np.log(radius)
    else:
       expr_prism1= x*np.log(y + radius)
       
    if x==0:
       expr_prism2= y*np.log(radius)
    else:
       expr_prism2= y*np.log(x + radius)   

    expr_prism3= -z*math.atan2(y*x,z*radius)
    full_expr_prism=expr_prism1+expr_prism2+expr_prism3    

    if DEBUG==1:
        print('--------------',file=file)
        print('',file=file)
        print("The x, y and z coordinates used in this prism_func2() call are: x=", str("{0:<10.2f}".format(x)),"y=", str("{0:<10.2f}".format(y)),"and z=",str("{0:<10.2f}".format(z)), file=file)
        print('Prism func2 part A = '+str(expr_prism1),file=file)
        print('Prism func2 part B = '+str(expr_prism2),file=file)
        print('Prism func2 part C = '+str(expr_prism3),file=file)
        print('The sum of all 3 parts = '+str(full_expr_prism),file=file)
        print('',file=file)
        print('--------------',file=file)
        print('',file=file)
    return full_expr_prism

def prism_part_calc(xmax, xmin, ymax, ymin, zmax, zmin):
    if np.abs(ymin)<0.000001:   ymin=0.00001
    if np.abs(xmin)<0.000001:   xmin=0.00001
    if DEBUG==1:    print('doing part 1 ',file=file)
    part1=prism_func2(xmax, ymax, zmin)
    if DEBUG==1:    print('doing part 2 ',file=file)
    part2=prism_func2(xmin, ymax, zmin)
    if DEBUG==1:    print('doing part 3 ',file=file)
    part3=prism_func2(xmax, ymin, zmin)
    if DEBUG==1:    print('doing part 4 ',file=file)
    part4=prism_func2(xmin, ymin, zmin)
    if DEBUG==1:    print('doing part 5 ',file=file)
    part5=prism_func2(xmax, ymax, zmax)
    if DEBUG==1:    print('doing part 6 ',file=file)
    part6=prism_func2(xmin, ymax, zmax)
    if DEBUG==1:    print('doing part 7 ',file=file)
    part7=prism_func2(xmax, ymin, zmax)
    if DEBUG==1:    print('doing part 8 ',file=file)
    part8=prism_func2(xmin, ymin, zmax)
    
    if DEBUG==1:
        print('#################################################',file=file)
        print('                Prism part results',file=file)
        print('#################################################',file=file)
        print('',file=file)
        print("Prisma part results for coordinates x'max: "+str(xmax)+" x'min: "+str(xmin)+" y'max: "+str(ymax)+" y'min: "+str(ymin)+" z'max: "+str(zmax)+" z'min: "+str(zmin)+" :",file=file)
        print('',file=file)
        print('part 1 results (pos) '+str(part1),file=file)
        print('part 2 results (neg) '+str(part2),file=file)
        print('part 3 results (neg) '+str(part3),file=file)
        print('part 4 results (pos) '+str(part4),file=file)
        print('part 5 results (neg) '+str(part5),file=file)
        print('part 6 results (pos) '+str(part6),file=file)
        print('part 7 results (pos) '+str(part7),file=file)
        print('part 8 results (neg) '+str(part8),file=file)
        partpos=part1+part4+part6+part7
        partneg=part2+part3+part5+part8
        print('',file=file)
        print('sum pos parts ='+str(partpos),file=file)
        print('sum neg parts ='+str(partneg),file=file)
        print('',file=file)    
    return part1, part2, part3, part4, part5, part6, part7, part8


def prism_part_sum(prism_parts):
    #Based on prism_case() results, the sum is part1 - part2 - part3 + part4 - part5 + part6 + part7 - part8
    part1, part2, part3, part4, part5, part6, part7, part8 = prism_parts
    sum_prism= part1 - part2 - part3 + part4 - part5 + part6 + part7 - part8
    if DEBUG==1:
        print('',file=file)
        print('Sum of all the prism segments :'+ str(sum_prism),file=file)
        print('',file=file)    
    return sum_prism

def sub_prism_sum(subPrismEff):
    fullPrismEff=sum(subPrismEff)
    if DEBUG==1:
        print('',file=file)
        print('The initial prism was divided into sub-prisms.', file=file)
        print('Sum of all the sub-prisms :'+ str(fullPrismEff),file=file)
        print('',file=file)   
    return fullPrismEff

def prism_case1(x2, x1, y2, y1, z2, z1):
    prism_parts_results=prism_part_calc(x2, x1, y2, y1, z2, z1)
    prism_results=prism_part_sum(prism_parts_results)
    return prism_results

def reflect_yaxis(x2,x1):
    xx2=np.abs(x1)
    xx1=np.abs(x2)
    if DEBUG==1:
        print('',file=file)
        print('Reflect x across y-Axis',file=file)
        print('New x coordinates are : xmax = ',str(xx2),'and xmin = ',str(xx1),file=file)
        print('',file=file)
    return xx2,xx1

def reflect_xaxis(y2,y1):
    yy2=np.abs(y1)
    yy1=np.abs(y2)
    if DEBUG==1:
        print('',file=file)
        print('y Axis reflected',file=file)
        print('New y coordinates are : ymax = ',str(yy2),'and ymin = ',str(yy1),file=file)
        print('',file=file)
    return yy2,yy1
            
######################################
#PARAMETERS AND CONSTANTS
######################################
    
#Create an output file to store the console output
#orig_stdout = sys.stdout
#f_console = open('out.txt', 'w')
#sys.stdout = f_console
np.seterr(all='raise')
np.seterr(all='print')
start_time = time.time()
    
#Gavitational constant in m3 kg-1 s-2 (or N m2 kg-2)
G=6.67408*10**(-11)




########################################
#INPUT CHECK
########################################
#Input checks common to both routines

#np.set_printoptions(precision=4) #Set precision of the values
np.set_printoptions(suppress=True) #Disables scientific notation

#DEBUG FLAG CHECK
if DEBUG ==1:
    print("Debug mode is ON multiple files with intermediate calculation results for each stations will be created.")
#If the input for "routine" is different from BGPOLY or PRISMA, raise an error and exit
if routine != 'PRISMA' and routine != 'BGPOLY':
    sys.exit('The routine name is wrong, please using either "BGPOLY" or "PRISMA" as inputs, with quotations marks and full capital letters.')
#Check the total number of MP
total_stations=len(MPs_all_data)
#Case exception for only 1 MP 
if number_stations==1 and total_stations!=4:
    sys.exit("The total number of stations is not the same as indicated : input total is "+str(number_stations)+" , total in file is "+str(total_stations)+" please check your input")
elif number_stations != 1 and number_stations != total_stations:
    if MPs_all_data.shape==(4,):
        sys.exit("The total number of stations is not the same as indicated : input total is "+str(number_stations)+" , total in file is 1 please check your input")
    else:
        sys.exit("The total number of stations is not the same as indicated : input total is "+str(number_stations)+" , total in file is "+str(total_stations)+" please check your input")
        
#Separate the MP numbers and the MP (x,y,z) LCS coordinates into two different arrays
if number_stations==1:
    MPs=MPs_all_data[7:10] #Array with all the MPs coordinates
    station_name=MPs_all_data[0] #Array with all the stations names
    residual_ano=MPs_all_data[9]
else:
    MPs=MPs_all_data[:,7:10] #Array with all the MPs coordinates
    station_name=MPs_all_data[:,0] #Array with all the stations names
    residual_ano=MPs_all_data[:,9]
#Check that the MPs are within the user's defined study region
if number_stations==1:
    if MPs[0].max() > area_MP[0]:
        sys.exit("The following stations are out of boundary along the x axis, check your station or check your defined study area max x axis value.")
    if MPs[0].min() < area_MP[1]:
        sys.exit("The following stations are out of boundary along the x axis, check your station or check your defined study area min x axis value.")
    if MPs[1].max() > area_MP[2]:
        sys.exit("The following stations are out of boundary along the y axis, check your station or check your defined study area max y axis value.")
    if MPs[1].min() < area_MP[3]:
        sys.exit("The following stations are out of boundary along the y axis, check your station or check your defined study area min y axis value.")
    if MPs[2].max() > area_MP[4]:
        sys.exit("The following stations are out of boundary along the z axis, check your station or check your defined study area max z axis value.")
    if MPs[2].min() < area_MP[5]:
        sys.exit("The following stations are out of boundary along the z axis, check your station or check your defined study area min z axis value.")  
else:
    if MPs[:,0].max() > area_MP[0]:
        outofbx= np.where(MPs[:,0]==MPs[:,0].max())
        sys.exit("The following stations are out of boundary along the x axis, check station(s) "+str(outofbx[0])+" or check your defined study area max x axis value.")
    if MPs[:,0].min() < area_MP[1]:
        outofbx= np.where(MPs[:,0]==MPs[:,0].min())
        sys.exit("The following stations are out of boundary along the x axis, check station(s) "+str(outofbx[0])+" or check your defined study area min x axis value.")
    if MPs[:,1].max() > area_MP[2]:
        outofbx= np.where(MPs[:,1]==MPs[:,1].max())
        sys.exit("The following stations are out of boundary along the y axis, check station(s) "+str(outofbx[0])+" or check your defined study area max y axis value.")
    if MPs[:,1].min() < area_MP[3]:
        outofbx= np.where(MPs[:,1]==MPs[:,1].min())
        sys.exit("The following stations are out of boundary along the y axis, check station(s) "+str(outofbx[0])+" or check your defined study area min y axis value.")
    if MPs[:,2].max() > area_MP[4]:
        outofbx= np.where(MPs[:,2]==MPs[:,2].max())
        sys.exit("The following stations are out of boundary along the z axis, check station(s) "+str(outofbx[0])+" or check your defined study area max z axis value.")
    if MPs[:,2].min() < area_MP[5]:
        outofbx= np.where(MPs[:,2]==MPs[:,2].min())
        sys.exit("The following stations are out of boundary along the z axis, check station(s) "+str(outofbx[0])+" or check your defined study area min z axis value.")    
#If the folder where the run will be saved doesn't exist, it will be created 
if not os.path.exists(run_folder_path):
    os.makedirs(run_folder_path)    
#Input checks for PRISMA's routine
if routine=='PRISMA':
    print("Selected routine is 'PRISMA'")
    print("Data will be saved in the folder '"+str(run_folder_path)+"'")
    #Check the total number of prisms
    total_prisms=len(pri_data.T)
    if prisms_number==1 and pri_data.shape!=(7,):
        sys.exit("The total number of prisms isn't "+str(prisms_number)+" it is "+str(total_prisms)+" please check your input")
    elif prisms_number!=1 and total_prisms != prisms_number:
        if pri_data.shape==(7,):
            sys.exit("The total number of prisms isn't "+str(prisms_number)+" it is 1 please check your input")
        else:
            sys.exit("The total number of prisms isn't "+str(prisms_number)+" it is "+str(total_prisms)+" please check your input")
    elif prisms_number!=1 and pri_data.shape==(7,):
        sys.exit("The total number of prisms isn't "+str(prisms_number)+" it is 1 please check your input")
    #Check that no prisms are inside each other or that they have gaps between them
        ############To be added #############
    #Check if MPs are inside a prism using the elevation
    for i,j in itertools.product(range(prisms_number), range(number_stations)):
        if number_stations==prisms_number==1:
            if z1s <= MPs[2] :
                sys.exit("1 MP and 1 prism At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
        elif number_stations==1 and prisms_number!=1:
            if z1s[i] <= MPs[2] :
                sys.exit("1 MP At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
        elif number_stations!=1 and prisms_number==1:
            if z1s <= MPs[:,2][j] :
                sys.exit("1 prism At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
        else:
            if z1s[i] <= MPs[:,2][j] :
                sys.exit("At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
    print("No MP below the top of prism, all input is correct.")
    #Write input information in the output file
    run_day_time=datetime.datetime.now()
    
    #Getting the name of the model
    with open(PRISMA_input_file, 'r') as file_input:
        lines = file_input.read().splitlines()
        modelname=lines[0].split(sep=" ")[2]
    
    if modelname != model_name:
        sys.exit("The model name in the input file is "+modelname+", which is not the same as the name entered in the command file, "+model_name+", please check your input.")
    
    output_file_name=str(run_folder_path)+'/'+str(modelname)+'_PRISM_output.txt'
    
    file_prisma=open(str(output_file_name),'w')
    print("#PRISMA routine run date and start time:", file=file_prisma)
    print("#"+str("%s"%run_day_time), file=file_prisma)
    print("", file=file_prisma)
    print("#Run folder path",file=file_prisma)
    print("#"+str(run_folder_path), file=file_prisma)
    print("", file=file_prisma)
    print("#PRISMA input file name and path", file=file_prisma)
    print("#"+str(PRISMA_input_file), file=file_prisma)
    print("", file=file_prisma)
    print("#Study region xmax (m)  xmin(m)   ymax (m)  ymin(m)  zmax (m)  zmin(m) ", file=file_prisma)
    print("#"+str(area_MP), file=file_prisma)
    print("", file=file_prisma)
#Write the fixed header of the input file and the data from the control file, in the output file
    with open(PRISMA_input_file, 'r') as file_input:
        #Make a list with each of the lines of the input file
        lines = file_input.read().splitlines()
        #Write the header lines, except the gravity stations header, in the new output file
        for i in range(length_header-1):
            print(str(lines[i]), file=file_prisma)
    
#Header for the output data
    if prisms_number==1:
        print("#MP Station \t| CH coordinates(m) E,\tN,\tElev\t| Distance(m) along\taway from the profile\t| Distance check \t|" +
          " LCS coordinates(m) x,\ty,\tz\t| Bouguer anomaly(mGal)\t| Regional anomaly(mGal)\t| Residual anomaly(mGal)\t| Prism gravity effect (mGal)", file=file_prisma)
    else:
        print("#MP Station \t| CH coordinates(m) E,\tN,\tElev\t| Distance(m) along\taway from the profile\t| Distance check \t|" +
          " LCS coordinates(m) x,\ty,\tz\t| Bouguer anomaly(mGal)\t| Regional anomaly(mGal)\t| Residual anomaly(mGal)\t|"+
          " All prisms gravity effect (mGal)\t| Invididual prisms gravity effect (mGal)", file=file_prisma)
        
###### NEED TO UPDATE BGPOLY HEADER PRINT        
#Input checks for BGPoly's routine    
if routine=='BGPOLY':
    print("Selected routine is 'BGPOLY'")
    print("Data will be saved in the folder '"+str(run_folder_path)+"'")

    #Plot routine to check the input data in 3-D
      ###### To be added ######
        
    #Check and return the number of height lines. It must be an odd number GRAVI3D will exit if it is even
    #Select unique height lines values to store them
    height_lines=[]
    for i in range(len(z_poly)):
        if z_poly[i-1] != z_poly[i] : 
            height_lines.append(z_poly[i])
    #Sort height lines in an increasing order       
    height_lines=np.asarray(height_lines)
    height_lines.sort()
    #Check the total number of polygons
    if len(height_lines) != number_polygons:
        sys.exit("The total number of polygons isn't "+str(number_polygons)+" it is "+str(len(height_lines))+" please check your input")
    #Check if the total height lines number is odd or even
    if len(height_lines) % 2 != 0 :
        print('The number of height lines is '+str(len(height_lines)))
    else :
        print('The number of height lines is even '+str(len(height_lines))+', it has to be odd, the interpolation will not work otherwise.')
        sys.exit("The amount of height lines isn't odd")
    #Check if the total number of density values is correct
    if len(height_lines) != len(rho_poly):
        sys.exit("The total number of density values is "+str(len(rho_poly))+" it should be "+str(len(height_lines))+" please check your input")
    #Check if MPs are inside the polygons using the elevation
   #for i,j in itertools.product(range(len(height_lines)), range(total_stations)):
   #     if height_lines[i] <= MPs[:,2][j] :
   #         sys.exit("At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
    print("No MP below the top of prism, all input is correct.")          
    #Write input information in the output file
    run_day_time=datetime.datetime.now()
    file_bgpoly=open(str(run_folder_path)+'\BGpoly full results.txt','w')
    print("#BGPoly routine run date and start time:", file=file_bgpoly)
    print("#"+str("%s"%run_day_time), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Run folder path",file=file_bgpoly)
    print("#"+str(run_folder_path), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#BGPoly input file name and path", file=file_bgpoly)
    print("#"+str(BGPOLY_input_file), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Number of stations", file=file_bgpoly)
    print("#"+str(number_stations), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Study region xmax (m)  xmin(m)   ymax (m)  ymin(m)  zmax (m)  zmin(m) ", file=file_bgpoly)
    print("#"+str(area_MP), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Polygon coordinates and density file names", file=file_bgpoly)
    print("#"+str(BGPOLY_input_file), file=file_bgpoly)
    print("#"+str(density_poly_file_name), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Number of polygons", file=file_bgpoly)
    print("#"+str(number_polygons), file=file_bgpoly)
    print("", file=file_bgpoly)
    
#Header for the output data
    print("#MP Station | MP coordinates x, y, z | gravity effect (mGal)", file=file_bgpoly)

#######################################
#MAIN ROUTINE - BGPOLY
#######################################
#Store unique height line values, check that the number of height lines is odd and return to the use
if routine=='BGPOLY':        
    
    nb_MPs=number_stations
    iso_rho=np.asarray(rho_poly)    
    poly_coor=poly_data.T[:,:3] #Store the x,y,z coordinates of the polygons
    
    V=np.zeros(shape=(nb_MPs,len(height_lines))) # Store the height lines gravity effect
    #Makes sures that the numbers of column is the same as the number of interpolations needed, to store the value  
    num_lines=int(len(height_lines)//2) #For the interpolation - see above
    body_effect=np.zeros(shape=(nb_MPs,num_lines)) #Store the integration results
    gravity_value=np.zeros((len(body_effect))) #Store the disturbing body gravity effect

    for i in range(nb_MPs):
        
        #FOR TESTING - clear isoline effect list and store intermediate height line calculations for each new MP
        iso_effect_store=[] 
               
        #Announce which point is being processed
        if DEBUG==1:
            if nb_MPs == 1:
                file=open(str(run_folder_path)+'\BGPoly for MP'+str("{:.0f}".format(station_name))+'.txt','w')
                print(' ',file=file)
                print(' ',file=file)
                print('__________________________________________________',file=file)
                print(' ',file=file)
                print('The MP coordinates are '+str(MPs),file=file)
                print('__________________________________________________',file=file)
            else:
                file=open(str(run_folder_path)+'\BGPoly for MP'+str("{:.0f}".format(station_name[i]))+'.txt','w')
                print(' ',file=file)
                print(' ',file=file)
                print('__________________________________________________',file=file)
                print(' ',file=file)
                print('The MP coordinates are '+str(MPs[i]),file=file)
                print('__________________________________________________',file=file)
            
            #FOR TESTING - Return MP new coordinates, should always be [0,0,0]
            MPs_trans=np.ones(shape=(nb_MPs,3))
            if nb_MPs == 1:
                coor_transform(MPs,MPs,MPs_trans[i])
            else:
                coor_transform(MPs[i],MPs[i],MPs_trans[i])
            
            print(' ',file=file)
            print('The MP new coordinates are '+str(MPs_trans[i]),file=file)
            print(' ',file=file)
        
    #Create an array to store each height line coordinates
        
        for j in range(len(height_lines)):
            isoline=[]
            for l in range(len(poly_coor)):
                isoline,poly_list=poly_line(poly_coor[l],height_lines[j],isoline)
            isoline=np.asarray(isoline)
       
    #Transform the coordinates of the height line
            iso_trans=np.zeros(shape=(len(isoline),3))
            if nb_MPs == 1:
                for l in range(len(isoline)):
                    coor_transform(MPs,isoline[l],iso_trans[l])
            else:
                for l in range(len(isoline)):
                    coor_transform(MPs[i],isoline[l],iso_trans[l])
            x_iso=iso_trans[:,0]
            y_iso=iso_trans[:,1]
            z_iso=iso_trans[:,2]

    #Return the transformed coordinates of the height lines
            if DEBUG==1:
                print(' ',file=file)
                print('################################ START OF THE HEIGHT LINE CALCULATIONS ################################',file=file)
                print(' ',file=file)
                print('               -------               ',file=file)
                print(' ', file=file)
                print('The height line '+str(height_lines[j])+'(m) coordinates are:', file=file)
                print(' ', file=file)
                print(str(isoline), file=file)
                print(' ', file=file) 
                print(' ',file=file) 
                print('The height line '+str(height_lines[j])+'(m) new local coordinates are: ',file=file) 
                print(' ',file=file)
                print(str(iso_trans),file=file)
        
    #Search for values that are equal to 0 and change their value to 0.001
            for c in range(len(iso_trans)):
                for t in range(len(iso_trans[0])):
                    if iso_trans[c,t] == 0:
                        iso_trans[c,t]=0.001
                    if DEBUG==1:
                        print('A value 0 as been replaced by 0.001, in BGPoply input at position ['+str(c)+','+str(t)+']',file=file)
                    
    #Array to store the p,q,f,m,.... variables for each height line
            r_all=np.zeros(len(x_iso))
            sum_ri_all=np.zeros(len(x_iso)-1)
            p_all=np.zeros(len(x_iso)-1)
            m_all=np.zeros(len(x_iso)-1)
            q_all=np.zeros(len(x_iso)-1)
            f_all=np.zeros(len(x_iso)-1)
        
    #Calcul of the parameters r, sum_ri, p, m, q, f; 
    #if statement to allow the loop to do the calculation between the last and first element
            for l in range(len(x_iso)):
                r_all[l]=r_def(x_iso[l],y_iso[l])
            
            for l in range(len(x_iso)-1):
                sum_ri_all[l]=sum_ri_def(x_iso[l],x_iso[l+1],y_iso[l],y_iso[l+1])
            
    #The for loop duplicated, because r_all and sum_ri_all have to be calculated first
            for l in range(len(x_iso)-1):             
                p_all[l]=p_def(x_iso[l],x_iso[l+1],y_iso[l],y_iso[l+1],sum_ri_all[l])
                m_all[l]=m_def(x_iso[l],x_iso[l+1],y_iso[l],y_iso[l+1],r_all[l],r_all[l+1])
                q_all[l]=q_def(x_iso[l],x_iso[l+1],y_iso[l],y_iso[l+1],sum_ri_all[l],r_all[l])
                f_all[l]=f_def(x_iso[l],x_iso[l+1],y_iso[l],y_iso[l+1],sum_ri_all[l],r_all[l+1])
            
            height_effect=np.array([p_all,q_all,f_all,r_all,m_all])
            
            if DEBUG==1:
                print(' ',file=file)
                print('The values of p,q,f,r,m for the height line '+str(height_lines[j])+'(m) are: ',file=file)
                print(' ',file=file)
                print('sum_ri '+str(sum_ri_all), file=file)
                print('p '+str(height_effect[0]),file=file)
                print('q '+str(height_effect[1]),file=file)
                print('f '+str(height_effect[2]),file=file)
                print('r '+str(height_effect[3]),file=file)
                print('m '+str(height_effect[4]),file=file)
                print(' ',file=file)

    #Assign each variable it's array
            p=height_effect[0]
            q=height_effect[1]
            f=height_effect[2]
            r=height_effect[3]
            m=height_effect[4]
        
    #Calcul of the gravity anomaly of each height line
            G_E=np.ones(len(x_iso)-1)
            W=np.zeros(len(x_iso)-1)
            S=np.zeros(len(x_iso)-1)
        
            for l in range(len(height_effect[1])) :
    #Calculation of the terms inside of the sum for each points
                W[l]=W_sign(m[l])
                S[l]=S_sign(p[l])
                G_E[l]=iso_eff(z_iso[0], x_iso[l], x_iso[l+1], y_iso[l], y_iso[l+1], r[l], r[l+1], q[l], p[l], f[l],
                       W_sign(m[l]), S_sign(p[l]))
            if DEBUG==1:
                print(' ',file=file)
                print('W= '+str(W),file=file)
                print('S= '+str(W),file=file)
                print('Result (no unit) of the calculations between each points of height '+str(height_lines[j])+'(m): ',file=file)
                print(' ',file=file)
                print(G_E,file=file)
                print(' ',file=file)
        
            rho=iso_rho[j] #Select current isoline density
            
            #Use to check intermediary results of height lines effect  
            iso_effect_store.append(G_E)
            
            #Calcul the effect of the height lines on the MP, in s-2
            #Stores the results of each lamina for each MP : V[0,0] is the effect of the 1st lamina for the 1st MP
            V[i,j]=G*rho*sum(G_E)
            
            if DEBUG==1:
                print(' ',file=file)
                if nb_MPs == 1:
                    print('Height line '+str(height_lines[j])+'(m) effect (in s-2) on this MP:'+str(MPs)+':',file=file)
                else:
                    print('Height line '+str(height_lines[j])+'(m) effect (in s-2) on this MP:'+str(MPs[i])+':',file=file)
                print(' ',file=file)
                print(V[i,j],file=file)
                  

        if DEBUG==1:
            print(' ', file=file)
            print(' ', file=file)
            print('Effect of all polygons segment (no unit) on this MP:', file=file)
            print(' ', file=file)
            print(*iso_effect_store, sep="\n", file=file)
            print(' ', file=file)
    
            print(' ',file=file)
            if nb_MPs == 1:
                print('Effect of all height lines (in s-2) on this MP:'+str(MPs)+':',file=file)
                print(' ',file=file)
                print(V[i],file=file)
                print(' ',file=file)
        
                print(' ', file=file)
                print('End of '+str(MPs)+' calculations',file=file)
            else:
                print('Effect of all height lines (in s-2) on this MP:'+str(MPs[i])+':',file=file)
                print(' ',file=file)
                print(V[i],file=file)
                print(' ',file=file)
        
                print(' ', file=file)
                print('End of '+str(MPs[i])+' calculations',file=file)
            file.close()
        
    #Calculation of the gravity effect of the whole body
    
    #Does the integration and interpolation every 3 line
        for j,l in zip(range(1,len(height_lines),2), range(num_lines)):
            body_effect[i,l]=body_func(V[i,j-1],V[i,j],V[i,j+1],height_lines[j-1],height_lines[j],height_lines[j+1])
        
    #print(l,' '+str(MPs[i])+' height calc '+str(j)) #check that the loop is working properly
            
    #Sum the interpolated gravity effect of all the isolines and multiply by 10**3 to convert Gal to mGal
        gravity_value[i]=sum(body_effect[i])*10**5
    #Write the station coordinates and the disturbing body gravity effect on it
        np.set_printoptions(formatter={'float': '{: .4f}'.format})
        if nb_MPs == 1:
            print(str("{0:<10.0f}".format(i))," ",str("{0:<10.2f}".format(MPs[0]))," ",str("{0:<10.2f}".format(MPs[1]))," ",str("{0:<10.2f}".format(MPs[2]))," ","{0:<10.4f}".format(gravity_value[i]),file=file_bgpoly)
        else:
            print(str("{0:<10.0f}".format(i))," ",str("{0:<10.2f}".format(MPs[i,0]))," ",str("{0:<10.2f}".format(MPs[i,1]))," ",str("{0:<10.2f}".format(MPs[i,2]))," ","{0:<10.4f}".format(gravity_value[i]),file=file_bgpoly)
    
    if DEBUG==1:
        file_bgpoly.close()


#######################################
#MAIN ROUTINE - PRISMA
#######################################
if routine=='PRISMA':

    nb_MPs=number_stations
    nb_prisms=prisms_number
    
    prism_results=np.zeros(shape=(nb_MPs,nb_prisms))
    prism_effect=np.zeros(shape=(nb_MPs,nb_prisms))
    
    sum_prism_effect=np.zeros(nb_MPs)
    np.set_printoptions(formatter={'float': '{:.2f}'.format})
    
    #Select one MP at a time to do the calculation
    for i in range(nb_MPs):
        #Announce which point is being processed
        if DEBUG==1:
            
            if nb_MPs == 1:
                file=open(str(run_folder_path)+'\Prisma for MP'+str("{:.0f}".format(station_name[0]))+'.txt','w')
                print(' ',file=file)
                print(' ',file=file)
                print('__________________________________________________',file=file)
                print(' ',file=file)
                print('The MP coordinates are '+str(MPs),file=file)
                print('__________________________________________________',file=file)
            else:
                file=open(str(run_folder_path)+'\Prisma for MP'+str("{:.0f}".format(station_name[i]))+'.txt','w')
                print(' ',file=file)
                print(' ',file=file)
                print('__________________________________________________',file=file)
                print(' ',file=file)
                print('The MP coordinates are '+str(MPs[i]),file=file)
                print('__________________________________________________',file=file)
            
            #FOR TESTING - Return MP new coordinates, should always be [0,0,0]
            MPs_trans=np.ones(shape=(nb_MPs,3))
            if nb_MPs == 1:
                coor_transform(MPs,MPs,MPs_trans[i])
            else:
                coor_transform(MPs[i],MPs[i],MPs_trans[i])
            
            print(' ',file=file)
            print('The MP new coordinates are '+str(MPs_trans[i]),file=file)
            print(' ',file=file)
    
#Select one prism at a time
#Create an array with the coordinate couples of the selected prism
        for j in range(nb_prisms):
            if nb_prisms==1:
                coor_prism=np.array([[x2s, x1s],[y2s, y1s],[z2s, z1s]])
                rho_prism=rho_pri
                if DEBUG==1:
                    print('',file=file)
                    print('################################ START OF THE PRISM '+str(j+1)+' CALCULATIONS ################################',file=file)
                    print('',file=file)
                    print('Prism(s) coordinates: xmax: '+str(x2s)+', xmin: '+str(x1s)+
                          ', ymax: '+str(y2s)+', ymin: '+str(y1s)+
                          ', zmax: '+str(z2s)+', zmin: '+str(z1s),file=file)
                    print('Prism(s) density: '+str(rho_pri),file=file)
                    print('',file=file)
            else:
                coor_prism=np.array([[x2s[j],x1s[j]],[y2s[j],y1s[j]],[z2s[j],z1s[j]]])
                rho_prism=rho_pri[j]
                if DEBUG==1:
                    print('',file=file)
                    print('################################ START OF THE PRISM '+str(j+1)+' CALCULATIONS ################################',file=file)
                    print('',file=file)
                    print('Prism(s) coordinates: xmax: '+str(x2s[j])+', xmin: '+str(x1s[j])+
                          ', ymax: '+str(y2s[j])+', ymin: '+str(y1s[j])+
                          ', zmax: '+str(z2s[j])+', zmin: '+str(z1s[j]),file=file)
                    print('Prism(s) density: '+str(rho_pri[j]),file=file)
                    print('',file=file)
        
#Transform prisms coordinates with the MPs
            coor_trans=np.zeros(shape=(2,3))
            for k in range(len(coor_trans)):
                if nb_MPs == 1:
                    coor_transform(MPs,coor_prism.T[k],coor_trans[k])
                else:
                    coor_transform(MPs[i],coor_prism.T[k],coor_trans[k])
            coor_trans=coor_trans.T
                    
            if DEBUG==1:        
                print('',file=file)
                print('New local prism(s) coordinates: \n'
                      "x'max, x'min:",str(coor_trans[0])+'\n'
                      "y'max, y'min:",str(coor_trans[1])+'\n'
                      "z'max, z'min:",str(coor_trans[2]),file=file)
                
#Assign the new coordinates to new variables         
            x2_prism=coor_trans[0,0]
            x1_prism=coor_trans[0,1]
            y2_prism=coor_trans[1,0]
            y1_prism=coor_trans[1,1]
            z2_prism=coor_trans[2,0]
            z1_prism=coor_trans[2,1]

#Check the signs of the normalized coordinates and use the proper prism_case1() input
            if x1_prism > x2_prism:            
                xx1=x2_prism
                x2_prism=x1_prism
                x1_prism=xx1
                if DEBUG==1:
                    print('x1_prism was greater than x2_prism',file=file)
                    print('x-coordinates changed',file=file)
                    
            if y1_prism > y2_prism:
                yy1=y2_prism
                y2_prism=y1_prism
                y1_prism=yy1
                if DEBUG==1:
                    print('y1_prism was greater than y2_prism',file=file)
                    print('y-coordinates changed',file=file)
                    
#For x1 and x2, y1 and y2 of same signs, all positive
            if x1_prism >=0 :
                if y1_prism >=0:
                    if DEBUG==1:         
                        print('',file=file)
                        print('CASE 1',file=file)
                        print('x and y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                 
                    prism_results[i][j]=prism_case1(x2_prism, x1_prism, y2_prism, y1_prism, z2_prism, z1_prism)
                    #continue
                
 #For x1 and x2, y1 and y2 of same signs, x positive and y negative                      
                elif y2_prism <= 0:  
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 5b - Reflect the prism',file=file)
                        print('x signs are positive, y signs are negative, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)   
                    
                    #Reflect coordinates of the prism
                    yy2_prism, yy1_prism=reflect_xaxis(y2_prism,y1_prism)
                
                    prism_results[i][j]=prism_case1(x2_prism, x1_prism, yy2_prism, yy1_prism, z2_prism, z1_prism)
                    #continue
                
#For y1 and y2 of different signs and x1 and x2 of same sign
#The prism is separated in two, each part effect is calculated separately and then their effect is summed
                else : #y1_prism < 0 and y2_prism > 0
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 3a - Two sub-prisms, one reflected',file=file)
                        print('x signs are positive, y have different signs, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                
                    #Reflect coordinates of the prism
                    yy2_prism, yy1_prism=reflect_xaxis(0,y1_prism)
                    
                    #Effect of the first half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 3a sub-prism 1 (case1) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p1=prism_case1(x2_prism, x1_prism, y2_prism, 0, z2_prism, z1_prism)
                    
                    #Effect of the second half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 3a sub-prism 2 (case5b) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p2=prism_case1(x2_prism, x1_prism, yy2_prism, yy1_prism, z2_prism, z1_prism)
                    
                    #Sum for the total prism effect                    
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                #continue   
                 
            elif x1_prism < 0 and x2_prism <= 0:
                if y1_prism >=0 and y2_prism > 0:
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 5a - Reflect the prism',file=file)
                        print('x signs are negative, y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                        
                    #Reflect coordinates of the prism
                    xx2_prism, xx1_prism= reflect_yaxis(x2_prism,x1_prism)

                    prism_results[i][j]=prism_case1(xx2_prism, xx1_prism, y2_prism, y1_prism, z2_prism, z1_prism)
                    #continue
                
                elif y1_prism < 0 and y2_prism <= 0:
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 5c - Reflect the prism',file=file)
                        print('x and y signs are negative, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                    
                    #Reflect coordinates of the prism
                    xx2_prism, xx1_prism= reflect_yaxis(x2_prism,x1_prism)
                    yy2_prism, yy1_prism=reflect_xaxis(y2_prism,y1_prism)

                    prism_results[i][j]=prism_case1(xx2_prism, xx1_prism, yy2_prism, yy1_prism, z2_prism, z1_prism)
                    #continue
                
#For y1 and y2 of different signs and x1 and x2 of same sign
#The prism is separated in two, each part effect is calculated separately and then their effect is summed
                else : #y1_prism < 0 and y2_prism >0
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 3b - Two sub-prisms, one reflected',file=file)
                        print('x signs are negative, y have different signs, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                        
                    #Reflect coordinates of the prism
                    xx2_prism, xx1_prism= reflect_yaxis(x2_prism,x1_prism)
                    yy2_prism, yy1_prism=reflect_xaxis(0.0,y1_prism)
                    
                    #Effect of the first half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 3b sub-prism 1 (case 5a) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p1=prism_case1(xx2_prism, xx1_prism, y2_prism, 0, z2_prism, z1_prism)
                    
                    #Effect of the second half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 3b sub-prism 2 (case 5c) calculations ---',file=file)
                        print('',file=file)                                       
                    prism_results_p2=prism_case1(xx2_prism, xx1_prism, yy2_prism, yy1_prism, z2_prism, z1_prism)                 
                    
                    #Sum for the total prism effect 
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])                    
                #continue
            
#For x1 and x2 of different signs and y1 and y2 of same signs
#The prism is separated in two, each part effect is calculated separately and then their effect is summed
            elif x1_prism < 0 and x2_prism > 0:
                if y1_prism >=0 and y2_prism > 0:
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 2a - Two sub-prisms, one reflected twice, one once',file=file)
                        print('x have different signs, y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                        
                    #Reflect coordinates of the prism
                    xx2_prism, xx1_prism= reflect_yaxis(0.0,x1_prism)    
                    
                    #Effect of the first half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 2a sub-prism 1 (case1) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p1=prism_case1(x2_prism, 0, y2_prism, y1_prism, z2_prism, z1_prism)

                    #Effect of the second half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 2a sub-prism 2 (case 5a) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p2=prism_case1(xx2_prism, xx1_prism, y2_prism, y1_prism, z2_prism, z1_prism)
                    
                    #Sum for the total prism effect
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                    #continue
                    
                elif y1_prism < 0 and y2_prism <= 0:
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 2b - Two sub-prisms, one reflected twice, one once',file=file)
                        print('x have different signs, y signs are negative, x y and values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                        
                    #Reflect coordinates of the prism
                    xx2_prism, xx1_prism= reflect_yaxis(0.0,x1_prism)
                    yy2_prism, yy1_prism=reflect_xaxis(y2_prism,y1_prism)
                    
                    #Effect of the first half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 2b sub-prism 1 (case 5b) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p1=prism_case1(x2_prism, 0, yy2_prism, yy1_prism, z2_prism, z1_prism)

                    #Effect of the second half of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 2b sub-prism 2 (case 5c) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p2=prism_case1(xx2_prism, xx1_prism, yy2_prism, yy1_prism, z2_prism, z1_prism)

                    #Sum for the total prism effect
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                    #continue

#For x1 and x2, y1 and y2 of different signs
#The prism is separated in four, each part effect is calculated separately and then their effect is summed
                else : ##x1_prism < 0 and x2_prism > 0 and y1_prism < 0 and y2_prism > 0
                    if DEBUG==1:
                        print('',file=file)
                        print('CASE 4 - Four sub-prisms',file=file)
                        print('x and y have different signs, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                        print('',file=file)
                        
                    #Reflect coordinates of the prism
                    xx2_prism, xx1_prism= reflect_yaxis(0.0,x1_prism)
                    yy2_prism, yy1_prism=reflect_xaxis(0.0,y1_prism)
                        
                    #Effect of the first fourth of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 4 sub-prism 1 (case1) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p1=prism_case1(x2_prism, 0, y2_prism, 0, z2_prism, z1_prism)

                    #Effect of the second fourth of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 4 sub-prism 2 (case 5a) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p2=prism_case1(xx2_prism, xx1_prism, y2_prism, 0, z2_prism, z1_prism)
                    
                    #Effect of the third fourth of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 4 sub-prism 3 (case 5b) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p3=prism_case1(x2_prism, 0, yy2_prism, yy1_prism, z2_prism, z1_prism)
                    
                    #Effect of the last fourth of the prism
                    if DEBUG==1:
                        print('',file=file)
                        print('--- CASE 4 sub-prism 4 (case 5c) calculations ---',file=file)
                        print('',file=file)
                    prism_results_p4=prism_case1(xx2_prism, xx1_prism, yy2_prism, yy1_prism, z2_prism, z1_prism)
                    
                    #Sum for the total prism effect
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2,prism_results_p3,prism_results_p4])
                #continue
                    
            else :
                sys.exit("No matching case found, program stops for MP=",str(i), "and prism-nr=",str(j))
                    
            #Calculate and store the gravity effect of each prism on the station         
            prism_effect[i][j]=prism_results[:][i][j]*rho_prism*G*10**5
            if DEBUG==1 :
                print("Prism gravity effect = "+str(prism_effect[i][j]), file=file)
                print('',file=file)
                print('################################ END OF THE PRISM '+str(j+1)+' CALCULATIONS ################################',file=file)
            #Calculate and store the disturbing body gravity effect by summing the prisms effects
            sum_prism_effect[i]=sum(prism_effect[i])
        
        if DEBUG==1:
            file.close()
        
        #Write the stations coordinates and both the disturbing body effect and the individual prisms effects
        #For one station and one prism
        if nb_MPs==1 and nb_prisms==1:
            #Gravity stations, CH-coordinates and distance along/away from the gravity profile
            print(str("{0:<10.0f}".format(station_name)),"\t", "\t".join(str("{0:<10.3f}".format(mpchcoor)) for mpchcoor in MPs_all_data[1:6]),"\t", MP_distance_check,"\t",
                  #LCS coordinates used in this code
                  str("{0:<10.2f}".format(MPs[0])),"\t",str("{0:<10.2f}".format(MPs[1])),"\t",str("{0:<10.2f}".format(MPs[2])),"\t",
                  #Bouguer anomaly, regional anomaly and residual anomaly values
                  "\t".join(str("{0:<10.3f}".format(anomalies)) for anomalies in MPs_all_data[10:]),"\t",
                  #Modelled gravity effect for each station
                  "{0:<10.4f}".format(sum_prism_effect[i]),file=file_prisma)
            
        #For one station and multiple prisms
        elif nb_MPs==1 and nb_prisms!=1:
            print(str("{0:<10.0f}".format(station_name)),"\t", "\t".join(str("{0:<10.3f}".format(mpchcoor)) for mpchcoor in MPs_all_data[1:6]),"\t", MP_distance_check,"\t",
                  str("{0:<10.2f}".format(MPs[0])),"\t",str("{0:<10.2f}".format(MPs[1])),"\t",str("{0:<10.2f}".format(MPs[2])),"\t",
                  "\t".join(str("{0:<10.3f}".format(anomalies)) for anomalies in MPs_all_data[10:]),"\t",
                  #Modelled gravity effect for each station, whole body effect and single prism effect
                  "{0:<10.4f}".format(sum_prism_effect[i]),"\t", "\t".join(str("{0:<10.2f}".format(prism)) for prism in prism_effect[i]),file=file_prisma)
            
        #For multiple stations and one prism
        elif nb_MPs!=1 and nb_prisms==1:
            print(str("{0:<10.0f}".format(station_name[i])),"\t", "\t".join(str("{0:<10.3f}".format(mpchcoor)) for mpchcoor in MPs_all_data[i,1:6]),"\t", MP_distance_check[i],"\t",
                  str("{0:<10.2f}".format(MPs[i,0])),"\t",str("{0:<10.2f}".format(MPs[i,1])),"\t",str("{0:<10.2f}".format(MPs[i,2])),"\t",
                  "\t".join(str("{0:<10.3f}".format(anomalies)) for anomalies in MPs_all_data[i,10:]),"\t",
                  "{0:<10.4f}".format(sum_prism_effect[i]),file=file_prisma)
            
        #For multiple stations and multiple prisms
        else:
            print(str("{0:<10.0f}".format(station_name[i])),"\t", "\t".join(str("{0:<10.3f}".format(mpchcoor)) for mpchcoor in MPs_all_data[i,1:6]),"\t", MP_distance_check[i],"\t",
                  str("{0:<10.2f}".format(MPs[i,0])),"\t",str("{0:<10.2f}".format(MPs[i,1])),"\t",str("{0:<10.2f}".format(MPs[i,2])),"\t",
                  "\t".join(str("{0:<10.2f}".format(anomalies)) for anomalies in MPs_all_data[i,10:]),"\t",
                  "{0:<10.2f}".format(sum_prism_effect[i]),"\t", "\t".join(str("{0:<10.2f}".format(prism)) for prism in prism_effect[i]),file=file_prisma)
        
        #Write the rest of the input file in the output
    with open(PRISMA_input_file, 'r') as file_input:
        lines = file_input.read().splitlines()
        #Write the drillings data and the header for the prisms data
        for i in range(length_header+number_stations , prism_data_header):
            print(str(lines[i]), file=file_prisma)
        if nb_prisms==1:
            print(str("{0:<10.2f}".format(x2s)),"\t",str("{0:<10.2f}".format(x1s)),"\t",
                  str("{0:<10.2f}".format(y2s)),"\t",str("{0:<10.2f}".format(y1s)),"\t",
                  str("{0:<10.2f}".format(z2s)),"\t",str("{0:<10.2f}".format(z1s)),"\t",
                  str("{0:<10.2f}".format(rho_pri)),"\t",file=file_prisma)
        else:
        #Write the prisms data as used in this file - It should be the same as the input
            for j in range(nb_prisms):
                print(str("{0:<10.2f}".format(x2s[j])),"\t",str("{0:<10.2f}".format(x1s[j])),"\t",
                      str("{0:<10.2f}".format(y2s[j])),"\t",str("{0:<10.2f}".format(y1s[j])),"\t",
                      str("{0:<10.2f}".format(z2s[j])),"\t",str("{0:<10.2f}".format(z1s[j])),"\t",
                      str("{0:<10.2f}".format(rho_pri[j])),"\t",file=file_prisma)
  
            
            
    #Write individual prisms effects in a separate file - optional
    #np.savetxt(str(run_folder_path)+'\Prisma individual prisms results.txt', realest_effect, delimiter=' ', header='Invididual prisms gravity effect (mGal)')

#Close the output file after all the results are written
    file_prisma.close()

#Make a copy of the output file in the Display folder
shutil.copy2(str(output_file_name), "Display")
    
print("--- Gravi3D run finished in %s seconds ---" % (time.time() - start_time))

