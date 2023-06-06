# -*- coding: utf-8 -*-
"""
@author: Bandou
"""

##################
# USER INPUTS
#################

#Path to the folder where the results will be stored, the folder will be created if it doesn't exist
run_folder_path="Data_files"

#Switch to run either 'BGPOLY' or 'PRISMA' - This input is case sensitive
routine='PRISMA'

#Name of the model this input is written for
model_name="Bern1_model1"

#Path to the PRISMA input file
PRISMA_input_file="Data_files/"+str(model_name)+"_PRISM_input.txt"

#Study region with MPs relative to origin of coordinate system xmax (m)  xmin(m)   ymax (m)  ymin(m)  zmax (m)  zmin(m)
#Z axis is positive DOWNWARD
area_MP= 10000.0, -8000.0, 10000.0, -8000.0, 4000.0, -1000.0  

#This command executes Gravi3D 
exec(open("../../Gravi3D_ver10.4.8.py").read())


#Path to the DEM profile file
DEM_file="Preparation/DEM_bern1_PRISM.txt"
#Path to the bedrock profile file
bedrock_file="Preparation/bedrock_bern1_PRISM.txt"

#Plot all the prisms on the cross-section plots, even those that do not cross the gravity profile
#'yes' or 'no' input
all_prisms='yes'

#Project the data parallel to the valley direction, onto the gravity profile, for the cross-section plot
#Only use this if you have a large angle between your gravity profile and LCS x-axis
#'yes' or 'no' input
perpendicular_plot="no"

#Axes limits for the plot in meters and y axis top in mGal
xlimits_plot=[-100,3800]
ylimits_plot_top=[-2.75,0.25]
ylimits_plot_bot=[200,660]

#Axes limits for the 2D map plot in kilometers (GIS profile map boundaries)
xlimits_2D=[2596.065360,2601.254035]
ylimits_2D=[1196.735201,1201.752615]

xlimits_lcs_2D=[-2.5,2.5]
ylimits_lcs_2D=[-2.5,2.5]

#Label position for the density contrast values
xrho=-50
yrho=-2.1

#Label position for the directions on the cross-section plot
left_label='SW'
left_side=-2.25

right_label='NE'
right_side=-2.25

#Model name label height
name_height=620

#This command executes Gravi3D 
exec(open("../../Display_PRISMA_results.py", encoding='utf-8').read())