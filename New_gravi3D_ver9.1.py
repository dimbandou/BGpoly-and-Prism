# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 12:53:39 2020

@author: Dimitri
"""

import numpy as np
import math
import itertools
import sys
import os
import time
import datetime
import fileinput


###########################
# LOADING USER INPUTS FILES
###########################
#Load the stations coordinates
MPs=np.loadtxt(MPs_file_name)
#Load the polygons/isolines/height lines coordinates for BGPOLY
if routine=='BGPOLY':
    poly_data=x_poly, y_poly, z_poly=np.loadtxt(poly_file_name, unpack=True)
#Load the polygons/isolines/height lines density values for BGPOLY
    rho_poly=np.loadtxt(density_poly_file_name)
#Load the prisms coordinates and density values for PRISMA
if routine=='PRISMA':
    pri_data=x2s, x1s, y2s, y1s, z2s, z1s, rho_pri=np.loadtxt(prism_file_name, unpack=True)


#FLAG FOR DEBUGGING MODE 0=OFF, 1=ON 
#CHANGE THIS FLAG ONLY IF YOU WANT ALL INTERMEDIATE CALCULATIONS RESULTS
#IT WILL CREATE INDIVIDUAL FILES FOR EACH STATIONS/MP 
DEBUG=1

#########################################
#FUNCTIONS DEFINITION
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

#For each MPs, normalize the coordinates so each MP is considered the origin.
#Make it so the height_line array is read and only x and y are normalized in it
#Normalize so the coordinate system is centered around the MPs considered as the origin
#this function takes an array with the coordinates of a MP and an array with the coordinates of an height line polygon
#it returns an array with the normalized coordinates of the height line
def normalize(mp,isoline,iso_norm):
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
        print('replaced arcsin value', file=file)
    else :
        G_E_arcsin= - np.arcsin(arcsin_exp1) + np.arcsin(arcsin_exp2)
        
    print('iso eff arcsin1 = '+str(arcsin_exp1),file=file)    
    print('iso eff arcsin2 = '+str(arcsin_exp2),file=file)
    print('diff arcsin '+str(G_E_arcsin), file= file)    
    return G_E_arcsin

def bgpoly_arccos(x_iso1, x_iso2, y_iso1, y_iso2, r, r_1):
    #Avoiding nan for np.arccos(x) because of float format where arccos_exp = 1.000000002
    arccos_exp=(x_iso1/r)*(x_iso2/r_1) + (y_iso1/r)*(y_iso2/r_1)
    
    if 1 < arccos_exp <= 1.000001 :
        G_E_arccos = 1
    elif -1.0001 < arccos_exp < -1:
        G_E_arccos = -1
    else :
        G_E_arccos = arccos_exp
    print('old arccos expression '+str(arccos_exp), file= file) 
    print('new arccos expression '+str(G_E_arccos), file= file)    
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

#Function to calculate the gravity effect of one prism, from Nagy, 1966.
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

def prism_func(x,y,z):
# Full expression is: x*np.log(y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) \
#                     + y*np.log(x + np.sqrt(np.abs(x)*np.abs(x)+np.abs(y)*np.abs(y)+z*z)) \
#                     - z*np.arcsin((z*z+y*y+y*np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) \
#                     / ((y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))*np.sqrt(np.abs(y)*np.abs(y)+z*z)))
# The full expression is separated into the following parts    
    expr_prism1= x*np.log(y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) 
    expr_prism2= y*np.log(x + np.sqrt(np.abs(x)*np.abs(x)+np.abs(y)*np.abs(y)+z*z)) 
    expr_prism3_1=(z*z+y*y+y*np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))
    expr_prism3_2=((y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))*np.sqrt(np.abs(y)*np.abs(y)+z*z))
    prism3_diff=np.abs(expr_prism3_1-expr_prism3_2)
#The full arcsin expression is: z*np.arcsin((z*z+y*y+y*np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))/((y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))*np.sqrt(np.abs(y)*np.abs(y)+z*z)))
    #Avoiding np.arcsin(x) with x>1 because of float format
    if prism3_diff <= 0.0001 :
        expr_prism3=z*np.arcsin(1)
    else:
        expr_prism3=z*np.arcsin(expr_prism3_1/expr_prism3_2)
        
    full_expr_prism=expr_prism1+expr_prism2-expr_prism3
    if DEBUG==1:
        print('',file=file)
        print('Prisma expression part 1 = '+str(expr_prism1),file=file)
        print('Prisma expression part 2 = '+str(expr_prism2),file=file)
        print('Prisma expression part 3 = '+str(expr_prism3),file=file)
        print('Prisma expression part 3_1 = '+str(expr_prism3_1),file=file)
        print('Prisma expression part 3_2 = '+str(expr_prism3_2),file=file)
        print('',file=file)
    return full_expr_prism

def prism_case(xmax, xmin, ymax, ymin, zmax, zmin):
    part1=prism_func(xmax, ymax, zmin)
    part2=prism_func(xmin, ymax, zmin)
    part3=prism_func(xmax, ymin, zmin)
    part4=prism_func(xmin, ymin, zmin)
    part5=prism_func(xmax, ymax, zmax)
    part6=prism_func(xmin, ymax, zmax)
    part7=prism_func(xmax, ymin, zmax)
    part8=prism_func(xmin, ymin, zmax)
    if DEBUG==1:
        print('',file=file)
        print('Prisma case results for coordinates xmax: '+str(xmax)+' xmin: '+str(xmin)+' ymax: '+str(ymax)+' ymin: '+str(ymin)+' zmax: '+str(zmax)+' zmin: '+str(zmin)+' :',file=file)
        print('',file=file)
        print('part 1 results '+str(part1),file=file)
        print('part 2 results '+str(part2),file=file)
        print('part 3 results '+str(part3),file=file)
        print('part 4 results '+str(part4),file=file)
        print('part 5 results '+str(part5),file=file)
        print('part 6 results '+str(part6),file=file)
        print('part 7 results '+str(part7),file=file)
        print('part 8 results '+str(part8),file=file)
        print('',file=file)    
    return part1, part2, part3, part4, part5, part6, part7, part8


def prism_part_sum(prism_parts):
    #Based on prism_case() results, the sum is part1 - part2 - part3 + part4 - part5 + part6 + part7 - part8
    sum_prism= prism_parts[0] - prism_parts[1] - prism_parts[2] + prism_parts[3] - prism_parts[4] + prism_parts[5] + prism_parts[6] - prism_parts[7]
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
######################################
#PARAMETERS AND CONSTANTS
######################################
    
#Create an output file to store the console output
orig_stdout = sys.stdout
f_console = open('out.txt', 'w')
sys.stdout = f_console
np.seterr(all='raise')
np.seterr(all='print')
start_time = time.time()
    
#Gavitational constant in m3 kg-1 s-2 (or N m2 kg-2)
G=6.67408*10**(-11)




########################################
#INPUT CHECK
########################################
#Input checks common to both routines

#DEBUG FLAG CHECK
if DEBUG ==1:
    print("Debug mode is ON multiple files with intermediate calculation results for each stations will be created.")
#If the input for "routine" is different from BGPOLY or PRISMA, raise an error and exit
if routine != 'PRISMA' and routine != 'BGPOLY':
    sys.exit('The routine name is wrong, please using either "BGPOLY" or "PRISMA" as inputs, with quotations marks and full capital letters.')
#Check the total number of MP
total_stations=len(MPs)
if number_stations != total_stations:
    sys.exit("The total number of stations is not the same as indicated : input total is "+str(sum(number_stations))+" , total in file is "+str(total_stations)+" please check your input")
#Check that the MPs are within the user's defined study region
if MPs[:,0].max() > area_MP[0]:
    outofbx= np.where(MPs[:,0]==MPs[:,0].max())
    sys.exit("The following stations are out of boundary along the x axis, check station(s) "+str(outofbx[0])+" or check your defined study area max x axis value.")
if MPs[:,0].min() < area_MP[1]:
    outofbx= np.where(MPs[:,0]==MPs[:,0].min())
    sys.exit("The following stations are out of boundary along the x axis, check station(s) "+str(outofbx[0])+" or check your defined study area min x axis value.")
if MPs[:,1].max() > area_MP[2]:
    outofbx= np.where(MPs[:,1]==MPs[:,1].max())
    sys.exit("The following stations are out of boundary along the y axis check station(s) "+str(outofbx[0])+" or check your defined study area max y axis value.")
if MPs[:,1].min() < area_MP[3]:
    outofbx= np.where(MPs[:,1]==MPs[:,1].min())
    sys.exit("The following stations are out of boundary along the y axis check station(s) "+str(outofbx[0])+" or check your defined study area min y axis value.")
if MPs[:,2].max() > area_MP[4]:
    outofbx= np.where(MPs[:,1]==MPs[:,1].max())
    sys.exit("The following stations are out of boundary along the z axis check station(s) "+str(outofbx[0])+" or check your defined study area max z axis value.")
if MPs[:,2].min() < area_MP[5]:
    outofbx= np.where(MPs[:,1]==MPs[:,1].min())
    sys.exit("The following stations are out of boundary along the z axis check station(s) "+str(outofbx[0])+" or check your defined study area min z axis value.")    
#If the folder where the run will be saved doesn't exist, it will be created 
if not os.path.exists(run_folder_name):
    os.makedirs(run_folder_name)    
#Input checks for PRISMA's routine
if routine=='PRISMA':
    print("Selected routine is 'PRISMA'")
    print("Data will be saved in the folder '"+str(run_folder_name)+"'")
    #Check the total number of prisms
    total_prisms=len(x2s)
    if total_prisms!= prisms_number:
        sys.exit("The total number of prisms isn't "+str(prisms_number)+" it is "+str(total_prisms)+" please check your input")
    #Check that no prisms are inside each other or that they have gaps between them
        ############To be added #############
    #Check if MPs are inside a prism using the elevation
    for i,j in itertools.product(range(total_prisms), range(total_stations)):
        if z1s[i] <= MPs[:,2][j] :
            sys.exit("At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
    print("No MP below the top of prism, all input is correct.")
    #Write input information in the output file
    run_day_time=datetime.datetime.now()
    file_prisma=open(str(run_folder_name)+'\Prisma full results.txt','w')
    print("#Run date and start time:", file=file_prisma)
    print("#"+str("%s"%run_day_time), file=file_prisma)
    print("", file=file_prisma)
    print("#Run name",file=file_prisma)
    print("#"+str(run_folder_name), file=file_prisma)
    print("", file=file_prisma)
    print("#Routine", file=file_prisma)
    print("#"+str(routine), file=file_prisma)
    print("", file=file_prisma)
    print("#Stations coordinates file name", file=file_prisma)
    print("#"+str(MPs_file_name), file=file_prisma)
    print("", file=file_prisma)
    print("#Distributed stations", file=file_prisma)
    print("#"+str(grid_MP), file=file_prisma)
    print("", file=file_prisma)
    print("#Number of stations", file=file_prisma)
    print("#"+str(number_stations), file=file_prisma)
    print("", file=file_prisma)
    print("#Study region xmax (m)  xmin(m)   ymax (m)  ymin(m)  zmax (m)  zmin(m) ", file=file_prisma)
    print("#"+str(area_MP), file=file_prisma)
    print("", file=file_prisma)
    print("#Prism coordinates and density file name", file=file_prisma)
    print("#"+str(prism_file_name), file=file_prisma)
    print("", file=file_prisma)
    print("#Number of prisms", file=file_prisma)
    print("#"+str(prisms_number), file=file_prisma)
    print("", file=file_prisma)
    
#Header for the output data
    print("#MP coordinates x, y, z | All prisms gravity effect (mGal) | Invididual prisms gravity effect (mGal)", file=file_prisma)    

#Input checks for BGPoly's routine    
if routine=='BGPOLY':
    print("Selected routine is 'BGPOLY'")
    print("Data will be saved in the folder '"+str(run_folder_name)+"'")

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
    total_poly=len(height_lines)
    if total_poly!= number_polygons:
        sys.exit("The total number of polygons isn't "+str(number_polygons)+" it is "+str(total_poly)+" please check your input")
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
    for i,j in itertools.product(range(len(height_lines)), range(total_stations)):
        if height_lines[i] <= MPs[:,2][j] :
            sys.exit("At this moment, this code will stop if any MP is below the top of a prism. Please check your MP number "+str(j))
    print("No MP below the top of prism, all input is correct.")          
    #Write input information in the output file
    run_day_time=datetime.datetime.now()
    file_bgpoly=open(str(run_folder_name)+'\BGpoly full results.txt','w')
    print("#Run date and start time:", file=file_bgpoly)
    print("#"+str("%s"%run_day_time), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Run name",file=file_bgpoly)
    print("#"+str(run_folder_name), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Routine", file=file_bgpoly)
    print("#"+str(routine), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Stations coordinates file name", file=file_bgpoly)
    print("#"+str(MPs_file_name), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Distributed stations", file=file_bgpoly)
    print("#"+str(grid_MP), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Number of stations", file=file_bgpoly)
    print("#"+str(number_stations), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Study region xmax (m)  xmin(m)   ymax (m)  ymin(m)  zmax (m)  zmin(m) ", file=file_bgpoly)
    print("#"+str(area_MP), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Polygon coordinates and density file names", file=file_bgpoly)
    print("#"+str(poly_file_name), file=file_bgpoly)
    print("#"+str(density_poly_file_name), file=file_bgpoly)
    print("", file=file_bgpoly)
    print("#Number of polygons", file=file_bgpoly)
    print("#"+str(number_polygons), file=file_bgpoly)
    print("", file=file_bgpoly)
    
#Header for the output data
    print("#MP coordinates x, y, z | gravity effect (mGal)", file=file_bgpoly)

#######################################
#MAIN ROUTINE - BGPOLY
#######################################
#Store unique height line values, check that the number of height lines is odd and return to the use
if routine=='BGPOLY':
        
    
    iso_rho=[]
    nb_MPs=len(MPs)
    os.makedirs(run_folder_name, exist_ok=True)
    poly_coor=poly_data.T[:,:3] #Store the x,y,z coordinates of the polygons
    iso_rho=np.asarray(rho_poly)
    V=np.zeros(shape=(nb_MPs,len(height_lines))) # Store the height lines gravity effect
    #Makes sures that the numbers of column is the same as the number of interpolations needed, to store the value  
    num_lines=int(len(height_lines)//2) #For the interpolation - see above
    body_effect=np.zeros(shape=(nb_MPs,num_lines)) #Store the integration results
    gravity_value=np.zeros((len(body_effect))) #Store the disturbing body gravity effect

    for i in range(nb_MPs):
        
        #FOR TESTING - clear isoline effect list and store intermediate height line calculations for each new MP
        iso_effect_store=[] 
        
        
        file=open(str(run_folder_name)+'\BGPoly for MP'+str(i)+'.txt','w')
        #Announce which point is being processed
        print(' ',file=file)
        print(' ',file=file)
        print('__________________________________________________',file=file)
        print(' ', file=file)
        print('The MP coordinates are '+str(MPs[i]),file=file)
        print('__________________________________________________',file=file)
        
        #FOR TESTING - Return MP new coordinates, should always be [0,0,0]
        MPs_norm=np.ones(shape=(len(MPs),3))
        normalize(MPs[i],MPs[i],MPs_norm[i]) #Temporary, for testing purposes
        
        print(' ',file=file)
        print('The MP new coordinates are '+str(MPs_norm[i]),file=file)
        print(' ',file=file)
        
    #Create an array to store each height line coordinates
        
        for j in range(len(height_lines)):
            isoline=[]
            for l in range(len(poly_coor)):
                isoline,poly_list=poly_line(poly_coor[l],height_lines[j],isoline)
            isoline=np.asarray(isoline)
       
    #Normalize the coordinates of the height line
            iso_norm=np.zeros(shape=(len(isoline),3))
            for l in range(len(isoline)):
                normalize(MPs[i],isoline[l],iso_norm[l])
            x_iso=iso_norm[:,0]
            y_iso=iso_norm[:,1]
            z_iso=iso_norm[:,2]

    #Return the normalized coordinates of the height lines
            print(' ',file=file)
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
            print(str(iso_norm),file=file)
        
    #Search for values that are equal to 0 and change their value to 0.001
            for c in range(len(iso_norm)):
                for t in range(len(iso_norm[0])):
                    if iso_norm[c,t] == 0:
                        iso_norm[c,t]=0.001
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
        
            print(' ',file=file)
            print('The values of p,q,f,r,m for the height line '+str(height_lines[j])+'(m) are: ',file=file)
            print(' ',file=file)
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
            print(' ',file=file)
            print('Height line '+str(height_lines[j])+'(m) effect (in s-2) on this MP:'+str(MPs[i])+':',file=file)
            print(' ',file=file)
            print(V[i,j],file=file)
                  

        
        print(' ', file=file)
        print(' ', file=file)
        print('Effect of all polygons segment (no unit) on this MP:', file=file)
        print(' ', file=file)
        print(*iso_effect_store, sep="\n", file=file)
        print(' ', file=file)
    
        print(' ',file=file)
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
        gravity_value[i]=sum(body_effect[i])*10**3
    #Write the station coordinates and the disturbing body gravity effect on it
        print(str(MPs[i,0])," ",str(MPs[i,1])," ",str(MPs[i,2])," ", gravity_value[i],file=file_bgpoly)
        
    file_bgpoly.close()


#######################################
#MAIN ROUTINE - PRISMA
#######################################
if routine=='PRISMA':
    #Normalize prism coordinates and set each MP as the origin
    #store all the prisms in one array
    nb_MPs=len(MPs)
    nb_prisms=x2s.size

    prism_results=np.zeros(shape=(nb_MPs,nb_prisms))
    prism_effect=np.zeros(shape=(nb_MPs,nb_prisms))
    
    sum_prism_effect=np.zeros(nb_MPs)
    
    
    #Select one MP at a time to do the calculation
    for i in range(nb_MPs):
        #Announce which point is being processed
        
        file=open(str(run_folder_name)+'\Prisma for MP'+str(i)+'.txt','w')
        
        print(' ',file=file)
        print(' ',file=file)
        print('__________________________________________________',file=file)
        print(' ',file=file)
        print('The MP coordinates are '+str(MPs[i]),file=file)
        print('__________________________________________________',file=file)
    
#Select one prism at a time
#Create an array with the coordinate couples of the selected prism
        for j in range(nb_prisms):
            coor_prism=np.array([[x2s[j],x1s[j]],[y2s[j],y1s[j]],[z2s[j],z1s[j]]])
            rho_prism=rho_pri[j]
            print('',file=file)
            print('Prism(s) coordinates: '+str(x2s[j])+' '+str(x1s[j])+' '+str(y2s[j])+' '+str(y1s[j])+' '+str(z2s[j])+' '+str(z1s[j]),file=file)
            print('Prism(s) density: '+str(rho_pri[j]),file=file)
            print('',file=file)
        
#Normalize prisms coordinates with the MPs
            coor_norm=np.zeros(shape=(2,3))
            for k in range(len(coor_norm)):
                normalize(MPs[i],coor_prism.T[k],coor_norm[k])
            coor_norm=coor_norm.T
        
#Check for any 0 in the normalized coordinates
            for k,l in itertools.product(range(len(coor_norm)), range(2)):
                if coor_norm[k,l]==0:
                    coor_norm[k,l]=0.001
            print('',file=file)
            print('New local prism(s) coordinates: '+str(coor_norm),file=file)
            print('',file=file)
        
            x2_prism=coor_norm[0,0]
            x1_prism=coor_norm[0,1]
            y2_prism=coor_norm[1,0]
            y1_prism=coor_norm[1,1]
            z2_prism=coor_norm[2,0]
            z1_prism=coor_norm[2,1]

#Check the signs of the normalized coordinates  and use proper equation
#For x1 and x2, y1 and y2 of same signs
            if x1_prism >=0 and x2_prism >=0:
                if y1_prism >=0 and y2_prism >=0:
                    prism_parts_xy_pos=prism_case(x2_prism, x1_prism, y2_prism, y1_prism, z2_prism, z1_prism)
                    
                    prism_results[i][j]=prism_part_sum(prism_parts_xy_pos)
                    
                    print('',file=file)
                    print('x and y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
                
                elif y1_prism < 0 and y2_prism < 0:
                    prism_parts_xy_pos_neg=prism_case(x2_prism, x1_prism, -y1_prism, -y2_prism, z2_prism, z1_prism)

                    prism_results[i][j]=prism_part_sum(prism_parts_xy_pos_neg)
                    
                    print('',file=file)
                    print('x signs are positive, y signs are negative, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
#For y1 and y2 of different signs and x1 and x2 of same sign
#The prism is separated in two, each part effect is calculated separately and then their effect is summed
                else :
                    #Effect of the first half of the prism
                    prism_parts_xy_pos_dif_p1=prism_case(x2_prism, x1_prism, y2_prism, 0, z2_prism, z1_prism)

                    prism_results_p1=prism_part_sum(prism_parts_xy_pos_dif_p1)
                    
                    #Effect of the second half of the prism
                    prism_parts_xy_pos_dif_p2=prism_case(x2_prism, x1_prism, -y1_prism, 0, z2_prism, z1_prism)
                                        
                    prism_results_p2=prism_part_sum(prism_parts_xy_pos_dif_p2)
                    
                    #Sum for the total prism effect                    
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                                        
                    print('',file=file)
                    print('x signs are positive, y have different signs, y values are '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
                                
            elif x1_prism < 0 and x2_prism < 0:
                if y1_prism >=0 and y2_prism >=0:
                    prism_parts_xy_neg_pos=prism_case(-x1_prism, -x2_prism, y2_prism, y1_prism, z2_prism, z1_prism)

                    prism_results[i][j]=prism_part_sum(prism_parts_xy_neg_pos)

                    print('',file=file)
                    print('x signs are negative, y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
                        
                elif y1_prism < 0 and y2_prism < 0:
                    prism_parts_xy_neg=prism_case(-x1_prism, -x2_prism, -y1_prism, -y2_prism, z2_prism, z1_prism)

                    prism_results[i][j]=prism_part_sum(prism_parts_xy_neg)

                    print('',file=file)
                    print('x and y signs are negative, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
#For y1 and y2 of different signs and x1 and x2 of same sign
#The prism is separated in two, each part effect is calculated separately and then their effect is summed
                else :
                    #Effect of the first half of the prism
                    prism_parts_xy_neg_dif_p1=prism_case(-x1_prism, -x2_prism, y2_prism, 0, z2_prism, z1_prism)

                    prism_results_p1=prism_part_sum(prism_parts_xy_neg_dif_p1)
                    
                    #Effect of the second half of the prism
                    prism_parts_xy_neg_dif_p2=prism_case(-x1_prism, -x2_prism, -y1_prism, 0, z2_prism, z1_prism)
                                        
                    prism_results_p2=prism_part_sum(prism_parts_xy_neg_dif_p2)                    
                    
                    #Sum for the total prism effect 
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                    

                    print('',file=file)
                    print('x signs are negative, y have different signs, y values are '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
                    
#For x1 and x2 of different signs and y1 and y2 of same signs
#The prism is separated in two, each part effect is calculated separately and then their effect is summed
            else :
                if y1_prism >=0 and y2_prism >=0:
                    #Effect of the first half of the prism
                    prism_parts_xy_dif_pos_p1=prism_case(x2_prism, 0, y2_prism, y1_prism, z2_prism, z1_prism)
                                       
                    prism_results_p1=prism_part_sum(prism_parts_xy_dif_pos_p1)

                    #Effect of the second half of the prism
                    prism_parts_xy_dif_pos_p2=prism_case(-x1_prism, 0, y2_prism, y1_prism, z2_prism, z1_prism)
                                        
                    prism_results_p2=prism_part_sum(prism_parts_xy_dif_pos_p2)
                    
                    #Sum for the total prism effect
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                    
                    print('',file=file)
                    print('x have different signs, y signs are positive, x values are '+str(x2_prism)+' '+str(x1_prism),file=file)
                    print('',file=file)
                    
                elif y1_prism < 0 and y2_prism < 0:
                    #Effect of the first half of the prism
                    prism_parts_xy_dif_neg_p1=prism_case(x2_prism, 0, -y1_prism, -y2_prism, z2_prism, z1_prism)
                                       
                    prism_results_p1=prism_part_sum(prism_parts_xy_dif_neg_p1)

                    #Effect of the second half of the prism
                    prism_parts_xy_dif_neg_p2=prism_case(-x1_prism, 0, -y1_prism, -y2_prism, z2_prism, z1_prism)
                                        
                    prism_results_p2=prism_part_sum(prism_parts_xy_dif_neg_p2)

                    #Sum for the total prism effect
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2])
                    

                    print('',file=file)
                    print('x have different signs, y signs are negative, x values are '+str(x2_prism)+' '+str(x1_prism),file=file)
                    print('',file=file)
#For x1 and x2, y1 and y2 of different signs
#The prism is separated in four, each part effect is calculated separately and then their effect is summed
                else :
                    #Effect of the first fourth of the prism
                    prism_parts_xy_dif_p1=prism_case(x2_prism, 0, y2_prism, 0, z2_prism, z1_prism)
                                       
                    prism_results_p1=prism_part_sum(prism_parts_xy_dif_p1)

                    #Effect of the second fourth of the prism
                    prism_parts_xy_dif_p2=prism_case(-x1_prism, 0, y2_prism, 0, z2_prism, z1_prism)
                                       
                    prism_results_p2=prism_part_sum(prism_parts_xy_dif_p2)
                    
                    #Effect of the third fourth of the prism
                    prism_parts_xy_dif_p3=prism_case(x2_prism, 0, -y1_prism, 0, z2_prism, z1_prism)
                                       
                    prism_results_p3=prism_part_sum(prism_parts_xy_dif_p3)
                    
                    #Effect of the last fourth of the prism
                    prism_parts_xy_dif_p4=prism_case(-x1_prism, 0, -y1_prism, 0, z2_prism, z1_prism)
                                       
                    prism_results_p4=prism_part_sum(prism_parts_xy_dif_p4)
                    
                    #Sum for the total prism effect
                    prism_results[i][j]=sub_prism_sum([prism_results_p1,prism_results_p2,prism_results_p3,prism_results_p4])
                    
                    
                    print('',file=file)
                    print('x and y have different signs, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                    print('',file=file)
                    
            #Calculate and store the gravity effect of each prism on the station         
            prism_effect[i][j]=prism_results[:][i][j]*rho_pri[j]*G*10**3
            print("Prism gravity effect = "+str(prism_effect[i][j]), file=file)
            #Calculate and store the disturbing body gravity effect by summing the prisms effects
            sum_prism_effect[i]=sum(prism_effect[i])
            
        file.close()
        
        #Write the stations coordinates and both the disturbing body effect and the individual prisms effects
        print(str(MPs[i,0])," ",str(MPs[i,1])," ",str(MPs[i,2])," ",sum_prism_effect[i]," ", "    ".join(str(prism) for prism in prism_effect[i]),file=file_prisma)
        
    #Write individual prisms effects in a separate file - optional
    #np.savetxt(str(run_folder_name)+'\Prisma individual prisms results.txt', realest_effect, delimiter=' ', header='Invididual prisms gravity effect (mGal)')

    file_prisma.close()
    
print("--- Gravi3D run finished in %s seconds ---" % (time.time() - start_time))
