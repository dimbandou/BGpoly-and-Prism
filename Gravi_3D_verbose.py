# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:37:51 2019

@author: Bandou
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:11:38 2018

@author: Bandou
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib import cm
from scipy.interpolate import griddata
import sys
import os
import time
import csv

"""
                                    ### GRAVI 3D ###
                                ### CODE DESCRIPTION ###

This code purpose is to model the geometry of a subsurface body, using its gravity effect on measurement points (MPs);
It can be used with a grid of points or with a set of MPs, depending on the user input data.

To do so it allows the user to chose between two methods :
    - The 'BGPoly' method, from Talwani & Ewing 1960
    - The 'Prisma' method, from Nagy 1966

The first method will allow the user to use height lines corresponding to the geometry of the body.
Calculation of the gravity effect is done by summing the effect of each height line on the MPs,
and by interpolating between each height line.
Require the density for each height line.

the second method let the user calculate the gravity effect on the MPs, by approximating the shape of the subsurface bodies,
with prisms. Calculation of the gravity effect is done by summing the effect of each prisms on the MPs.
This method doesn't involve interpolation.
Each prism is made of 8 points, in a (x,y,z) cartesian system. 
Require the lithology to be known and the density of each layer.

The code works the following way: 
    
    -Ask the user for grid or set of points;
        function(s) involved: data_type
        
    -Read the input files: MPs coordinates;
        function(s) involved: None, check main routine 
        
    -Ask the user if they want to use the BGPoly or Prisma method;
        function(s) involved: which_routine
        
    -Check if the array is sorted properly : starting from topmost left point, and sorted clockwise;
        function(s) involved: Modify poly_line, topmost, clockwise
        
    -Return message to user and ask for prompt to sort the array if it's not the case;
        function(s) involved: TBD
            
    -Recenter the height line or prism coordinates around the MPs';
        function(s) involved: normalize
            
    -Calculation of the gravity effect depending on the method choosed;
        function(s) involved for BGPoly: param_val, S_sign, W_sign and check main routine
        function(s) involved for Prisma: calc_T, calc_c, prism, prism_x, prism_y, prism_xy


To do : Function to check if prism are crossing an axis (x or y or both) and use the proper equation accordingly,
while notifying the user.
find a way to input the prism effect calculation, so it is not too lenghty and tedious




"""

orig_stdout = sys.stdout
f = open('out.txt', 'w')
sys.stdout = f
np.seterr(all='raise')
np.seterr(all='print')
start_time = time.time()
##################
# USER INPUTS
#################

#Calling data, need to know how the it is presented in each files
#Measurement points might not be needed if we use the mesh grid

#Import measurement points (MPs) data longitude, latitude, altitude
MPs=np.loadtxt("measurement_data.txt")

#Import Polygons coordinates 
poly_data=x_poly, y_poly, z_poly, rho_poly=np.loadtxt("Sphere_poly.txt", unpack=True)

#Import Prisms coordinates
pri_data=x2s, x1s, y2s, y1s, z2s, z1s, rho_pri=np.loadtxt("prisms_data.txt", unpack=True)

#Create a folder for the runs in the current working directory
run_dir="Test 5"

#Origin of the coordinate system
coor_origin= np.array([600000,200000])

#Gavitational constant in m3 kg-1 s-2 (or N m2 kg-2)
G=6.67408*10**(-11)

#########################################
#FUNCTIONS DEFINITION
#########################################

#Ask the user if they want to use measurement points or a mesh grid for the model
def data_type():
    data_type=input("Please specify data type ""MPs"" or a ""Grid"": ")
    if data_type=="MPs":
        x_mp, y_mp, z_mp=np.loadtxt("measurement_data.txt", unpack=True)
    elif data_type=="Grid":
        x_mp, y_mp, z_mp=input("Specifiy number of points for x y z: ").split()
    return data_type,x_mp, y_mp, z_mp

def which_routine():
    which_routine=input("Which routine, ""BGPoly"" or ""Prisma"" ?")
    return which_routine

def number_height_lines():
    nbhl=input("How many height lines ? An odd number is required")
    return nbhl

#Check if coordinates have n as their height and store them in isoline if it's the case
def poly_line(poly_list,n,isoline):
    #Its purpose is to create an array of (x,y,z) points for a specific height
    #For all the points available search if z is equal to n
    if poly_list[2]==n:
    #in which case store it in isoline
        isoline.append(poly_list)
    #else :
      #  print('the height value select is not present select another or check your values')
        #sys.exit()
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


#Re-arrange height lines points, starting from the topmost left point to facilitate the clokwise function
def topmost(vect):
    #Sort the coordinates from x min to x max regarless of y value
    vect=vect[np.lexsort(vect[:,::-1].T)]
    #Search in isoline if set of coordinates have a same x value
    for i in range(len(vect)-1):
        b=vect[i][1]
        c=vect[i+1][1]
        #The coordinate with the highest y will be the 1st value
        if vect[i][0]==vect[i+1][0] and vect[i][1]<vect[i+1][1] :
           vect[i][1]=c
           vect[i+1][1]=b
    origin=vect[0]
    return origin

#Take an isoline as the entry parameter
#Sort the coordinates by the angles in a clockwise manner
#the origin is the topmost left point of the isoline and is determined for each isoline by a for loop
#the reference vector (which length is 1) is pointing toward the top refv=[0,1]
#this function return the angle AND the length of the vector, to differenciate in case of a same angle
#This function shall be used as the key argument within sorted() to sort the isolines clockwise
def clockwise(isoline_point):
    #Reference vector to sort the isolines coordinates
    refv=[0,1]
    # Vector between point and the origin: v = p - o
    vector = [isoline_point[0]-origin[0], isoline_point[1]-origin[1]]
    # Length of vector: ||v||
    lenvector = math.hypot(vector[0], vector[1])
    # If length is zero there is no angle
    if lenvector == 0:
        return -math.pi, 0
    # Normalize vector: v/||v||
    normalized = [vector[0]/lenvector, vector[1]/lenvector]
    dotprod  = normalized[0]*refv[0] + normalized[1]*refv[1]     # x1*x2 + y1*y2
    diffprod = refv[1]*normalized[0] - refv[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    # Negative angles represent counter-clockwise angles so we need to subtract them 
    # from 2*pi (360 degrees)
    if angle < 0:
        return 2*math.pi+angle, lenvector
    # I return first the angle because that's the primary sorting criterium
    # but if two vectors have the same angle then the shorter distance should come first.
    return angle, lenvector

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
#Function to calculate the gravity effect of height lines
def iso_eff(z_iso, x_iso1, x_iso2, y_iso1, y_iso2, r, r_1, q, p, f, W, S):
    G_E=W*np.arccos((x_iso1/r)*(x_iso2/r_1) + (y_iso1/r)*(y_iso2/r_1)) \
        - np.arcsin(z_iso*q*S/(p*p + z_iso*z_iso)**(1/2)) + np.arcsin(z_iso*f*S/(p*p + z_iso*z_iso)**(1/2))
    return G_E

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

def easy_prism(x,y,z):
#    real_prism=x*np.log(y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) \
#                + y*np.log(x + np.sqrt(np.abs(x)*np.abs(x)+np.abs(y)*np.abs(y)+z*z)) \
#               - z*np.arcsin((z*z+y*y+y*np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) \
#                / ((y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))*np.sqrt(np.abs(y)*np.abs(y)+z*z)))
    
    real_prism1= x*np.log(y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z)) 
    real_prism2= y*np.log(x + np.sqrt(np.abs(x)*np.abs(x)+np.abs(y)*np.abs(y)+z*z)) 
    real_prism3_1=(z*z+y*y+y*np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))
    real_prism3_2=((y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))*np.sqrt(np.abs(y)*np.abs(y)+z*z))
    prism3_diff=np.abs(real_prism3_1-real_prism3_2)
#    real_prism3=z*np.arcsin((z*z+y*y+y*np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))/((y + np.sqrt(np.abs(x)*np.abs(x) + np.abs(y)*np.abs(y)+z*z))*np.sqrt(np.abs(y)*np.abs(y)+z*z)))
    if prism3_diff <= 0.001 :
        real_prism3=z*np.arcsin(1)
    else:
        real_prism3=z*np.arcsin(real_prism3_1/real_prism3_2)
        
    real_prism=real_prism1+real_prism2-real_prism3
    print('')
    print('real_prism1 = '+str(real_prism1),file=file)
    print('real_prism2 = '+str(real_prism2),file=file)
    print('real_prism3 = '+str(real_prism3),file=file)
    print('real_prism3_1 = '+str(real_prism3_1),file=file)
    print('real_prism3_2 = '+str(real_prism3_2),file=file)
    print('',file=file)
    return real_prism

#Will have to be careful, depending on the number of points any operation in the mesh grid could take a lot of time
#Need to limit loop see itertools docmentation
    
#######################################
#MAIN ROUTINE - BGPOLY
#######################################
#Store unique height line values, check that the number of height lines is odd and return to the use
height_lines=[]
iso_rho=[]
for i in range(len(z_poly)):
    if z_poly[i-1] != z_poly[i] :
        height_lines.append(z_poly[i])
        iso_rho.append(rho_poly[i])
height_lines=np.asarray(height_lines)
height_lines.sort() #Makes sure height lines are sorted in an increasing order
iso_rho=np.asarray(iso_rho)

if len(height_lines) % 2 != 0 :
    print('The number of height lines is '+str(len(height_lines)))
else :
    print(' ')
    print('WARNING')
    print(' ')
    print('The number of height lines is even, it has to be odd, the interpolation will not work otherwise.')
    print(' ')
    print('WARNING')
    sys.exit("The amount of height lines isn't odd")
    
os.makedirs(run_dir, exist_ok=True)
V=np.zeros(shape=(len(MPs),len(height_lines)))
V_T=np.zeros(shape=(len(MPs),len(height_lines)))
body_effect=np.zeros(shape=(3,2))
poly_coor=poly_data.T[:,:3] #Store the x,y,z coordinates of the polygons
for i in range(len(MPs)):
    #Announce which point is being processed
    file=open(str(run_dir)+'\BGPoly for MP'+str(i)+'.txt','w')
    
    print(' ',file=file)
    print(' ',file=file)
    print('__________________________________________________',file=file)
    print(' ')
    print('The MP coordinates are '+str(MPs[i]),file=file)
    print('__________________________________________________',file=file)
    
#Create an array to store each height line coordinates
    for j in range(len(height_lines)):
        isoline=[]
        for l in range(len(poly_coor)):
            isoline,poly_list=poly_line(poly_coor[l],height_lines[j],isoline)
        isoline=np.asarray(isoline)
        
    #Store all the coordinates with that height
        True_coor=isoline
        
    #Start to sort the height line array
        origin=topmost(isoline)
        
    #Clockwise
        isoline=np.asarray(sorted(isoline, key=clockwise))  
       ##### if isoline1[0,]!=isoline[0]:
        #####    sys.exit('height lines arrays are not properly sorted')
        
        if isoline[0,0]==isoline[1,0] and isoline[0,1]==isoline[1,1]:
            isoline=np.roll(isoline,-1,axis=0)
            
        if isoline[0,0]!=isoline[len(isoline)-1,0] and isoline[0,1]!=isoline[len(isoline)-1,1]:
            isoline=np.vstack((isoline,isoline[0]))
            
        for L in range(1,len(isoline)-2):
            if np.array_equal(isoline[L],isoline[L+1])==True:
                isoline=np.delete(isoline,[L+1],axis=0)
                
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
        print(' ',file=file) 
        print('The height line '+str(height_lines[j])+'(m) normalized coordinates are: ',file=file) 
        print(' ',file=file)
        print(str(iso_norm),file=file)
        
    #Search for values that are equal to 0 and change their value to 0.1
        for c in range(len(iso_norm)):
            for t in range(len(iso_norm[0])):
                if iso_norm[c,t] is 0:
                    iso_norm[c,t]=0.001
                    print('This value as been replaced'+str(c)+' '+str(t),file=file)
                    
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
        W=np.zeros(len(x_iso))
        S=np.zeros(len(x_iso))
        
        for l in range(len(height_effect[1])) :
    #Calculation of the terms inside of the sum for each points
            W[l]=W_sign(m[l])
            S[l]=S_sign(m[l])
            G_E[l]=iso_eff(z_iso[0], x_iso[l], x_iso[l+1], y_iso[l], y_iso[l+1], r[l], r[l+1], q[l], p[l], f[l],
                   W_sign(m[l]), S_sign(m[l]))
            
        print(' ',file=file)
        print('Result (no unit) of the calculations between each points of height '+str(height_lines[j])+'(m): ',file=file)
        print(' ',file=file)
        print(G_E,file=file)
        
        rho=rho_poly[j]
        
        #Calcul the effect of the height lines on the MP, in s-2
        #Stores the results of each lamina for each MP : V[0,0] is the effect of the 1st lamina for the 1st MP
        V[i,j]=G*rho*sum(G_E)

        print(' ',file=file)
        print('Height line '+str(height_lines[j])+'(m) effect (in s-2) on this MP'+str(MPs[i])+':',file=file)
        print(' ',file=file)
        print(V[i],file=file)
        
    G_E_store=np.zeros(shape=((len(height_lines),len(x_iso)-1)))
    
    for j in range(len(height_lines)):
        G_E_store[j]=G_E
        
        print(' ', file=file)
        print(' ', file=file)
        print('G_E_store', file=file)
        print(' ', file=file)
        print(G_E_store, file=file)
        print(' ', file=file)
        print(' ', file=file)
        
    print(' ', file=file)
    print('End of '+str(MPs[i])+' calculations',file=file)
    
#Calcul of the gravity effect of the whole body
#Depending on the number of height lines, interpolation will be done 3 height lines at a time
z=height_lines

#num_lines makes sures that the numbers of column is the same as the number of interpolations needed, to store the value
#and then sum everything
num_lines=int(len(height_lines)//2)
Body_effect=np.zeros(shape=(len(MPs),num_lines))
MPs_value=np.zeros((len(Body_effect)))
V_sph=np.zeros(len(MPs))  

#Does the integration and interpolation every 3 line
for i in range(len(MPs)):
    for j,l in zip(range(1,len(height_lines),2), range(num_lines)):
        Body_effect[i,l]=body_func(V[i,j-1],V[i,j],V[i,j+1],z[j-1],z[j],z[j+1])
        
#Result here is in Gal
#        print(l,' '+str(MPs[i])+' height calc '+str(j)) #check that the loop is working properly 

for i in range(len(Body_effect)):
    MPs_value[i]=sum(Body_effect[i])*10**3 #Multiplication by 10**3 to convert Gal to mGal
    
    


#######################################
#MAIN ROUTINE - PRISMA
#######################################

#Normalize prism coordinates and set each MP as the origin
#store all the prisms in one array
data_prisms=[x2s,x1s,y2s,y1s,z2s,z1s,rho]
prism_results=np.zeros(shape=(len(MPs),len(x2s)))
prism_results_2=np.zeros(shape=(len(MPs),len(x2s)))
real_results=np.zeros(shape=(len(MPs),len(x2s)))
realest_results=np.zeros(shape=(len(MPs),len(x2s)))
count=0
stor_x=np.zeros(shape=(len(MPs),4))
stor_y=np.zeros(shape=(len(MPs),4))
prism_effect=np.zeros(shape=(len(MPs),len(x2s)))

real_effect=np.zeros(shape=(len(MPs),len(x2s)))
realest_effect=np.zeros(shape=(len(MPs),len(x2s)))


#Select one MP at a time to do the calculation
for i in range(len(MPs)):
#Announce which point is being processed
    file=open(str(run_dir)+'\Prisma for MP'+str(i)+'.txt','w')
    
    print(' ',file=file)
    print(' ',file=file)
    print('__________________________________________________',file=file)
    print(' ')
    print('The MP coordinates are '+str(MPs[i]),file=file)
    print('__________________________________________________',file=file)
    
#Select one prism at a time
#Create an array with the coordinate couples of the selected prism
    for j in range(len(x2s)):
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
        print('Normalized prism(s) coordinates: '+str(coor_norm),file=file)
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
                p1_part1=easy_prism(x2_prism, y2_prism, z1_prism)
                p1_part2=easy_prism(x1_prism, y2_prism, z1_prism)
                p1_part3=easy_prism(x2_prism, y1_prism, z1_prism)
                p1_part4=easy_prism(x1_prism, y1_prism, z1_prism)
                p1_part5=easy_prism(x2_prism, y2_prism, z2_prism)
                p1_part6=easy_prism(x1_prism, y2_prism, z2_prism)
                p1_part7=easy_prism(x2_prism, y1_prism, z2_prism)
                p1_part8=easy_prism(x1_prism, y1_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results[i][j]=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('x and y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
                
            elif y1_prism < 0 and y2_prism < 0:
                p1_part1=easy_prism(x2_prism, -y2_prism, z1_prism)
                p1_part2=easy_prism(x1_prism, -y2_prism, z1_prism)
                p1_part3=easy_prism(x2_prism, -y1_prism, z1_prism)
                p1_part4=easy_prism(x1_prism, -y1_prism, z1_prism)
                p1_part5=easy_prism(x2_prism, -y2_prism, z2_prism)
                p1_part6=easy_prism(x1_prism, -y2_prism, z2_prism)
                p1_part7=easy_prism(x2_prism, -y1_prism, z2_prism)
                p1_part8=easy_prism(x1_prism, -y1_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results[i][j]=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x signs are positive, y signs are negative, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
#For y1 and y2 of different signs and x1 and x2 of same sign
            else :
                p1_part1=easy_prism(x2_prism, y2_prism, z1_prism)
                p1_part2=easy_prism(x1_prism, y2_prism, z1_prism)
                p1_part3=easy_prism(x2_prism, 0, z1_prism)
                p1_part4=easy_prism(x1_prism, 0, z1_prism)
                p1_part5=easy_prism(x2_prism, y2_prism, z2_prism)
                p1_part6=easy_prism(x1_prism, y2_prism, z2_prism)
                p1_part7=easy_prism(x2_prism, 0, z2_prism)
                p1_part8=easy_prism(x1_prism, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results_p1=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('realest results part 1  '+ str(realest_results_p1),file=file)
                print('',file=file)
                
                p2_part1=easy_prism(x2_prism, -y1_prism, z1_prism)
                p2_part2=easy_prism(x1_prism, -y1_prism, z1_prism)
                p2_part3=easy_prism(x2_prism, 0, z1_prism)
                p2_part4=easy_prism(x1_prism, 0, z1_prism)
                p2_part5=easy_prism(x2_prism, -y1_prism, z2_prism)
                p2_part6=easy_prism(x1_prism, -y1_prism, z2_prism)
                p2_part7=easy_prism(x2_prism, 0, z2_prism)
                p2_part8=easy_prism(x1_prism, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p2_1  '+str(p2_part1),file=file)
                print('p2_2  '+str(p2_part2),file=file)
                print('p2_3  '+str(p2_part3),file=file)
                print('p2_4  '+str(p2_part4),file=file)
                print('p2_5  '+str(p2_part5),file=file)
                print('p2_6  '+str(p2_part6),file=file)
                print('p2_7  '+str(p2_part7),file=file)
                print('p2_8  '+str(p2_part8),file=file)
                print('',file=file)
                
                realest_results_p2=p2_part1-p2_part2-p2_part3+p2_part4-p2_part5+p2_part6+p2_part7-p2_part8
                
                print('',file=file)
                print('realest results part 2  '+ str(realest_results_p2),file=file)
                print('',file=file)
                
                realest_results[i][j]=realest_results_p1+realest_results_p2
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x signs are positive, y have different signs, y values are '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
                
        elif x1_prism < 0 and x2_prism < 0:
            if y1_prism >=0 and y2_prism >=0:
                p1_part1=easy_prism(-x1_prism, y2_prism, z1_prism)
                p1_part2=easy_prism(-x2_prism, y2_prism, z1_prism)
                p1_part3=easy_prism(-x1_prism, y1_prism, z1_prism)
                p1_part4=easy_prism(-x2_prism, y1_prism, z1_prism)
                p1_part5=easy_prism(-x1_prism, y2_prism, z2_prism)
                p1_part6=easy_prism(-x2_prism, y2_prism, z2_prism)
                p1_part7=easy_prism(-x1_prism, y1_prism, z2_prism)
                p1_part8=easy_prism(-x2_prism, y1_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results[i][j]=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x signs are negative, y signs are positive, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
                
            elif y1_prism < 0 and y2_prism < 0:
                p1_part1=easy_prism(-x1_prism, -y1_prism, z1_prism)
                p1_part2=easy_prism(-x2_prism, -y1_prism, z1_prism)
                p1_part3=easy_prism(-x1_prism, -y2_prism, z1_prism)
                p1_part4=easy_prism(-x2_prism, -y2_prism, z1_prism)
                p1_part5=easy_prism(-x1_prism, -y1_prism, z2_prism)
                p1_part6=easy_prism(-x2_prism, -y1_prism, z2_prism)
                p1_part7=easy_prism(-x1_prism, -y2_prism, z2_prism)
                p1_part8=easy_prism(-x2_prism, -y2_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results[i][j]=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x and y signs are negative, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
#For y1 and y2 of different signs and x1 and x2 of same sign
            else :
                p1_part1=easy_prism(-x1_prism, y2_prism, z1_prism)
                p1_part2=easy_prism(-x2_prism, y2_prism, z1_prism)
                p1_part3=easy_prism(-x1_prism, 0, z1_prism)
                p1_part4=easy_prism(-x2_prism, 0, z1_prism)
                p1_part5=easy_prism(-x1_prism, y2_prism, z2_prism)
                p1_part6=easy_prism(-x2_prism, y2_prism, z2_prism)
                p1_part7=easy_prism(-x1_prism, 0, z2_prism)
                p1_part8=easy_prism(-x2_prism, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results_p1=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('realest results part 1  '+ str(realest_results_p1),file=file)
                print('',file=file)
                
                p2_part1=easy_prism(-x1_prism, -y1_prism, z1_prism)
                p2_part2=easy_prism(-x2_prism, -y1_prism, z1_prism)
                p2_part3=easy_prism(-x1_prism, 0, z1_prism)
                p2_part4=easy_prism(-x2_prism, 0, z1_prism)
                p2_part5=easy_prism(-x1_prism, -y1_prism, z2_prism)
                p2_part6=easy_prism(-x2_prism, -y1_prism, z2_prism)
                p2_part7=easy_prism(-x1_prism, 0, z2_prism)
                p2_part8=easy_prism(-x2_prism, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p2_1  '+str(p2_part1),file=file)
                print('p2_2  '+str(p2_part2),file=file)
                print('p2_3  '+str(p2_part3),file=file)
                print('p2_4  '+str(p2_part4),file=file)
                print('p2_5  '+str(p2_part5),file=file)
                print('p2_6  '+str(p2_part6),file=file)
                print('p2_7  '+str(p2_part7),file=file)
                print('p2_8  '+str(p2_part8),file=file)
                print('',file=file)
                
                realest_results_p2=p2_part1-p2_part2-p2_part3+p2_part4-p2_part5+p2_part6+p2_part7-p2_part8
                
                print('',file=file)
                print('realest results part 2  '+ str(realest_results_p2),file=file)
                print('',file=file)
                
                realest_results[i][j]=realest_results_p1+realest_results_p2
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x signs are negative, y have different signs, y values are '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
                
#For x1 and x2 of different signs and y1 and y2 of same signs
        else :
            if y1_prism >=0 and y2_prism >=0:
                p1_part1=easy_prism(x2_prism, y2_prism, z1_prism)
                p1_part2=easy_prism(0, y2_prism, z1_prism)
                p1_part3=easy_prism(x2_prism, y1_prism, z1_prism)
                p1_part4=easy_prism(0, y1_prism, z1_prism)
                p1_part5=easy_prism(x2_prism, y2_prism, z2_prism)
                p1_part6=easy_prism(0, y2_prism, z2_prism)
                p1_part7=easy_prism(x2_prism, y1_prism, z2_prism)
                p1_part8=easy_prism(0, y1_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results_p1=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('realest results part 1  '+ str(realest_results_p1),file=file)
                print('',file=file)
                
                p2_part1=easy_prism(-x1_prism, y2_prism, z1_prism)
                p2_part2=easy_prism(0, y2_prism, z1_prism)
                p2_part3=easy_prism(-x1_prism, y1_prism, z1_prism)
                p2_part4=easy_prism(0, y1_prism, z1_prism)
                p2_part5=easy_prism(-x1_prism, y2_prism, z2_prism)
                p2_part6=easy_prism(0, y2_prism, z2_prism)
                p2_part7=easy_prism(-x1_prism, y1_prism, z2_prism)
                p2_part8=easy_prism(0, y1_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p2_1  '+str(p2_part1),file=file)
                print('p2_2  '+str(p2_part2),file=file)
                print('p2_3  '+str(p2_part3),file=file)
                print('p2_4  '+str(p2_part4),file=file)
                print('p2_5  '+str(p2_part5),file=file)
                print('p2_6  '+str(p2_part6),file=file)
                print('p2_7  '+str(p2_part7),file=file)
                print('p2_8  '+str(p2_part8),file=file)
                print('',file=file)
                
                realest_results_p2=p2_part1-p2_part2-p2_part3+p2_part4-p2_part5+p2_part6+p2_part7-p2_part8
                
                print('',file=file)
                print('realest results part 2  '+ str(realest_results_p2),file=file)
                print('',file=file)
                
                realest_results[i][j]=realest_results_p1+realest_results_p2
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x have different signs, y signs are positive, x values are '+str(x2_prism)+' '+str(x1_prism),file=file)
                print('',file=file)
                
            elif y1_prism < 0 and y2_prism < 0:
                p1_part1=easy_prism(x2_prism, -y1_prism, z1_prism)
                p1_part2=easy_prism(0, -y1_prism, z1_prism)
                p1_part3=easy_prism(x2_prism, -y2_prism, z1_prism)
                p1_part4=easy_prism(0, -y2_prism, z1_prism)
                p1_part5=easy_prism(x2_prism, -y1_prism, z2_prism)
                p1_part6=easy_prism(0, -y1_prism, z2_prism)
                p1_part7=easy_prism(x2_prism, -y2_prism, z2_prism)
                p1_part8=easy_prism(0, -y2_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results_p1=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('realest results part 1  '+ str(realest_results_p1),file=file)
                print('',file=file)
                
                p2_part1=easy_prism(-x1_prism, -y1_prism, z1_prism)
                p2_part2=easy_prism(0, -y1_prism, z1_prism)
                p2_part3=easy_prism(-x1_prism, -y2_prism, z1_prism)
                p2_part4=easy_prism(0, -y2_prism, z1_prism)
                p2_part5=easy_prism(-x1_prism, -y1_prism, z2_prism)
                p2_part6=easy_prism(0, -y1_prism, z2_prism)
                p2_part7=easy_prism(-x1_prism, -y2_prism, z2_prism)
                p2_part8=easy_prism(0, -y2_prism, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p2_1  '+str(p2_part1),file=file)
                print('p2_2  '+str(p2_part2),file=file)
                print('p2_3  '+str(p2_part3),file=file)
                print('p2_4  '+str(p2_part4),file=file)
                print('p2_5  '+str(p2_part5),file=file)
                print('p2_6  '+str(p2_part6),file=file)
                print('p2_7  '+str(p2_part7),file=file)
                print('p2_8  '+str(p2_part8),file=file)
                print('',file=file)
                
                realest_results_p2=p2_part1-p2_part2-p2_part3+p2_part4-p2_part5+p2_part6+p2_part7-p2_part8
                
                print('',file=file)
                print('realest results part 2  '+ str(realest_results_p2),file=file)
                print('',file=file)
                
                realest_results[i][j]=realest_results_p1+realest_results_p2
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                print('',file=file)
                print('x have different signs, y signs are negative, x values are '+str(x2_prism)+' '+str(x1_prism),file=file)
                print('',file=file)
#For x1 and x2, y1 and y2 of different signs
            else :
                p1_part1=easy_prism(x2_prism, y2_prism, z1_prism)
                p1_part2=easy_prism(0, y2_prism, z1_prism)
                p1_part3=easy_prism(x2_prism, 0, z1_prism)
                p1_part4=easy_prism(0, 0, z1_prism)
                p1_part5=easy_prism(x2_prism, y2_prism, z2_prism)
                p1_part6=easy_prism(0, y2_prism, z2_prism)
                p1_part7=easy_prism(x2_prism, 0, z2_prism)
                p1_part8=easy_prism(0, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p1_1  '+str(p1_part1),file=file)
                print('p1_2  '+str(p1_part2),file=file)
                print('p1_3  '+str(p1_part3),file=file)
                print('p1_4  '+str(p1_part4),file=file)
                print('p1_5  '+str(p1_part5),file=file)
                print('p1_6  '+str(p1_part6),file=file)
                print('p1_7  '+str(p1_part7),file=file)
                print('p1_8  '+str(p1_part8),file=file)
                print('',file=file)
                
                realest_results_p1=p1_part1-p1_part2-p1_part3+p1_part4-p1_part5+p1_part6+p1_part7-p1_part8
                
                print('',file=file)
                print('realest results part 1  '+ str(realest_results_p1),file=file)
                print('',file=file)
                
                p2_part1=easy_prism(-x1_prism, y2_prism, z1_prism)
                p2_part2=easy_prism(0, y2_prism, z1_prism)
                p2_part3=easy_prism(-x1_prism, 0, z1_prism)
                p2_part4=easy_prism(0, 0, z1_prism)
                p2_part5=easy_prism(-x1_prism, y2_prism, z2_prism)
                p2_part6=easy_prism(0, y2_prism, z2_prism)
                p2_part7=easy_prism(-x1_prism, 0, z2_prism)
                p2_part8=easy_prism(0, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p2_1  '+str(p2_part1),file=file)
                print('p2_2  '+str(p2_part2),file=file)
                print('p2_3  '+str(p2_part3),file=file)
                print('p2_4  '+str(p2_part4),file=file)
                print('p2_5  '+str(p2_part5),file=file)
                print('p2_6  '+str(p2_part6),file=file)
                print('p2_7  '+str(p2_part7),file=file)
                print('p2_8  '+str(p2_part8),file=file)
                print('',file=file)
                
                realest_results_p2=p2_part1-p2_part2-p2_part3+p2_part4-p2_part5+p2_part6+p2_part7-p2_part8
                
                print('',file=file)
                print('realest results part 2  '+ str(realest_results_p2),file=file)
                print('',file=file)
                
                p3_part1=easy_prism(x2_prism, -y1_prism, z1_prism)
                p3_part2=easy_prism(0, -y1_prism, z1_prism)
                p3_part3=easy_prism(x2_prism, 0, z1_prism)
                p3_part4=easy_prism(0, 0, z1_prism)
                p3_part5=easy_prism(x2_prism, -y1_prism, z2_prism)
                p3_part6=easy_prism(0, -y1_prism, z2_prism)
                p3_part7=easy_prism(x2_prism, 0, z2_prism)
                p3_part8=easy_prism(0, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p3_1  '+str(p3_part1),file=file)
                print('p3_2  '+str(p3_part2),file=file)
                print('p3_3  '+str(p3_part3),file=file)
                print('p3_4  '+str(p3_part4),file=file)
                print('p3_5  '+str(p3_part5),file=file)
                print('p3_6  '+str(p3_part6),file=file)
                print('p3_7  '+str(p3_part7),file=file)
                print('p3_8  '+str(p3_part8),file=file)
                print('',file=file)
                
                realest_results_p3=p3_part1-p3_part2-p3_part3+p3_part4-p3_part5+p3_part6+p3_part7-p3_part8
                
                print('',file=file)
                print('realest results part 3  '+ str(realest_results_p3),file=file)
                print('',file=file)
                
                p4_part1=easy_prism(-x1_prism, -y1_prism, z1_prism)
                p4_part2=easy_prism(0, -y1_prism, z1_prism)
                p4_part3=easy_prism(-x1_prism, 0, z1_prism)
                p4_part4=easy_prism(0, 0, z1_prism)
                p4_part5=easy_prism(-x1_prism, -y1_prism, z2_prism)
                p4_part6=easy_prism(0, -y1_prism, z2_prism)
                p4_part7=easy_prism(-x1_prism, 0, z2_prism)
                p4_part8=easy_prism(0, 0, z2_prism)
                
                print('',file=file)
                print('Easy_prism results :',file=file)
                print('',file=file)
                print('p4_1  '+str(p4_part1),file=file)
                print('p4_2  '+str(p4_part2),file=file)
                print('p4_3  '+str(p4_part3),file=file)
                print('p4_4  '+str(p4_part4),file=file)
                print('p4_5  '+str(p4_part5),file=file)
                print('p4_6  '+str(p4_part6),file=file)
                print('p4_7  '+str(p4_part7),file=file)
                print('p4_8  '+str(p4_part8),file=file)
                print('',file=file)
                
                realest_results_p4=p4_part1-p4_part2-p4_part3+p4_part4-p4_part5+p4_part6+p4_part7-p4_part8
                
                print('',file=file)
                print('realest results part 4  '+ str(realest_results_p4),file=file)
                print('',file=file)
                
                realest_results[i][j]=realest_results_p1+realest_results_p2+realest_results_p3+realest_results_p4
                
                print('',file=file)
                print('Realest results : '+str(realest_results[i][j]),file=file)
                print('',file=file)
                
                print('',file=file)
                print('x and y have different signs, x and y values are '+str(x2_prism)+' '+str(x1_prism)+' and '+str(y2_prism)+' '+str(y1_prism),file=file)
                print('',file=file)
                
#Calculate the gravity effect for each points for the large and small prism
    print(' ')
for i,j in itertools.product(range(len(x2s)),range(len(MPs))):
    realest_effect[j][i]=realest_results[:][j][i]*rho_pri[i]*G*10**3
with open(str(run_dir)+' Prisma 1 results.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(realest_effect[:,0]))
with open(str(run_dir)+' Prisma 2 results.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(realest_effect[:,1]))
with open(str(run_dir)+' Prisma 3 results.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(realest_effect[:,2]))
with open(str(run_dir)+' Prisma 4 results.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(realest_effect[:,3]))
with open(str(run_dir)+' Prisma 5 results.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(realest_effect[:,4]))
with open(str(run_dir)+' Prisma 6 results.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(realest_effect[:,5]))
    
print("--- %s seconds ---" % (time.time() - start_time))
