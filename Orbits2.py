# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 11:38:12 2020

@author: Akshat Shah
"""

from scipy.integrate import ode,odeint
import numpy as np
import matplotlib.pyplot as plt
from math import *
import time
from matplotlib.animation import FuncAnimation

G = 6.67408e-11
Me = 5.9722e+24
Re = 6.3781e+06

'''This function creates a first stage orbit in 2D using given parameters. 
r_init - initial position vector of rocket; on earth's surface
v_init - initial velocity vector of rocket; usually 0 
init_mass - initial total mass of rocket
t_vector - time vector of length N at which orbit of rocket is calculated
v_exh - exhaust velocity of rocket
dmdt - dm/dt of rocket i.e amount of fuel used per second
orientation_vec - a matrix of length N,2 containing orientation vector at each time step
engine_status - a vector of length N with entries True(engine On at that time) or False(engine Off)
second_stage - True (second satge is present) or False (absent)
second_stage_mass - Mass of second stage
seperation time - time at which first stage seperates from second stage
'''

def create_orbit_fs(r_init, v_init, t_vector, init_mass, v_exh, dmdt, orientation_vec, engine_status,second_stage,second_stage_mass,seperation_time):
    time_step = np.diff(t_vector)[0]
    length = np.shape(t_vector)[0]
    ori_length = np.shape(orientation_vec)[0]
    eng_length = np.shape(engine_status)[0]
    r_i = r_init
    v_i = v_init
    mass_i = init_mass
    rx =[]
    ry=[]
    vx=[]
    vy=[]
    estat=[]
    if second_stage == False:
        second_stage_mass = 0
    else:
        pass
    fuel_meter = 0.9*(mass_i-second_stage_mass) '''90% of first stage mass is fuel'''
    if length == ori_length and length == eng_length:
        for i in range(0,length-1):
            t_vec = [t_vector[i],t_vector[i+1]]
            ori_vec = [orientation_vec[i]][0]
            eng = [engine_status[i]][0]
            if t_vector[i]  > seperation_time:
                mass_i = mass_i - second_stage_mass
            else:
                pass
            def gravity_eqn(r, t):
                mod_re=np.sqrt((r[0]**2)+(r[1]**2))
                mod_t = np.sqrt((ori_vec[0]**2)+(ori_vec[1]**2))
                if eng == 1 and fuel_meter > 0 :
                    dr0=r[2]
                    dr1=r[3]
                    dr2=-(G*Me*r[0]/(mod_re**3)) + eng*v_exh*dmdt*ori_vec[0]/(mod_t*(fuel_meter+0.1*mass_i) 
                    dr3=-(G*Me*r[1]/(mod_re**3)) + eng*v_exh*dmdt*ori_vec[1]/(mod_t*(fuel_meter+0.1*mass_i))
                else:
                    dr0=r[2]
                    dr1=r[3]
                    dr2=-(G*Me*r[0]/(mod_re**3)) 
                    dr3=-(G*Me*r[1]/(mod_re**3))
                return [dr0,dr1,dr2,dr3]    
            r_solve = odeint(gravity_eqn, [r_i[0],r_i[1],v_i[0],v_i[1]],t_vec)
            rx.append(r_solve[1][0])
            ry.append(r_solve[1][1])
            vx.append(r_solve[1][2])
            vy.append(r_solve[1][3])
            if eng ==1 and fuel_meter > 0 :
                estat.append('On')
            else:
                estat.append('Off')
            r_i = r_solve[1][0:2]
            v_i = r_solve[1][2:4]
            fuel_meter = fuel_meter - eng*dmdt*time_step
    else: 
        print('Error')
    global rx_new
    rx_new=[]
    global ry_new
    ry_new=[]
    global vx_new
    vx_new=[] 
    global vy_new  
    vy_new=[]
    global t_new
    t_new=[]
    t_new_int=[]
    global estat_new
    estat_new=[]
    for x,y,v_x,v_y,e,t in zip(rx,ry,vx,vy,estat,t_vector[0:length]):
        if ((x**2) + (y**2)) > (Re-0.001)**2 :
            rx_new.append(x)
            ry_new.append(y)
            vx_new.append(v_x)
            vy_new.append(v_y)
            estat_new.append(e)
            t_new.append(t)
            t_new_int.append(int(t))
    global t_sep
    t_sep = seperation_time
    t_sep_index = t_new_int.index(t_sep)
    global r_ss_init
    r_ss_init = (rx_new[t_sep_index],ry_new[t_sep_index]) 
    global v_ss_init
    v_ss_init = (vx_new[t_sep_index],ry_new[t_sep_index])  
    global s_mass
    s_mass = second_stage_mass
                       
def plot_orbit():        
    xdata=[]
    ydata=[] 
    fig,ax = plt.subplots()
    earth = plt.Circle((0, 0), Re, color='blue')
    def animate(i):
        try:
            ax.cla()
            xdata.append(rx_new[i])
            ydata.append(ry_new[i])
            time = int(t_new[i])
            velocity = int(np.sqrt((vx_new[i]**2)+(vy_new[i]**2)))
            altitude = int((-1*Re + np.sqrt((rx_new[i]**2)+(ry_new[i]**2)))/1000)
            stat_e = estat_new[i]
            landing = (rx_new[-1],ry_new[-1])
            ax.plot(xdata,ydata,color='k')
            ax.plot(landing[0],landing[1],'ro')
            ax.add_artist(earth)
            text1='Time Elapsed = {a} s \n Velocity = {b} m/s \n Altitude = {c} km \n Engine Status = {d} \n Estimated Landing at {e}'.format(a=time,b=velocity,c=altitude,d=stat_e,e=landing)
            ax.set_title(text1 ,loc='center')
            ax.set_xlim([-0.1*Re,0.1*Re])
            ax.set_ylim([0.99*Re,1.7*Re])
        except IndexError:
            pass
            
    ani=FuncAnimation(fig,animate,interval=1)            
    
    plt.show()


            