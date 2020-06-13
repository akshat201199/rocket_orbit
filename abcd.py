# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 12:03:30 2020

@author: Akshat Shah
"""

import orbit_plotter as orpl
import Orbits2 as orpl2
from math import *
import numpy as np

#orpl.create_orbit([0,6.3781e+06], [0,0], np.linspace(0,200,1000), 550000,2770,1935,[0.005,1],True)
#orpl.plot_orbit()

t_vector = np.linspace(0,1500,1000)
orientation_vect = np.full((1000,2),[0,1])
eng_stat = np.full((1000,1),1)


orpl2.create_orbit_fs([0,6.3781e+06], [465,0], t_vector, 395000,2770,1535,orientation_vect,eng_stat,False,36537,291)
orpl2.plot_orbit()

