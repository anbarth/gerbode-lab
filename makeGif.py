import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)
import imageio

dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\peanut_experiment\full peanut'
images = []
# make all frames
for filename in os.listdir(dir):   
    if filename[-3:]=='png':
        images.append(imageio.imread("../peanut_experiment/full peanut/"+filename))
# make the gif!
imageio.mimsave('waaa.gif', images)