import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)
import imageio

dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\grainsplits\r70_103\snowflakes'
images = []
# make all frames
for filename in os.listdir(dir):   
    if filename[-3:]=='png':
        images.append(imageio.imread("../grainsplits/r70_103/snowflakes/"+filename))
# make the gif!
imageio.mimsave('myfakeevent.gif', images)