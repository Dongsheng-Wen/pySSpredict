
import ternary
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc, font_manager
ticks_font = font_manager.FontProperties(family='serif', style='normal',
                                         size=24, weight='normal', stretch='normal')
# strength/property plotter -2d 
class plotter_2d:

    def __init__(self,
                 data,
                 x_name,
                 y_name):

        self.X_col = x_name
        self.Y_col = y_name
        self.X = data[self.X_col]
        self.Y = data[self.Y_col]

    def plot(self):
        



# strength/property plotter -3d
class plotter_3d:

    def __init__(self,X,Y,Z,property):

        self.X = 
        self.Y = 
        self.Z = 
        self.H = property



    def plot(self):

        self.fig, self.tax = ternary