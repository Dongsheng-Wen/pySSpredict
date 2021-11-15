
import ternary
import matplotlib.pyplot as plt
import matplotlib as mpl 
import numpy as np
import pandas as pd
from matplotlib import rc, font_manager

# strength/property plotter -2d 
class plotter_2d:

    def __init__(self,
                 data,
                 X='T',
                 Y='tau_y'):
        self.ticks_font = font_manager.FontProperties(family='serif', style='normal',
                                         size=24, weight='normal', stretch='normal')
        self.X = data[self.X_col]
        self.Y = data[self.Y_col]
        plt.draw()
        self.fig = plt.figure(figsize=(12,9))
        self.ax = self.f.subplots(1,1)

        self.ax.plot(self.X,self.Y,'b-')
        self.ax.tick_params(direction='in',axis='both',length=10,right=True)
        for labelx1 in (ax.get_xticklabels()):
            labelx1.set_fontproperties(self.ticks_font)
    
        for labely1 in (ax.get_yticklabels()):
            labely1.set_fontproperties(self.ticks_font)

        self.ax.set_xlabel(X,fontproperties=self.ticks_font)
        self.ax.set_ylabel(Y,fontproperties=self.ticks_font)
        self.ax.legend(fontsize=24,frameon=False)
        

    def show(self):
        
        plt.show()


# strength/property plotter -3d
class pyss_plotter_3d:

    def __init__(self,data,X='comp(psA)',Y='comp(psB)',Z='comp(psC)',H='tau_y',scale=100):

        self.ticks_font = font_manager.FontProperties(family='serif', style='normal',
                                         size=24, weight='normal', stretch='normal')
        self.data = data # pd.dataframe
        self.X_name = X 
        self.Y_name = Y 
        self.Z_name = Z
        self.H_name = H 
        self.X = np.round(self.data[X],3) # 1d
        self.Y = np.round(self.data[Y],3) # 1d 
        self.Z = np.round(self.data[Z],3) # 1d 
        self.H = np.round(self.data[H],3) # 1d
        self.points = list(zip(self.X,self.Y,self.Z))
        self.hm_dict = dict(zip(self.points,self.H))
        self.scale = scale 
        self.vmin = min(self.H)
        self.vmax = max(self.H)
                # default plot settings
        self.fig, self.tax = ternary.figure(scale=self.scale)
        self.fig.set_figwidth(16)
        self.fig.set_figheight(12)
        self.fig.set_dpi(300)

        self.tax.get_axes().axis('off')
        self.tax.clear_matplotlib_ticks()
        self.tax.resize_drawing_canvas()

        self.tax.boundary(linewidth=2.0)
        self.tax.gridlines(color="blue", multiple=5)
        self.tax.ax.axis("equal")
        self.cmap='Wistia'

        fs = 24
        left_corner = self.Z_name
        right_corner = self.X_name
        top_corner = self.Y_name 
        self.tax.right_corner_label(right_corner, fontsize=fs,offset=0.02) 
        self.tax.top_corner_label(top_corner, fontsize=fs,offset=0.15)
        self.tax.left_corner_label(left_corner, fontsize=fs,offset=0.08)
        left_axis_label = self.Z_name
        right_axis_label = self.Y_name 
        bottom_axis_label = self.X_name
        self.tax.left_axis_label(left_axis_label, fontsize=fs, offset=0.16)
        self.tax.right_axis_label(right_axis_label, fontsize=fs, offset=0.16)
        self.tax.bottom_axis_label(bottom_axis_label, fontsize=fs, offset=0.16)
        self.tax.ticks(axis='lbr', multiple=10, linewidth=1, offset=0.025,fontsize=fs)


        self.tax.heatmap(self.hm_dict,cmap=self.cmap,
            colorbar=False,vmin=self.vmin,vmax=self.vmax)
        self.h = self.tax.get_figure()

        self.norm = mpl.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.cbar = self.h.colorbar(mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap))
        self.cbar.set_label(label=self.H_name,size=fs)
        self.cbar.ax.tick_params(labelsize=fs)
        plt.close()

        

    def set_corner_labels(self,left_corner='Z',right_corner='X',top_corner='Y',fontsize=24):

        self.tax.right_corner_label(right_corner, fontsize=fontsize,offset=0.02) 
        self.tax.top_corner_label(top_corner, fontsize=fontsize,offset=0.15)
        self.tax.left_corner_label(left_corner, fontsize=fontsize,offset=0.08)

    def set_axes_labels(self,left_axis_label='Z',right_axis_label='Y',bottom_axis_label='X',fontsize=24,offset=0.16):
        self.tax.left_axis_label(left_axis_label, fontsize=fontsize, offset=offset)
        self.tax.right_axis_label(right_axis_label, fontsize=fontsize, offset=offset)
        self.tax.bottom_axis_label(bottom_axis_label, fontsize=fontsize, offset=offset)

    def set_ticks(self,axis='lbr', multiple=10, linewidth=1, offset=0.025,fontsize=24):
    
        self.tax.ticks(axis=axis, multiple=multiple, linewidth=linewidth, 
                        offset=offset,fontsize=fontsize)

    def show(self):
        return self.tax.get_figure()


