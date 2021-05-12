#!/usr/bin/python
'''
 You can use this script with python 3(.6.x) to plot your atomic structure and its bonds within your preferred range.
 Tip #0: In order to remove C-H bonds from your plot, choose (vmin > 1.1), also H atoms are to be plotted with a white marker by default.
 It's only been tested on smaller/aromatic systems, so be careful and in case of running into trouble please send me an email. Go science!
 Shayan Edalatmanesh (shyn.ein@gmail.com)
'''

import os
import numpy as np
from pylab import genfromtxt
import math
import matplotlib as mpl
# mpl.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
# import cmocean   # a set of standard colormaps - needs to be installed separately
#===========================================================================================================
# This part is defined by user
# Assuming structure is periodic along the X (/Y) axis:

fname = 'start.xyz'  # The name of (and the path to) your file
vmin = 1.30      # color bar min
vmax = 1.50      # color bar max
bondWeight = 6   # Thickness of the bonds:
showLengths = True    # Bond length values to be printed on top of each bond? (True/False)
blFont = 12      # fontsize
#===========================================================================================================

def findBonds(rmin, rmax):
    rm = float(rmin)
    rr = float(rmax)
    bonds = []
    radiz = []
    n = len(xs)
    for i in np.arange(n):
        for j in np.arange(i):
            dx = xs[j] - xs[i]
            dy = ys[j] - ys[i]
            dz = zs[j] - zs[i]
            r = math.sqrt(dx * dx + dy * dy + dz * dz)
            if (r < rr) and (rm < r):
                bonds.append((i, j))
                radiz.append(r)
            # plt.arrow(xs[i], ys[i], xs[j]-xs[i], ys[j]-ys[i], head_width=0.0, head_length=0.0,normalize_data=scaled, lw= 1.0,ls='solid',zorder=2 )
    return bonds, radiz


def plotBonds(bb):
    for (q1, q2) in bb:
        dx = xs[q2] - xs[q1]
        dy = ys[q2] - ys[q1]
        dz = zs[q2] - zs[q1]
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        rnd = round(r, 2)
        # fig = plt.figure()
        # ax  = fig.add_axes([0.1, 0.1, 0.7, 0.85]) # [left, bottom, width, height]
        # axc = fig.add_axes([0.85, 0.10, 0.05, 0.85])
        cmap = plt.cm.jet
        ### cmap = cmocean.cm.balance
        #  cmap = cmocean.cm.cmap_d
        cNorm = colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
        colorVal = scalarMap.to_rgba(r)
        # cbar = scalarMap.set_array(r)
        # cb1 = mpl.colorbar.ColorbarBase(axc, cmap=cmap, norm=cNorm,orientation='vertical')
        plt.arrow(xs[q1], ys[q1], xs[q2] - xs[q1], ys[q2] - ys[q1], head_width=0.0, head_length=0.0, fc='none',
                  color=colorVal, lw=bondWeight, ls='-', alpha=0.99)  # ,zorder=2
        if (showLengths == True):
            plt.text((xs[q1] + dx / 2) - 0.4, (ys[q1] + dy / 2) + 0.1, str(rnd), fontsize=blFont, fontweight='bold',
                     bbox=props)
        # plt.colorbar(im)
        # plt.text()
        # mpl.colorbar.ColorbarBase(cmap=cmap, norm=cNorm,orientation='vertical')
    #  print (q1,q2)


def pltcolor(lst):
    cols = []
    for l in lst:
        if l == 'C':
            cols.append('dimgrey')
        elif l == 'H':
            cols.append('w')
        elif l == 'Au':
            cols.append('goldenrod')
        else:
            cols.append('g')
    return cols

def conversion(num):
    try:
        num = float(num)
        return True
    except Exception:
        return False


#===========================================================================================================
# Flow starts here:

fil0 = genfromtxt(fname, skip_header=2)
data,target = np.array_split(np.loadtxt(fname, skiprows=2, dtype=str), [-1], axis=1)
xs = fil0[:, 1]
ys = fil0[:, 2]
zs = fil0[:, 3]

print("X: min & max: ", min(xs), max(xs))
print("Y: min & max: ", min(ys), max(ys))
print("XYZ file read, now proceeding with the calculations.")

# plotting window (coordinates are in Angstroms, if you wanna do it manually)
xStart = min(xs) - 1.5
xStop = max(xs) + 1.5
yStart = min(ys) - 1.5
yStop = max(ys) + 1.5

deltaX = max(xs) - min(xs)
deltaY = max(ys) - min(ys)

bonds, radiz = findBonds(vmin, vmax)
# print("DEBUG: Bonds are between atoms: ", bonds, "\n")
print("Total no. of bonds: ", len(radiz))
# radiz2 = sorted(radiz)
# print("DEBUG: Sorted bond lengths (Ang): ", radiz2, "\n")

print("Average bond length: ", sum(radiz) / len(radiz), " (Ang)")

numbers = [float(num) for num in radiz if conversion(num)]
minim = 99999999999999
maximum = -99999999999999
if numbers:
    if max(numbers) > maximum:
        maximum = max(numbers)
    if min(numbers) < minim:
        minim = min(numbers)

print("The shortest bond: " + str(minim))
print("The longest bond: " + str(maximum))
# print('DEBUG: Delta d:', maximum - minim)

smallRadiz = [x for x in radiz if x <= 1.3]
if len(smallRadiz) == 0:
    print("No triple bonds were found!")
else:
    print('I have found', len(smallRadiz), 'triple bonds.')
    print('The triple bonds are:', smallRadiz)
    print("Average length of the triple bonds (Ang) = ", sum(smallRadiz) / len(smallRadiz))

#===========================================================================================================
# Plotting section:

print(" Working on output files, please wait... ")
fname2 = os.path.splitext(fname)[0]
fig = plt.figure(figsize=(deltaX / 2, deltaY / 2))
plt.subplot(1, 1, 1)
# plt.title("C-C Bonds: "+str(fname2),fontsize = 20)
plt.xlim(xStart, xStop)
plt.ylim(yStart, yStop)
plt.xticks([])
plt.yticks([])
cols = pltcolor(data[:, 0])
plt.scatter(fil0[:, 1], fil0[:, 2], zorder=1, marker='o', alpha=1.0, facecolors=cols, edgecolors=cols)
props = dict(boxstyle='round', facecolor='none', edgecolor='black', alpha=0.01)
plotBonds(bonds)
ax = fig.add_axes([0.92, .12, .02, .75])  # position of the colorbar [left, bottom, width, height]
cmap = plt.cm.jet
# cmap = cmocean.cm.balance # standard colors - need to install and import cmocean to use it.
cNorm = colors.Normalize(vmin=vmin, vmax=vmax)
cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=cNorm)  # set min, max of colorbar
cbar.ax.tick_params(labelsize=25)
#===========================================================================================================
# Saving the output:
print("Saving the output image as: ", fname2 + "_.png")
plt.savefig(fname2 + "_.png", format='png', dpi=600)  # ,bbox_inches='tight'
plt.show()
plt.close()
print("ALL DONE!")
