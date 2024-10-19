# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 13:23:12 2024

@authors: Very fun people
"""

import numpy as np
import imageio as iio
import matplotlib.pyplot as plt
from skimage.transform import radon

### Physical Parameters 
N = 50 # Number of equally spaced nails on the base
R = 300 # Radius of circular base in mm; centre to nails
stringWidth = 1 # String thickness in mm

### Input Image Path
ImgPath = "../TestImages/BlackCircle.png"

###################################################################
#%% ################## Radon Transform Image ######################
###################################################################
# Adapted from https://scikit-image.org/docs/stable
        # /auto_examples/transform/plot_radon_transform.html

image = np.array(iio.imread(ImgPath)); # Read image from file

theta = np.linspace(0.0, 180.0, N, endpoint=False) # Number of angles considered for the radon transform
sinogram = radon(image, theta=theta) # Radon transform

# Start of showing image and radon transform side by side
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
ax1.set_title("Original")
ax1.imshow(image, cmap=plt.cm.Greys_r)
dx, dy = 0.5 * 180.0 / N, 0.5 / sinogram.shape[0]
ax2.set_title("Radon transform\n(Sinogram)")
ax2.set_xlabel("Projection angle (deg)")
ax2.set_ylabel("Projection position (pixels)")
ax2.imshow(
    sinogram,
    cmap=plt.cm.Greys_r,
    extent=(-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),
    aspect='auto',
)
fig1.tight_layout()
plt.show()
# End of showing image and radon transform side by side

###################################################################
#%% ################## Get Images of All Lines #####################
###################################################################

nailAngles = np.linspace(0, 2*np.pi, num=N, endpoint=False) # Polar angle coordinates of the nails
nailCoors = R * np.array([np.cos(nailAngles), np.sin(nailAngles)]).T

# Start of plotting nail positions
fig2 = plt.figure()
plt.scatter(nailCoors[:,0], nailCoors[:,1])
ax3 = plt.gca()
ax3.set_aspect('equal')
ax3.set_title('Nail coordinates in world measurements')
plt.show()
# End of plotting nail positions

pixelSize = 2*R / max(image.shape) # Size of a pixel in mm
pixelCoorsX = np.kron(np.ones([image.shape[0],1]), np.array(range(image.shape[0]))) +0.5 # X coordinates of pixel centres in pixels
pixelCoorsY = np.kron(np.ones([image.shape[1],1]), np.array(range(image.shape[1]))).T +0.5 # Y coordinates of pixel centres in pixels

pixelStringSize = stringWidth / pixelSize # Width of string in pixels

stringImgs = np.empty([N,N, image.shape[0],image.shape[1]]) # Empty placeholders for all the inputs

for n1 in range(N):
    for n2 in range(n1+1, N): # Images are symmetric about the same two nails
        nail1_world, nail2_world = nailCoors[n1, :], nailCoors[n2, :] # Real world coordinates with origin at base centre
        nail1 = nail1_world / pixelSize + max(image.shape)/2 # Image coordinates with origin bottom left
        nail2 = nail2_world / pixelSize + max(image.shape)/2
        
        dist = np.abs((nail2[1]-nail1[1])*pixelCoorsX - (nail2[0]-nail1[0])*pixelCoorsY + nail2[0]*nail1[1] - nail2[1]*nail1[0]) / np.sqrt((nail2[1]-nail1[1])**2 + (nail2[0]-nail1[0])**2)
            # https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
            # Distance of each pixel to the string
        
        stringImgs[n1, n2, :,:] = 255 # Fill with white
        
        stringTouched_idx = np.where(dist < pixelStringSize) # Pixels that should show the string
        stringImgs[n1, n2, stringTouched_idx[0],stringTouched_idx[1]] = dist[stringTouched_idx[0],stringTouched_idx[1]]/pixelStringSize * 255 # Linear darkening of the pixel
        
outsideNails_idx = np.where(np.sqrt((pixelCoorsX-np.mean(pixelCoorsX))**2 + (pixelCoorsY-np.mean(pixelCoorsY))**2) >= np.mean(pixelCoorsX)) # Pixels outside of nails circle
stringImgs[:, :, outsideNails_idx[0],outsideNails_idx[1]] = 255 # Make them white

stringImgs[0, 0, :,:] = 255 # Transform at least one fully white later

#plt.figure()
#plt.imshow(np.sum(stringImgs[0, 1:N, :,:] -255,0)+255, cmap='gray', vmin=0, vmax=255)
##plt.imshow(stringImgs[0, 1, :,:], cmap='gray', vmin=0, vmax=255)
#plt.show()

###################################################################
#%% ############### Transform Possible Lines ######################
###################################################################

stringSinograms = np.empty([N,N, sinogram.shape[0],N])

for n1 in range(N):
    for n2 in range(n1+1, N):
        stringSinograms[n1,n2,:,:] = radon(stringImgs[n1, n2, :,:], theta=theta) # Radon transform

stringSinograms[0,0,:,:] = radon(stringImgs[0, 0, :,:], theta=theta) # Radon transform blank

#%% Start of showing string and radon transform side by side
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
ax1.set_title("Original")
ax1.imshow(stringImgs[0, 0, :,:], cmap=plt.cm.Greys_r)
dx, dy = 0.5 * 180.0 / N, 0.5 / sinogram.shape[0]
ax2.set_title("Radon transform\n(Sinogram)")
ax2.set_xlabel("Projection angle (deg)")
ax2.set_ylabel("Projection position (pixels)")
ax2.imshow(
    stringSinograms[0,15,:,:],
    cmap=plt.cm.Greys_r,
    extent=(-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),
    aspect='auto',
)
fig2.tight_layout()
plt.show()
# End of showing string and radon transform side by side

###################################################################
#%% ############### Copy Duplicate Transforms #####################
###################################################################

for n1 in range(N):
    for n2 in range(n1+1):
        if n2<n1:
            stringSinograms[n1,n2,:,:] = stringSinograms[n2,n1,:,:]
        elif n1==n2:
            stringSinograms[n1,n2,:,:] = stringSinograms[0,0,:,:]
            
InputSpace = np.reshape(stringSinograms, [N*N, -1])

###################################################################
#%% ############### Find optimal path #####################
###################################################################
# Compute starting nail
#remainingSinogram = sinogram.copy()
#np.linalg.lstsq(InputSpace, remainingSinogram)


#%%
#import time
#start_time = time.time()
#zzz = radon(stringImgs[10, 15, :,:], theta=np.linspace(0, 180.0, 2*N, endpoint=False)) # Radon transform
#print("--- %s seconds ---" % (time.time() - start_time))


