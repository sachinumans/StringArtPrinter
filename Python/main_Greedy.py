# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 13:23:12 2024

@authors: Very fun people
"""

import numpy as np
import imageio as iio
import matplotlib.pyplot as plt
from skimage.measure import block_reduce
from scipy.optimize import nnls

### Physical Parameters 
N = 150 # Number of equally spaced nails on the base
R = 300 # Radius of circular base in mm; centre to nails
stringWidth = 1.5 # String thickness in mm

density_coeff = 0.15 # Multiplier for thread effect per pass

### Input Image Path
ImgPath = "../TestImages/Gunter_cropped_smaller_thick.png"

###################################################################
#%% ################## Read Image ######################
###################################################################

image = np.array(iio.imread(ImgPath)); # Read image from file
image = np.mean(image, axis=2)

#subsample image
image = block_reduce(image, block_size=(6,6), func=np.mean)

plt.ion()  # turning interactive mode on
#plt.imshow(image, cmap='gray', vmin=0, vmax=255)
#plt.pause(0.1)

###################################################################
#%% ################## Get Images of All Lines #####################
###################################################################

nailAngles = np.linspace(0, 2*np.pi, num=N, endpoint=False) # Polar angle coordinates of the nails
nailCoors = R * np.array([np.cos(nailAngles), np.sin(nailAngles)]).T

# Start of plotting nail positions
fig2, (ax1, ax2) = plt.subplots(1,2)
ax1.imshow(image, cmap='gray', vmin=0, vmax=255)
ax2.scatter(nailCoors[:,0], nailCoors[:,1])
ax2.set_aspect('equal')
ax2.set_title('Nail coordinates in world measurements')
plt.show()
plt.pause(0.1)
# End of plotting nail positions

pixelSize = 2*R / max(image.shape) # Size of a pixel in mm
pixelCoorsX = np.kron(np.ones([image.shape[0],1]), np.array(range(image.shape[0]))) +0.5 # X coordinates of pixel centres in pixels
pixelCoorsY = np.kron(np.ones([image.shape[1],1]), np.array(range(image.shape[1]))).T +0.5 # Y coordinates of pixel centres in pixels

pixelStringSize = max(stringWidth / pixelSize, 1) # Width of string in pixels

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

dictSize = np.int((N-1)*N/2)
stringImgs_flat = np.empty([stringImgs.shape[2]*stringImgs.shape[3], dictSize])
stringDict = np.empty([3, dictSize]) # Index to keep track of what collums correlate to what strings
idx = 0

for n1 in range(N):
    for n2 in range(n1+1, N): # Images are symmetric about the same two nails
        stringImgs_flat[:,idx] = np.reshape(stringImgs[n1, n2, :,:], [-1])
        stringDict[:, idx] = [idx, n1, n2]
        idx +=1

print("Done getting string images")
###################################################################
#%% ################ Show some #######################
###################################################################

#fig2 = plt.figure()
##plt.imshow(np.sum(stringImgs[10,:,:,:]-255, axis=0)+255, cmap='gray', vmin=0, vmax=255)
#plt.imshow(stringImgs[4,10,:,:], cmap='gray', vmin=0, vmax=255)
#plt.show()
#plt.pause(0.1)

###################################################################
#%% #################### Get Nail Order ###########################
###################################################################
RemainingImage_flat = np.reshape(image.copy(), [-1])

# Invert images
RemainingImage_flat_inv = np.abs(255-RemainingImage_flat)
stringImgs_flat_inv = np.abs(255-stringImgs_flat)


# Initialise
sol = nnls(stringImgs_flat_inv, RemainingImage_flat_inv)
approxImg = stringImgs_flat_inv @ sol[0]
plt.figure()
plt.title('Most perfect outcome, non-integer')
plt.imshow(np.reshape(255-approxImg,image.shape), cmap='gray', vmin=0, vmax=255)
plt.pause(0.1)

startingInput_idx = np.argmax(sol[0])
NailOrder = stringDict[[1,2], startingInput_idx]

img_fnorm = np.linalg.norm(RemainingImage_flat_inv, ord=2)
RemainingImage_flat_inv = RemainingImage_flat_inv - density_coeff* stringImgs_flat_inv[:, startingInput_idx]
newResidual = np.linalg.norm(RemainingImage_flat_inv, ord=2)
oldResidual = newResidual +1
nextInput_idx = -1

StringArt = np.ones(image.shape)*255
fig3, (ax1, ax2) = plt.subplots(1,2)
ax1.set_title('String Art')
s_a = ax1.imshow(StringArt, cmap='gray', vmin=0, vmax=255)
ax2.set_title('Remaining image')
r_i = ax2.imshow(StringArt, cmap='gray', vmin=0, vmax=255)
plotTkn = 0

print("Initialised greedy algorithm")

while newResidual/img_fnorm > 0.10:# and newResidual >= oldResidual+0.1:
    validNails_idx = np.array(np.where(np.logical_and(np.any(stringDict[[1,2], :] == NailOrder[-1], axis=0), stringDict[0, :] != nextInput_idx))).squeeze()
    sol = nnls(stringImgs_flat_inv[:,validNails_idx], RemainingImage_flat_inv)
    nextInput_idx = np.int(stringDict[0, validNails_idx[np.argmax(sol[0])]])
#    print(stringDict[[1,2], nextInput_idx])
    nextNail_idx = np.array(np.where(stringDict[[1,2], nextInput_idx] != NailOrder[-1])) +1
    nextNail = stringDict[nextNail_idx, nextInput_idx]
    
    NailOrder = np.append(NailOrder, nextNail)

    RemainingImage_flat_inv = RemainingImage_flat_inv - density_coeff* stringImgs_flat_inv[:, nextInput_idx]
    oldResidual = newResidual
    newResidual = np.linalg.norm(RemainingImage_flat_inv, ord=2)
#    print(newResidual)

    StringArt += (stringImgs[int(stringDict[1, nextInput_idx]), int(stringDict[2, nextInput_idx])] - 255) * density_coeff
    
    if plotTkn > 10:
        ax1.imshow(StringArt, cmap='gray', vmin=0, vmax=255)
        ax2.imshow(255-np.reshape(RemainingImage_flat_inv, image.shape), cmap='gray', vmin=0, vmax=255)
        plt.pause(0.1)
        pltTkn = 0
    else:
        plotTkn += 1


###################################################################
#%% ################ Show String Art #######################
###################################################################

fig4, (ax1, ax2) = plt.subplots(1,2)
ax1.set_title('Original')
ax2.set_title('String Art')

ax1.imshow(image, cmap='gray', vmin=0, vmax=255)
ax2.imshow(StringArt, cmap='gray', vmin=0, vmax=255)

plt.show()












