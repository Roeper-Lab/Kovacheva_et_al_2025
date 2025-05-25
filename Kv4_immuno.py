"""
Created on Sun Apr  4 21:51:20 2021

@author: LK

v1.3
"""

import sys
import os
import numpy as np
from matplotlib import pyplot as plt

from scipy.ndimage import label, binary_fill_holes
from scipy.stats import ks_2samp  # returns KS-stat, p-value


from skimage import io, color, measure
from skimage.filters import threshold_otsu
from skimage.morphology import erosion, square

def  Kv4_membr_Stat (membrArea, membrKv, mfROI = 1.2, mfBG = 0.6):
    # membrArea = img[ROI[ROI_no].slice]
    connectivity = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    
    thresh_Area = threshold_otsu(membrArea)
    bin_imgHigh = membrArea > thresh_Area * mfROI
    fill_imgHigh = binary_fill_holes(bin_imgHigh)
    bin_imgSmall = erosion(fill_imgHigh, square(10))  # -> cytoplasm & nucleus
    div = fill_imgHigh.copy()
    div [bin_imgSmall == True] = 0  # -> thresholded - erosion image = membrane area
        
    # membrKv = c2_img[ROI[ROI_no].slice]
    
    ## step 1 - membrane analysis ----------------------------
    
    membr_mask, num_membr_labels = label(div, structure=connectivity)
    membr_ROI = measure.regionprops(membr_mask, membrKv)
    
    # find biggest ROI 
    membr_biggest = 0
    membr_size = 0
    
    for i in range(0, num_membr_labels):
        if membr_size < membr_ROI[i].area:
            membr_size = membr_ROI[i].area
            membr_biggest = i
    
    # get intenisties of membrane signal
    tt = membr_ROI[membr_biggest].intensity_image > 0
    meanMembrKv4 = np.mean(membr_ROI[membr_biggest].intensity_image[tt==True])
    areaMembr = membr_ROI[membr_biggest].area
    # medianMembrKv4 = np.median(membr_ROI[membr_biggest].intensity_image[tt==True])
    
    
    # step 2 - cytoplasm (c) analysis ----------------------------

    c_bin = bin_imgHigh.copy()
    c_bin[div == 1] = 0
    c_mask, c_num = label(c_bin, structure=connectivity)
    c_ROI = measure.regionprops(c_mask, membrKv)
    
    # find biggest cluster 
    c_biggest = 0
    c_size = 0
    c_MeanInten = 0
    
    for i in range(0, c_num):
        if c_size < c_ROI[i].area:
            c_size = c_ROI[i].area
            c_biggest = i
            c_MeanInten = c_ROI[i].mean_intensity
    if not c_ROI:
        return meanMembrKv4, areaMembr, membr_ROI[membr_biggest].intensity_image, 0, 0 , []
    else:
        return meanMembrKv4, areaMembr, membr_ROI[membr_biggest].intensity_image, c_MeanInten, c_size, c_ROI[c_biggest].intensity_image



def Kv4analysis (path, path2, fileName, neuronSize = 400, mfROI = 1.2, mfBG = 0.6):
    # path = "g6_A03_r1s6_60x_left01c1.tif" 
    # path2 = "g6_A03_r1s6_60x_left01c2.tif"
    save_path = os.path.dirname(os.path.abspath(__file__))  # "C:/Users/LK/Dropbox (Personal)/6-OHDA-elife-Paper/Kv4.3/PythonAnalysis" # 
    name = fileName
    
       
    connectivity = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    
    # Generate 1000 random colors.
    # It is a list of tuples.
    # Ex. [(0, 0, 0), (1, 0, 1) ... (float, float, float)] -> Values are floats from 0-1
    clrs = np.random.random((1000, 3))
    
    # Step 1 - import images
    img =  io.imread(path, as_gray=True) # img_as_ubyte(io.imread(path, as_gray=True)) # - Convert the colors to 0-1 base.
    c2_img = io.imread(path2, as_gray=True)  # img_as_ubyte(io.imread(path2, as_gray=True))  #  - Converts the colors to 0-1 base.
    print(img.max())
    print (c2_img.max())
    
    
    
    # Step 2 - Generate Threshold using OTSU algo
    thresh = threshold_otsu(img)
    bin_img = img > thresh * mfROI
    bg_img = img < thresh * mfBG
    bg_img = erosion(bg_img, square(10))
    
    # Default is 4 connectivity - https://sites.ualberta.ca/~ccwj/teaching/image/morph/Figs/PNG/connectivity.png
    # find and color TH ROIs
    label_mask, num_labels = label(bin_img, structure=connectivity)
    ROI = measure.regionprops(label_mask, c2_img)
    colored_img = color.label2rgb(label_mask, bg_label=0, colors=clrs)
    
    
    # find and analyse BG ROIs
    bg_mask, bgArea_num = label(bg_img, structure=connectivity)
    bg_data = measure.regionprops(bg_mask, c2_img)
    colored_BGimg = color.label2rgb(bg_mask, bg_label=0, colors=clrs)
    
    """
    `clusters` in an array of `RegionProperties` objects.
    Each object has these values: (Debug with breakpoint and see the values)
    
    ['area', 'bbox', 'bbox_area', 'centroid', 'convex_area', 'convex_image', 'coords', 'eccentricity', 
    'equivalent_diameter', 'euler_number', 'extent', 'feret_diameter_max', 'filled_area', 'filled_image', 'image', 
    'inertia_tensor', 'inertia_tensor_eigvals', 'intensity_image', 'label', 'local_centroid', 'major_axis_length', 
    'max_intensity', 'mean_intensity', 'min_intensity', 'minor_axis_length', 'moments', 'moments_central', 'moments_hu', 
    'moments_normalized', 'orientation', 'perimeter', 'perimeter_crofton', 'slice', 'solidity', 'weighted_centroid', 
    'weighted_local_centroid', 'weighted_moments', 'weighted_moments_central', 'weighted_moments_hu', 
    'weighted_moments_normalized']
    """
    # Measure the data. Pass the original image as second param.
    # clusters = measure.regionprops(label_mask, img)
    
    # Measure the data. Pass the original image as second param.
    # ROI = measure.regionprops(label_mask, c2_img)
    
   
   # Show TH original TH - masks images 
    
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(15, 17))
    axes = ax.ravel()
    
    axes[0].imshow(img, cmap="gray")
    axes[0].set_title("Original Image")
    
    axes[1].imshow(c2_img, cmap="gray")
    axes[1].contour(label_mask, [0.5], linewidths=0.4, colors='r')
    axes[1].contour(bg_mask, [0.5], linewidths=0.4, colors='g')
    axes[1].set_title("TH-ROIs on Kv4.3")
    
    axes[2].imshow(colored_img)
    axes[2].set_title("TH ROIs")
    
    axes[3].imshow(colored_BGimg)
    axes[3].set_title("BG ROIs")
    
    for x in axes:
        x.axis("off")
    
    fig.suptitle(f'{name} - TH ROIs')   
    plt.savefig(f"{save_path}{name}_TH_ROIs-v2.jpg", dpi=450, format='jpg')
    plt.savefig(f"{save_path}{name}_TH_ROIs-v2.svg", dpi=450, format='svg')
    plt.tight_layout()
    plt.show()
     
   
    
   
    
   # ----------------------- begin neuron analysis -------------------------
    
    # min NEURON size - 400 pixels
    ROI_no = np.array([], dtype=int)
    
    # find ROIs above neuronSize
    for i in range(0, num_labels):
        if ROI[i].area > neuronSize:
            ROI_no = np.append(ROI_no, i)
            
    
    # shows asreas bigger than the neuronSize-value
    if np.size(ROI_no) > 1: 
        fig, ax = plt.subplots(nrows = 6, ncols=np.size(ROI_no), figsize=(3*np.size(ROI_no), 9))
        # axes = ax.ravel()
        nplot = 0
        meanMembr = np.array([], dtype=int)  # membran
        areaMembr = np.array([], dtype=int)
        
        meanC = np.array([], dtype=int)  # cytoplasma
        areaC = np.array([], dtype=int)
        
        for i in ROI_no:
            tempMean, tempArea, tempMembIntensity, tempCPmean, tempCParea, tempCProi = Kv4_membr_Stat (img[ROI[i].slice],
                                                                          c2_img[ROI[i].slice])
            meanMembr = np.append(meanMembr, tempMean)
            areaMembr = np.append(areaMembr, tempArea)
            
            if tempCParea > 50: # only bigger ones
                meanC = np.append(meanC, tempCPmean)
                areaC = np.append(areaC, tempCParea)
                ax[5][nplot//6].imshow(tempCProi)
                ax[2][nplot//6].contour(tempCProi, [0.5], linewidths=0.4, colors='g')
                
            
            ax[0][nplot//6].imshow(~img[ROI[i].slice], cmap="gray")
            ax[0][nplot//6].set_title(f"{i} - TH")
            ax[0][nplot//6].axis("off")
                    
            ax[1][nplot//6].imshow(ROI[i].image)
            ax[1][nplot//6].set_title("TH-mask")
            ax[1][nplot//6].axis("off")
                         
            ax[2][nplot//6].imshow(~c2_img[ROI[i].slice], cmap="gray")
            ax[2][nplot//6].set_title("Kv4")
            ax[2][nplot//6].contour(ROI[i].intensity_image, [0.5], linewidths=0.4, colors='r')
            ax[2][nplot//6].axis("off")
            
            ax[3][nplot//6].imshow(ROI[i].intensity_image)
            ax[3][nplot//6].set_title("Kv4-intens")
            ax[3][nplot//6].axis("off")
            
            ax[4][nplot//6].imshow(tempMembIntensity)
            ax[4][nplot//6].set_title("Kv4-membr")
            ax[4][nplot//6].axis("off")
            
            ax[5][nplot//6].set_title("Kv4-cytoplasm")
            ax[5][nplot//6].axis("off")
                    
            nplot += 6
        
        fig.suptitle(f'{name} - TH ROIs > {neuronSize}')   
        plt.tight_layout()
        plt.savefig(f"{save_path}{name}_Kv4_ROIs-v2.jpg", dpi=450, format='jpg')
        plt.savefig(f"{save_path}{name}_Kv4_ROIs-v2.svg", dpi=450, format='svg')
        plt.show()
        
    else: # ----------------------------------- single ROI ----------------------------------------------
    
        fig, ax = plt.subplots(nrows = 6, ncols=1, figsize=(4*np.size(ROI_no), 12))
        axes = ax.ravel()
       
        meanMembr = np.array([], dtype=int)
        areaMembr = np.array([], dtype=int)
        
        meanC = np.array([], dtype=int)
        areaC = np.array([], dtype=int)
        
        for i in ROI_no:
            tempMean, tempArea, tempMembIntensity, tempCPmean, tempCParea, tempCProi = Kv4_membr_Stat (img[ROI[i].slice],
                                                                          c2_img[ROI[i].slice])
            meanMembr = np.append(meanMembr, tempMean)
            areaMembr = np.append(areaMembr, tempArea)
            
            if tempCParea > 50: # only bigger ones
                meanC = np.append(meanC, tempCPmean)
                areaC = np.append(areaC, tempCParea)
                ax[5].imshow(tempCProi)
                ax[2].contour(tempCProi, [0.5], linewidths=0.4, colors='g')
                
            
            ax[0].imshow(img[ROI[i].slice], cmap="gray")
            ax[0].set_title(f"{i} - TH")
            ax[0].axis("off")
                    
            ax[1].imshow(ROI[i].image)
            ax[1].set_title("TH-mask")
            ax[1].axis("off")
                         
            ax[2].imshow(c2_img[ROI[i].slice], cmap="gray")
            ax[2].set_title("Kv4")
            ax[2].contour(tempMembIntensity, [0.5], linewidths=0.4, colors='r')
            ax[2].axis("off")
            
            ax[3].imshow(ROI[i].intensity_image)
            ax[3].set_title("Kv4-intens")
            ax[3].axis("off")
            
            ax[4].imshow(tempMembIntensity)
            ax[4].set_title("Kv4-membr")
            ax[4].axis("off")
            
            ax[5].set_title("Kv4-cytoplasm")
            ax[5].axis("off")
            
        fig.suptitle(f'{name} - TH ROIs > {neuronSize}')                    
        plt.tight_layout()
        plt.savefig(f"{save_path}{name}_Kv4_ROIs.jpg", dpi=450, format='jpg')
        plt.savefig(f"{save_path}{name}_Kv4_ROIs.svg", dpi=450, format='svg')
        plt.show()
        
    
    ## export t1 == TH-ROI Kv4-signal intensity classic
    f = open(f"{save_path}t1_{fileName}.csv",'w')
    p = 0
    for i in range(0, num_labels):
        if ROI[i].area > 50:  # remove very small areas
            f.write(f"{p},{ROI[i].area},{ROI[i].mean_intensity}\n")
            p +=1
    f.close()
    
    ## export t2 = TH-membrane ROI Kv4-signal intensity
    f = open(f"{save_path}t2_{fileName}.csv",'w')
    for i in range(0, np.size(meanMembr)):
        f.write(f"{i},{areaMembr[i]},{meanMembr[i]}\n")
    f.close()
    
    ## export t3 = TH-cytoplasm ROI Kv4-signal intensity
    f = open(f"{save_path}t3_{fileName}.csv",'w')
    for i in range(0, np.size(meanC)):
        f.write(f"{i},{areaC[i]},{meanC[i]}\n")
    f.close()
    
    ## export t4 = background ROIs Kv4 intensity
    f = open(f"{save_path}t4_{fileName}.csv",'w')
    for i in range(0, bgArea_num):
        f.write(f"{i},{bg_data[i].area},{bg_data[i].mean_intensity}\n")
    f.close()
    


search_path = os.path.dirname(os.path.abspath(__file__))
c = 0
for file in os.listdir(search_path):
    if  file.endswith("c1.tif") :
        # print(file) # file.endswith('6_C2.tif'): # 
        t1 = f"{search_path}/{file}" # f"{file}" # 
        t2 = f"{search_path}/{file[:-5]}2.tif" #  f"{file[:-5]}3.tif" # 
        print(t1)
        print(t2)
        Kv4analysis (t1, t2, file[:-6], neuronSize=2000, mfROI = 1.2, mfBG = 0.6)  # mf = mutiplication factor  
        c +=1

print(c)
