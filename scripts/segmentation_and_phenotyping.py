#!/bin/python3

# Import required packages
import numpy as np
import pandas as pd
import re
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from skimage import io, segmentation
from skimage.filters import threshold_otsu, threshold_triangle
from tifffile import TiffFile
from xml.etree import ElementTree
import argparse


#
# 0. DEFINING FUNCTIONS
#

# Generate function to get channel names and numbers
def get_channel_names(tif_path):

    channel_dict = {} # Create an empty dictionary to store the channel names and number
    channel_counter = 0 # Create a counter to keep track of the channel number

    # Open the tif file
    with TiffFile(tif_path) as tif:

        for page in tif.series[0].pages:
            
            # Add channel to counter
            channel_counter += 1
            # Adding channel name and number to dictionary
            channel_dict[ElementTree.fromstring(page.description).find('Name').text.split(" ")[0].upper()] = channel_counter
    
    return channel_dict


#
# 1. PARSING COMMAND LINE ARGUMENTS
#

# Initialize argparse
parser = argparse.ArgumentParser(description="Segmentation and phenotyping of IHC images")

# Add arguments
parser.add_argument("-p", "--paths_to_image", type=str, help="Paths to the mIHC image to be segmented and phenotyped separated by :")
parser.add_argument("-o", "--output_path", type=str, help="Path to the output folder to store cell phenotypes")
parser.add_argument("-s", "--segmentation_path", type=str, help="Path to the output folder to store segmentation masks")

# Parse arguments
args = parser.parse_args()

paths_to_images = args.paths_to_image.split(":")


for path_to_image in paths_to_images:

#
# 2. NUCLEI SEGMENTATION, MASK EXPANSION AND SIZE MEASUREMENT
#

    print(path_to_image)

    # Get channel names and numbers
    channel_dict = get_channel_names(path_to_image)

    # Load the image
    mIHC_image = io.imread(path_to_image)

    # Applying StarDist to segment the nuclei
    model = StarDist2D.from_pretrained('2D_versatile_fluo') # Load the pretrained model
    labels, labels_dic = model.predict_instances(normalize(mIHC_image[channel_dict["DAPI"]-1])) # Segment the nuclei based on the normalized DAPI channel

    # Expanding nuclear masks to include the cytoplasm
    expanded_labels = segmentation.expand_labels(labels, 2) # Expand the nuclear masks by 2 pixels

    # Store parameters of the segmeneted nuclei

    # Generate pandas dataframe to store the parameters
    num_of_cells = len(np.unique(labels)) - 1 # Number of cells (substract 1 to exclude background)
    cells_df = pd.DataFrame(columns=["Cell_ID", "X_coor", "Y_coor", "Nuclear_area"], index=range(num_of_cells))

    # Getting cell coordinates and nuclear area

    x_coor = [coor[0] for coor in labels_dic['points']] # X coordinates of the cells
    y_coor = [coor[1] for coor in labels_dic['points']] # Y coordinates of the cells
    nuclear_area = [np.sum(labels == i) for i in range(1, num_of_cells +1)] # Nuclear area of the cells

    # Store the parameters in the dataframe
    cells_df["Cell_ID"] = range(1, num_of_cells + 1)
    cells_df["X_coor"] = x_coor
    cells_df["Y_coor"] = y_coor
    cells_df["Nuclear_area"] = nuclear_area

    # Generate segmentation masks (original and expanded) as tiff objects

    # Generate output filename

    pattern = r'Core\[\d+,\d+,[A-Z]\]' # Define the regular expression pattern
    match = re.search(pattern, path_to_image) # Search for the pattern in the filename

    io.imsave("{}/{}_nuclear_masks.tiff".format(args.segmentation_path, match.group(0)), labels)
    io.imsave("{}/{}_whole_cell.tiff".format(args.segmentation_path, match.group(0)), expanded_labels)

#
# 3. INTENSITY MEASUREMENT
#

    # Determining intensities of the remaining channels for each cell and append to the dataframe

    # Create columns for each channel
    for channel in channel_dict.keys():
        if channel != "AUTOFLUORESCENCE": 
            cells_df[channel] = np.nan 

    # Get the intensity of each channel for each cell. Using the expanded labels for CD8, CD20, CD68 and CD4 and the original labels for DAPI and FOXP3 
    for i in range(1, num_of_cells + 1):
        for channel in channel_dict.keys():
            if channel == "AUTOFLUORESCENCE":
                pass 
            elif channel == "DAPI" or channel == "FOXP3":
                cells_df.loc[i-1, channel] = np.mean(mIHC_image[channel_dict[channel] - 1][labels == i])
            else:
                cells_df.loc[i-1, channel] = np.mean(mIHC_image[channel_dict[channel] - 1][expanded_labels == i])



#
# 4. PRELIMINARY PHENOTYPING
#

    # Setting a preliminary threshold for phenotyping
    prel_threshold = 1.2

    # Create new column in dataframe to store the phenotypes to store string values

    cells_df["Preliminary_phenotype"] = ""

    # Phenotyping based on the intensity of CD8, CD20, FOXP3, PAN.CK and CD4
    for cell in range(0, num_of_cells):
        #Assign phenotype to Others if the cell does not express any of the markers
        if all(cells_df.loc[cell, ["CD8", "CD20", "PAN-CK", "CD4", "CD68"]] < prel_threshold):
            cells_df.loc[cell, "Preliminary_phenotype"] = "Other"
        #Assign phenotype to the marker with the maximum intensity
        else:
            # Get the main marker
            main_marker = cells_df.loc[cell, ["CD8", "CD20", "PAN-CK", "CD4", "CD68"]].idxmax()

            # Determine if FOXP3 is expressed
            if cells_df.loc[cell, "FOXP3"] > prel_threshold:
                cells_df.loc[cell, "Preliminary_phenotype"] = main_marker + "_FOXP3"
            else:
                cells_df.loc[cell, "Preliminary_phenotype"] = main_marker



#
# 5. CHOOSING THRESHOLDING METHOD BASED ON THE PRELIMINARY PHENOTYPING
#

    # Generature a dictionary to store the estimated counts of each phenotype
    phenotype_counts = cells_df["Preliminary_phenotype"].value_counts().to_dict()

    # Get channel names to determine thresholds
    channels_to_threshold = list(channel_dict.keys())

    # Remove the channels that are not used for phenotyping
    channels_to_threshold.remove("DAPI")
    channels_to_threshold.remove("AUTOFLUORESCENCE")


    # Create a dictionary to store the thresholds
    channel_thresholds = {k:v for k, v in zip(channels_to_threshold, [0]*len(channels_to_threshold))}

    # Determine the thresholds based on the cell counts obtained in the preliminary phenotyping
    for channel in channels_to_threshold:
        
        # Getting the phenotypes that contain the channel name
        channel_phenotypes = [phenotype for phenotype in phenotype_counts.keys() if channel in phenotype]

        # Getting the sum of the cells identified positive for the marker
        positive_cells = sum([phenotype_counts[phenotype] for phenotype in channel_phenotypes])

        print(channel, positive_cells)

        #Using otsu thresholding for PAN-CK
        if channel == "PAN-CK" and positive_cells > 100:
            channel_thresholds[channel] = threshold_otsu(mIHC_image[channel_dict[channel] - 1])

        elif positive_cells > 100:
            channel_thresholds[channel] = threshold_triangle(mIHC_image[channel_dict[channel] - 1])

        elif positive_cells > 0:
            # Get the minimum intensity of the cells expressing the marker minus 20% of the intensity
            min_intensity = min([cells_df.loc[cell, channel] for cell in range(0, num_of_cells) if channel in cells_df.loc[cell, "Preliminary_phenotype"]])

            # Set the threshold to the minimum intensity
            channel_thresholds[channel] = min_intensity - 0.2 * min_intensity

        else:
            channel_thresholds[channel] = 1.2

    print(channel_thresholds)

    # Using 1.6 as the minimum threshold for pan-ck
    if channel_thresholds["PAN-CK"] < 1.2:
        channel_thresholds["PAN-CK"] = 1.2


    # Calculate otsu ond triangle thresholds for each channel and store them

    # Generate new dataframe to store the thresholds of each channel and core
    #thresholds_df = pd.DataFrame(columns=["Channel", "Otsu_threshold", "Triangle_threshold"], index=channels_to_threshold)


    # for channel in channels_to_threshold:

    #    # Calculate the otsu threshold
    #     otsu_threshold = threshold_otsu(mIHC_image[channel_dict[channel] - 1])

    #     # Calculate the triangle threshold
    #     triangle_threshold = threshold_triangle(mIHC_image[channel_dict[channel] - 1])

    #     # Adding thresholds to the dataframe
    #     thresholds_df.loc[channel, "Channel"] = channel
    #     thresholds_df.loc[channel, "Otsu_threshold"] = otsu_threshold
    #     thresholds_df.loc[channel, "Triangle_threshold"] = triangle_threshold


    # print(thresholds_df)
    # # Saving final_df and intensity dataframes as csv files
    # file_name = path_to_image.split("/")[-1].split("_component_data_")[0]
    # # Save the thresholds dataframe
    # thresholds_df.to_csv(args.output_path + file_name + "thresholds.csv", index=False)

#
# 6. DETERMINING PROPORTION OF PIXELS ABOVE THRESHOLD FOR EACH CHANNEL
#

    # Create new columns in the dataframe to store the proportion of pixels above the threshold for each channel
    for channel in channel_thresholds.keys():
        cells_df[channel + "_proportion"] = np.nan

    # Determine the proportion of pixels above the threshold for each channel and cell
    for cell in range(0, num_of_cells):
        for channel in channel_thresholds.keys():

            #If the channel is not nulcear using the expanded labels
            if channel != "FOXP3":
                cells_df.loc[cell, channel + "_proportion"] = np.sum(mIHC_image[channel_dict[channel] - 1][expanded_labels == cell + 1] > channel_thresholds[channel]) / np.sum(expanded_labels == cell + 1)

            #If the channel is nuclear using the original labels
            else:
                cells_df.loc[cell, channel + "_proportion"] = np.sum(mIHC_image[channel_dict[channel] - 1][labels == cell + 1] > channel_thresholds[channel]) / np.sum(labels == cell + 1)

#
# 7. DEFINITIVE PHENOTYPING 
#

    # Check if the cells are positiove for each marker. The mean intensity, and at least 50% of the pixels included in the mask should be above the threshold
    # If there are conflictiveresult, the one with the highest intensity will be chosen

    # Create new column in dataframe to store the phenotypes to store string values
    cells_df["Phenotype"] = ""

    # Phenotyping based on the intensity of CD8, CD20, FOXP3, PAN.CK and CD4


    for cell in range(0, num_of_cells):

        # Determine which markers are positive
        positive_markers = [marker for marker in channel_thresholds.keys() if cells_df.loc[cell, marker] > channel_thresholds[marker]]

        # DEfining boolean to check if FOXP3 is expressed
        pos_FOXP3 = False

        # Remove FOXP3 from the list of markers and setting pos_FOXP3 to True
        if "FOXP3" in positive_markers:
            positive_markers.remove("FOXP3")
            pos_FOXP3 = True 

        #Assigning phenotypes to cells

        # The "Others" phenotype is assigned if the cell does not express any of the markers
        if len(positive_markers) == 0: 

            cells_df.loc[cell, "Phenotype"] = "Other"

        # If a single marker is expressed, the phenotype is assigned to the marker
        elif len(positive_markers) == 1:

            # Determine if FOXP3 is expressed
            if positive_markers[0] in ["CD4", "CD8"] and pos_FOXP3:
                cells_df.loc[cell, "Phenotype"] = positive_markers[0] + "_FOXP3"
            else:
                cells_df.loc[cell, "Phenotype"] = positive_markers[0]

        # If multiple markers are expressed, the phenotype is assigned to the marker with the highest intensity
        else:

            # The main marker is the one with the highest intensity proportion of positive pixels of at least 40%

            # Filter the markers that have less than 30% of positive pixels
            positive_markers = [marker for marker in positive_markers if cells_df.loc[cell, marker + "_proportion"] > 0.3]
            

            # If there are no markers with at least 30% of positive pixels, assign the phenotype to Others, else assign the phenotype to the marker with the highest intensity
            if len(positive_markers) == 0:

                cells_df.loc[cell, "Phenotype"] = "Other"

            else:
                
                # Get the main marker
                main_marker = max(positive_markers, key=lambda x: cells_df.loc[cell, x])

                # Determine if FOXP3 is expressed
                if main_marker == "CD4" and pos_FOXP3:
                    cells_df.loc[cell, "Phenotype"] = main_marker + "_FOXP3"
                else:
                    cells_df.loc[cell, "Phenotype"] = main_marker

#
# 8. CHECKING FOR MISSASGNED "OTHER" CELLS
#

    # Reassign "Other" cells to pan-ck if their nucleus is biggern that 200 and the intensity of pan-ck is the highest among the markers
    for cell in range(0, num_of_cells):
        if cells_df.loc[cell, "Phenotype"] == "Other" and cells_df.loc[cell, "Nuclear_area"] > 170:
            if cells_df.loc[cell, "PAN-CK"] > cells_df.loc[cell, ["CD8", "CD20", "CD4"]].max():
                cells_df.loc[cell, "Phenotype"] = "PAN-CK"

# #
# # 9. PHENOTYPING CD68 CHANNEL USING CELLPOSE
# #
                    
#     # Perform cytoplasm segmentation using Cellpose for the CD68 channel
#     # Define the cellpose model
#     model = models.Cellpose(gpu=True, model_type='cyto')

#     # Define the channels to use for cellpose
#     my_channels = [0, 0] # Using 0,0 because only a cytoplamsatic channel will be used

#     # Running cellpose
#     masks, flows, styles, diams = model.eval(mIHC_image[channel_dict["CD68"]-1], diameter=20, channels=my_channels, do_3D=False, flow_threshold=0.45)

#     # Store the masks in a dataframe

#     # Determining the number of masks
#     number_of_masks = len(np.unique(masks)) - 1 # -1 to exclude the background

#     # Create a dataframe to store the cell ID, coordinates, area, intensity and phenotype
#     df_cd68 = pd.DataFrame(columns=["Cell_ID", "X_coor", "Y_coor", "Area", "Intensity", "Phenotype"], index=range(number_of_masks))
#     df_cd68 = df_cd68.astype({"Cell_ID": float, "X_coor": float, "Y_coor": float, "Area": float, "Intensity": float, "Phenotype": str})

#     # Filling the data frame
#     for mask in range(1, number_of_masks+1):

#         df_cd68.loc[mask-1, "Cell_ID"] = float(mask + num_of_cells) # Adding the number of cells to the mask ID to avoid overlapping with the other cell IDs
#         df_cd68.loc[mask-1, "X_coor"] = np.mean(np.where(masks == mask)[1])
#         df_cd68.loc[mask-1, "Y_coor"] = np.mean(np.where(masks == mask)[0])
#         df_cd68.loc[mask-1, "Area"] = len(np.where(masks == mask)[0])
#         df_cd68.loc[mask-1, "Intensity"] = np.mean(mIHC_image[channel_dict['CD68']-1, np.where(masks == mask)])
#         df_cd68.loc[mask-1, "Phenotype"] = "CD68"

#
# 10. MERGING AND SAVING THE DATAFRAMES
#
        
    # Merge the dataframes (cell_id, coordinates and phenotypes)
    #final_df = pd.concat([cells_df.loc[:, ["Cell_ID", "X_coor", "Y_coor", "Phenotype"]], df_cd68.loc[:, ["Cell_ID", "X_coor", "Y_coor", "Phenotype"]]])
    final_df = cells_df.loc[:, ["Cell_ID", "X_coor", "Y_coor", "Phenotype"]]


    # Saving final_df and intensity dataframes as csv files
    file_name = path_to_image.split("/")[-1].split("_component_data_")[0]
    
    # Saving dataframes
    final_df.to_csv(args.output_path + "/" + file_name + "_cells.csv", index=False)
    cells_df.to_csv(args.output_path + "/" + file_name + "_intensity.csv", index=False)
    # df_cd68.to_csv(args.output_path + file_name + "_cd68.csv", index=False)