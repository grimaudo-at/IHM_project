# IHM_project

# 2/13/2023
Today, initializing a new GitHub repository for the IHM project. Below is a list of this repository's current contents: 
 
1) IHM Project.Rproj -- The R project linked to GitHub
2) movement_exploratory.R -- An old script I used to explore some site-section level movement patterns between early and late hibernation. 
3) transmitter_ID_assignments.R -- An old script I used to randomly assign transmitters to bats. I don't think we actually ever used this in the field. 
4) transmitter_recover_summary.R -- An R script used to summarize the retrieval and success rate of all transmitters using the .csv file containing transmitter metadata found here: Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_metadata.csv
5) transmitter_data_cleaning_calibration.R -- An R script used to take the raw, uncleaned transmitter data and transform it into cleaned, calibrated data. 

My current objective is the cleaning and calibrating of all the transmitter data we downloaded last spring. There are 80 .csv files of transmitter data stored in "Data/IHM Project/Transmitter data" in my Dropbox folder. Each .csv contains several lines of metadata before the raw data actually starts. In this metadata are values for the model type and serial number of the logger, information I'd like to incorporate into the final transmitter dataframe. The ultimate objective is to have a single .csv file that contains all of the cleaned, transformed transmitter data, with the following columns of data: site, transmitter_id, serial_num, date, time, temp_c, model. I will accomplish this in the script "transmitter_data_cleaning_calibration.R", which will broadly follow the below workflow:

1) Read in all .csv files and store in list of dataframes.
2) From the metadata in each transmitter .csv, extract the model version and serial number. These will be stored in lists and re-incorporated into dataframes. Remove the rest of the metadata so that only raw data remains.
3) Once only the raw data remains, it needs to be coerced into a wider format. As it is read in, there is one column of raw data in the dataframe with repeating values of data, "C", and the temperature value. These need to be separated into 3 columns. 
4) Once the raw data is separated into its three columns for each dataframe, the model version and serial number data can be re-incorporated from the lists I made earlier. Because I am manipulating all of these dataframes while they are stored in a list, I think that I should be able to match raw data dataframes to metadata lists by using their position in the list.
5) Use the "transmitter_metadata.csv" metadata file to match in information on site, transmitter ID, and calibration group, matching on serial #. 
6) Merge all dataframes. 
7) Using data on when calibration trials took place, use each transmitter's data to calibrate them. Transmitters in group 3 did not receive the full calibration treatement and their calibration parameters will have to be estimated. Look at physical notebook for ideas on how to do this. 
8) Once all the data is calibrated, use date information in the "transmitter_metadata.csv" file to trim away all data recorded when the transmitter was off of bats. 

By the end of the day today, I had managed to get to the point where I was able to merge all transmitter dataframes. I am now in the process of calibration. I've calculated all available T and D values (see the Reeder paper for how calibration curves are calculated) for each logger, but I need to move the datasets from a long to wide format to actually construct the calibration curves. That's my objective for tomorrow. Qualitatively, I did notice that there is a lot of variation in how a logger's readings differed from the thermocouples'. **Some loggers were spot on whereas others were reading 2 or 3 degrees colder than the thermocouple, which is a substantial amount.** They never read warmer than the thermocouple. Hopefully the calibration method is enough to sufficiently correct for this error. 