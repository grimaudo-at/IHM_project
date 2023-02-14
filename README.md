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

By the end of the day today, I had managed to get to the point where I was able to merge all transmitter dataframes. I am now in the process of calibration. I've calculated all available T and D values (see the Reeder paper for how calibration curves are calculated) for each logger, but I need to move the datasets from a long to wide format to actually construct the calibration curves. That's my objective for tomorrow. Qualitatively, I did notice that there is a lot of variation in how a logger's readings differed from the thermocouples'. ***Some loggers were spot on whereas others were reading 2 or 3 degrees colder than the thermocouple, which is a substantial amount.*** They never read warmer than the thermocouple. Hopefully the calibration method is enough to sufficiently correct for this error. 

I am close to being able to calibrate the transmitters. Currently, I've devided the transmitter data up into the three calibration groups and have summarized both T and D values for each. However, these dataframes are in long format, which will make it difficult to run their T and D values through the calibration equations. I should first coerce the dataframes into wide format, with columns for all T and D values. 

# 2/14/2023
Today, I managed to coerce the calibration dataframes into wide format and calculate the quadratic equation parameters for group 1 and 2. I then calculated deviations between D1 and T1, D2 and T2, and D3 and T3 for group 1 and 2. I built two regressions to see how well the deviation between D1 and T1 could predict the deviations at T2 and T3. I used these regressions to predict what the derivation values were for group 3 at time 2 and 3 because only D1 was available for that group. The slopes of these regressions was not 1, indicating assuming a 1:1 change in deviations between T1, T2, and T3 is probably not appropriate. Then, assuming the thermocouple temperature at T2 and T3 was 23 and 37 C (temps used for calibration treatment, which was not available because a thermocouple was not available for this group), respectively, I used the regressions to estimate what the D2 and D3 values might have been for transmitters in this group. Then, using all of this data, I calibrated the temperature data. 

An important observation is that while some transmitters barely deviated from the thermocouple, others deviated significantly, sometimes up to several degrees. The calibration method, however performed very well at estimating the true temperature value. While group 3's calibration parameters had to be estimated from regression, I think I maximized the chance it was as close to the true values as possible. 

Once I calibrated all the raw data, I trimmed the datasets to not include any data from dates ***less than or equal to*** the date of deployment OR ***greater than or equal to*** the date of retrieval. 
