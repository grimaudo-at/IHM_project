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

# 2/15/2023
Today, I wrote the calibrated master transmitter dataframe to the following file: Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv. 

I have also created the R script entitled "arousal_classification.R". This script is meant to loop through all of the transmitter data and classify each datapoint as belonging to an arousal (1) or torpor bout (0). To do so, I first identified arousals in a broad-sweep using Reeder et al. 2012's definition. Under this definition, a datapoint is classified as an arousal if it falls within 10C of the logger's maximum value. I found that this indeed adequately identifies arousal events. However, in many cases, there is a "ramping up" or "cooling down" period on either side of the identified arousal event containing temperature data that falls outside of the 10C threshold but that is significantly warmer than the torpor bouts on either side of the arousal. Because I think it is inappropriate to allow these values to fall into the "torpor classification", I constructed a for-loop in this script to re-classify them as part of the arousal event. For those "ramping up" values that proceed the identified arousal, I classified them as part of the arousal if the temperature value at i-1 (i being the ramping up value) was at least 1C colder. Similarly, for the "cooling down" values following an arousal event, I classified them as part of the arousal if the temperature value at i+1 was at least 1C colder. I ran this loop twice, but the number of times can be extended. Running the loop twice extends the number of additional datapoints that can be incorporated into the arousal event by two on either side of the arousal. For example, if i was classified as an arousal under Reeder's definition, then both i-1 and i-2 can be classified as part of the arousal as long as they meet the conditions described above. 

In coming days, my goal is to take this new arousal dataset and summarize it, calculating each of the following metrics: 
a) total number of arousals
b) length of each arousal
c) total number of torpor bouts
d) length of each torpor bout
e) number of movements
f) change in temperature following movement
g) average temperature during each torpor bout
h) min/max temperature during each torpor bout/arousal
i) event number = order of the torpor/arousal in the sequence of bouts/arousals. 

# 2/16/2023
Today, I was again working in the "arousal_classification.R" script to identify all arousals and torpor bouts. After running the for-loop that does that (explained in previous entry), which takes about 1.5 hours to run over the whole dataset, I also constructed and ran a for-loop that identifies each torpor/arousal's "event number". The event number is each torpor bout or arousal's position in the sequence of torpor bouts/arousals, and it is how I am going to do much of the summarizing. Once I added these data to the dataframe, I re-wrote the Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv file. 

# 2/17/2023
Talked with Kate today about constructing simulations to explore role of bat movement in determining infection load over hibernation. In all three simulations, I will use a temperature response curve of the pathogen to simulate forward the early hibernation pathogen load on each bat, using its transmitter tempeature data. Simply, every load at time i will be a product of the load at time i-1 and r(temp[i-1]), r being the temperature-specific growth rate. Additionally, we will run this simulation on the temperature dataset collected from the site's psychrometer (which would be analagous to a bat that never moves) and on a dataset that assumes a constant average temperature over the course of hibernation. Comparing the simulated end-of-hibernation pathogen load values to our empirical data will elucidate which framework most closely predicts infection dynamics, with the hypothesized outcome that the simulation over bat movement data will fit best, indicating movement is important driver of infection severity. Comparing values of r over each simulation will further yield information on how bat movement specifically influences pathogen growth. There are three potential simulation approaches that I'll attempt in order:

1) Simulate pathogen load forward deterministically by pulling known, calculated r values from the thermal performance curve that Kate and Skylar fit.

2) Similar to the first method, but r in the equation will be substituted by its derivation, the parameters for which were also calculated by Kate and Skylar's fitting. 

3) Simulate pathogen load forward using the same equation as simulation 2, but iteratively fitting values of r as the simulation progresses. 

# 2/20/2023
First thing I did today was re-write the Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv file, which had changed slightly after removing erroneously high temperature values. After removing very high values, the arousal classification automation seems to have worked pretty well, save for a few instances where there are suspected arousals but they did not get high enough in temperature (or those temperatures were not detected during logging) to fall under the arousal classification. I think these instances might have to just be corrected by hand. 

Furthermore, I tried automating a process to identify movements, which I did under a new section entitled "summarizing event data" in the "arousal_classification.R" script. This section also now contains code on summarizing arousal and torpor data by several metrics. I attempted to identify movements based on the average temperatures of torpor bouts on either side of the arousal. If the mean temperature of arousal bout i **fell outside of *one standard deviation* of the average temperature of the previous arousal bout (i-1)**, then the arousal was classified as being accompanied by a movement. I then plotted this data and colored points based on if they occurred in the same location (a new location was assumed if the algorithm detected a movement). The algorithm to identify movements was able to identify movements in the most extreme cases, but it also failed to detect some small ones and even mis-classified several arousals as being movements. The latter case appeared to occur mostly if there was just natural cooling during winter, as this results in the mean temperature dropping without an actual movement. I suspect I will have to manually identify these cases and fix them by hand. I saved a .pdf of all the plots with points colored by locations and arousal bouts here, which overwrote an older version just colored by arousals and torpor:  "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_arousal_classification_plots.pdf"

Generally this algorithm for identifying movements does not seem to work very well. Tomorrow, I should explore alternate methods. 

# 2/21/2023
Today, I tried running the movement classification algorithm again, but this time using median instead of mean values during torpor bouts for identifying movements. Performed extremely similarly, maybe a little bit better, with still a ton of user correction necessary. I printed all of the new plots to: Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_arousal_classification_plots.pdf
