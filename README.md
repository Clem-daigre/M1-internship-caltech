
### RAW_DATA_FILES contains :

- NEIC_RAWDATA is the catalog from NEIC containing earthquakes of all depths, all magnitudes (only magnitudes above 4 for NEIC),occuring in the Tonga area since 2005.
- ISC_RAWDATA is the combined ISC/IDC/NEIC catalog containing earthquakes occuring in the Tonga area since 2005, of all depths, all magnitudes, unknown depth and unknown magnitude. UNKNOWN DEPTHS ARE WROTE AS 0 DEPTH.
- CMTSOLUTION_RAWDATA is the CMTSOLUTION format catalog from Global CMT containing all the events with Mw>6.5 occuring in the Tonga area since 2005, containing their magnitude info but not their depth.
- PSMECA_RAWDEPTHS is the PSMECA format catalogfrom Global CMT containing all the events with Mw>6.5 occuring in the Tonga area since 2005, containing their depths but not their magnitude.

### USED_DATA_FILES
contains data files actually used for tests, with various selection of depth/magnitude/area/time, relocalized depth or not...
also contains text files of the beta values computed by the code that will be used to draw maps.

### RESULTS
will contain figures automatically saved by the code.

### PREVIOUS_RESULTS
contain all figures obtained by a previous run of the entire code.

### STAN
contains the bayesian fitting codes in stan language. STAN_Utilities.py and LOO_Utilitie_function.py contain functions used to run the fitting code with the 2 models and plotting the result.
I didn't write STAN_Utilities but I did write LOO_Utilitie_function.py.


### PYTHON FILES
main_execute.py file run all the work. Keywords being true or false are used to select a piece of it.
If you have a MAC OS :
Packages used are all contained in the conda environnement called internship which is in the REQUIREMENTS.yml file. To run the code you thus need to write 'conda activate internship' in your terminal, then open your python and run main_execute.py.
Otherwise :
The dependencies of the conda environment doesn't work for other OS. The packages you need to install to run the code are : matplotlib, numpy, scipy, statistics, obspy, basemap, and pystan. 

main_execute.py is using the following functions and classes, which are also commented and explained within the py files :

- extract.py read raw catalog data and extract it into lists, compute the time in days of events since the start time (2005-01-01) thanks to the julian date, and write the lists in an output file.
- relocalize.py assign unknown depths events to a depth category, draw maps for examples.
- sort_data.py select events from the full data files given a time, magnitude, area, and depth window into a sorted data outfile.
- window_execute.py run all the process described in part 2.2 of the report and called significant change point identification code. It saves the final figure and prints the change point times. 
- window_steps.py do the same but with plotting and printing intermediate results at each step.
- betamap.py draw the map of significance of a given change thanks to the beta value. It saves a file containing the beta values of the grid and a map. This is
- plots.py contains functions plotting figures to test the triggering possibility, and plotting sum up figures of the main results.