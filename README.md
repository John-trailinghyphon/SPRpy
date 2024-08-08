# SPRpy - Analysis of MP-SPR measurements with python and dash

This program can be used to perform data analysis methods on multi-parameter surface plasmon resonance (MP-SPR) measurements
acquired using [Bionavis SPR instruments](https://www.bionavis.com/en/technology/why-choose-mpspr/). Specifically, fresnel modelling and exclusion height analysis of full angular scans is currently available.

Apart from launching SPRpy and selecting files in the file dialog windows, SPRpy is designed to be a fully interactive graphical user interface (GUI) running inside a web browser, thus requiring no programming knowledge to operate (of course, the core functions of the code may be adapted for your own programming scripts). 

Fresnel calculations are based on MATLAB implementations of the [transfer-matrix-method](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics)) 
by [Andreas Dahlin](https://www.adahlin.com/matlab-programs). The GUI elements are built using [Plotly Dash](https://dash.plotly.com/).

## Installation

SPRpy is written in Python 3.11. It is not yet compatible with python > 3.11 (dependency on package [bottleneck](https://pypi.org/project/Bottleneck/) <= python 3.11), and has not been tested on earlier versions of python (the source code is also available to clone from the [SPRpy github repository](https://github.com/John-trailinghyphon/SPRpy) instead). Python 3 releases can be found [here](https://www.python.org/downloads/). NOTE: It is recommended to check the box during installation that adds python.exe to the PATH environment variable (or see manual instructions [here](https://datatofish.com/add-python-to-windows-path/)) to allow running SPRpy scripts by double-clicking them in your file explorer. Alternatively, SPRpy can be set up as any python project in your favourite IDE.

SPRpy is available on [PyPI](https://pypi.org/project/SPRpy/) and can be installed using pip (after installing python):

Windows (in cmd or Powershell):
```python -m pip install SPRpy```

Linux/Mac (always remove "python -m"):
```pip install SPRpy```

To add a shortcut to the SPRpy folder to the desktop after installation, run the following command in the command prompt (windows only):
```SPRpy-desktop-shortcut```

To update to a newer version of SPRpy (overwriting the previously installed version), run the following command:
```python -m pip --no-cache-dir install --upgrade SPRpy```

To install an additional copy of a specific version of SPRpy, run the following command:
```python -m pip install --target ADD_FOLDER_PATH --ignore-installed SPRpy==X.Y.Z```

(change ADD_FOLDER_PATH to desired folder and change X.Y.Z to desired version). This may be necessary to properly open older SPRpy sessions no longer compatible with the latest release.

Note that SPRpy is designed to leverage parallel computing where applicable, thus its performance will be heavily influenced by the number of logical processors of your CPU and their individual clock speeds. While running exclusion height modelling calculations, one can typically expect to see 100 % CPU usage on a 12th generation Intel i7 with 10 logical processors with a runtime of a few minutes. Low-end laptops with weaker CPUs can experience significatly longer computation times in comparison.

## Running SPRpy

Before running SPRpy, you need to convert your MP-SPR measurement files (.spr2) to a specific .csv format. This can be achieved by running two separate scripts (simply double-click):
1) "SPRpy_X_cal.py", a script which generally only needs to be run once to convert the stepper motor values to angles for a particular Bionavis instrument (depending on its setup). Requires a full range scan at highest angular resolution (slow scan), along with its .spr2 file and corresponding exported .dto files from the Bionavis Viewer for each instrument wavelength. The script produces a .csv file that is used by the second script in 2).
2) "SPRpy_spr2_to_csv.py", a script that is used to convert measurements (.spr2) to a specific .csv format used by SPRpy. You will be prompted to select a .spr2 measurement file to convert (and X_cal.csv, file unless the script finds the default). One .csv file will be created for each wavelength in the same folder as the original file with the filename of the original + channel and wavelength information (NOTE! The appended part of the file name must not be changed, it is used by SPRpy).

A text configuration file can be found as "config.toml" that contain some settings that can be tuned. The path in "default_data_folder" can be set to a folder of your choice where you will be initially prompted when loading new data. To run SPRpy, double-click "SPRpy.py" from the SPRpy folder or run it inside a python interpreter.

SPRpy will first prompt you if you wish to load a previous session or start a new one. All sessions are initially created and stored in a subfolder as ...\\SPRpy\\SPRpy sessions\\SESSION EXAMPLE FOLDER. By default, each new session folder is generated with a name containing the date and time of its creation (thus giving it a unique name), but it can also be renamed to whatever you want inside the GUI while SPRpy is running. One can also rename or move the session folder using the file explorer when SPRpy is **not** running. However, its content structure or .pickle file names must not be changed! If you choose to load a previous session, you will be prompted to select a previous session.pickle file from a session folder. If you choose to start a new session, you will instead be prompted to select an initial SPRpy converted .csv measurement data file to load. NOTE! If you open the converted .csv files in a 3rd party program (like excel), it is recommended to **not** save them as the default .csv option as this may break the formatting (if this happens, rerun the SPRpy_spr2_to_csv.py conversion script for that measurement). Additional measurement files can be added later in the GUI workflow.

Next, the GUI will be initiated to a local IP address (by default http://127.0.0.1:8050/). Simply open this address in your browser to access the GUI. It is recommended to add this address as a bookmark for easy access. If you wish ro run multiple instances of SPRpy simultaneously, you can increment the host number in the config.toml file (e.g. http://127.0.0.2:8050/) and add the new IP address in your browser window before running SPRpy again. NOTE: It is a bad idea to open the same session file in two simultaneously running instances of SPRpy...

## Usage

### Workflow with sessions, sensors and analysis instances 

SPRpy is designed to perform the data analysis of multiple measurements within a working *session*. A session starts upon running SPRpy, with the optional choice to reload a previous session. Each session instance contains all the needed parameters to perform the analysis, along with any obtained results, and the session is automatically saved every time a change occurs. Thus, apart from when calculations are still running, SPRpy can be aborted at will and later be relaunched where a previous session may be reloaded to pick up where you left off. 

Within each session, any number of *sensor* instances can be added. A sensor instance contains information corresponding to the type of sensor chip that was used in a SPR measurement via the parameters describing its optical layer structure (layer thicknesses and wavelength dependent refractive index values). By default, an initial gold sensor will be instanced when starting a new SPRpy session. Additional sensor instances may be quick-added as gold, SiO2, platinum or palladium, with refractive index values matching the wavelength of the currently loaded measurement file. However, any custom layer structure may be built by interacting with the layers of the sensor table in the GUI. The grey button can be used to bring up a temporary table of refractive index values (this can even be customized in the config.toml file for your own materials). Refractive index values for various materials can be found at  [refractiveindex.info](https://refractiveindex.info/). 

Separate to the sensor instances are *analysis* instances for each type of analysis method (currently fresnel modelling (FM) and exclusion height determination (EH)). These keep track of selected additional model specfic parameters and the results. When a new analysis instance is added, it draws its data and optical parameters from the currently loaded measurement file and the currently selected sensor. However, when selecting previous analysis instances with mismatching data paths to the currently loaded measurement file (indicated by green data trace instead of blue), rerunning calculations will pull data from the initial path (this can fail if the folders or file along the path has been moved or renamed since).

Sharing sensor or analysis instances between different sessions via their .pickle files is currently not supported.

### Session name and session log

A log field is available to write comments that are saved to the session. The session name can also be renamed here.

### File and sensor controls

Next you can find the file and sensor controls. The file controls are used to load new data files and define optical parameters for the sensor that is used. NOTE! In general, all calculations and changes that are executed will pull the data from the currently loaded measurement file and parameter values from the currently loaded sensor in the sensor table (except during batch processing). 

When adding new sensors the refractive index values will be updated according to the wavelength from the filename of the current measurement file and the refractive index of the bulk medium is always calculated based on the TIR angle of the currently loaded measurement data. Copying a sensor will update the wavelength and channel name based on the currently loaded measurement file, but will **NOT** update the refractive index values if a different wavelength measurement was loaded (i.e. only from copy sensors with the same wavelength). 

New layers can be added to the sensor table with the "Add layer" button and layers can also be removed using the crosses next to the label column. The sensor table values can be directly edited, however, note that "Save edited values" must be clicked to commit any user edited values in the sensor table before they take effect (including changing number of layers). Clicking outside of the table cells while editing a value will abort the editing. Instead, select a neighboring cell when finished typing by pressing enter/tab/arrow keys or clicking to make the table accept what was typed. For fresnel fitting, the main variable to be fitted is selected using the sensor table by highlighting a value (red tint) and clicking on the green button "Select variable to fit" (no need for clicking the red save edited values button afterwards). The sensor may also be renamed to what you wish by clicking "Rename sensor", however it will always have a unique identifier SX, where X will be number of created sensors for this session. It is recommended to include some short unique identifier based on the type of layer that can be easily associated with a measurement file (like Au1, Glass2, PEG3 etc.).

### Response quantification

The first tab amoung the analysis options shows two figures for the loaded measurement file. The first left figure shows the angular trace for the last scan of the measurement file by default. The currently presented trace of this figure corresponds to the data that is used for fresnel modelling. Below it is a button that allows for loading more angular scan traces to the figure from other .csv files. Additionally, a second button can be used to plot theoretical fresnel model traces based on the values present in the currently selected sensor data table.

The right plot shows the full measurement sensorgram, including the SPR and TIR angle traces. When hovering the cursor over particular data points in the sensorgram, the left figure updates with the corresponding angular trace (unless the hover lock toggle is activated). This may be used to select exactly which angular trace to use for fresnel modelling from the loaded measurement. Additionally, clicking a datapoint in the sensorgram will create a new offset in the Y-axis at this timepoint. Clicking the legend of any of the data traces will hide it (this goes for all figures in general).

### Fresnel modelling

The second analysis tab may be used for performing fresnel model fits of angular scan measurements against any of the variables in the sensor data table. The SPRpy implementation of fresnel modelling operates under the following assumptions: 
- Each layer is assumed to be homogenous in X,Y and Z (good assumption for non-particulates and non-swollen films)
- The refractive index of each material layer is assumed to be equivalent to its dry bulk counterpart (good assumption in air)
- The incoming light is monochromatic
- The incoming light has no beam divergence
- The incoming light is fully p-polarized (can be adjusted with the p-factor setting)

Additionally, larger deviations in observed reflectivity and theoretically predicted reflectivity is assumed to be due to non-ideal optical conditions of the instrument and/or sensor glass. Perfectly compensating for these discrepancies in a physically accurate manner by tuning all available parameters within reasonable tolerances can be challenging and time-consuming. Fortunately, the reflectivity on its own generally carries little to no information about the thickness and refractive index of the sensor surface adlayer(s) of interests. Thus, in SPRpy, fitting is focused around the SPR minimum and a constant reflectvity offset of the fresnel coefficients is by default fitted together with the chosen layer variable. Additionally, small corrections to the prism extinction coefficient, *k*, can be simultaneously fitted (also by default) to compensate for peak broadening around the SPR angle minimum with negligible influence on the adlayer properties (works for adlayer materials with k ~ 0). 
#### Experimental procedure

1) Measure the angular reflectivity trace (including the TIR region) of the cleaned sensor (or use a previous representative one, depending on required accuracy).
2) Add the layer of interest or perform desired surface treatment and measure the sample.
#### Modelling procedure

1) Load the desired measurement file and create a new sensor instance using the file and sensor controls. Alternatively, if you remeasure the same sensor after multiple treatments or layers, a previous sensor can be selected and used as is or after adding a new layer to it.
2) Choose which parameter should be fitted. For clean metal sensors, choose the extinction coefficient, _k_, of the metal layer (the plasmonic metal layer thickness should usually not be altered). The thickness of the Cr layer may also need to be manually tuned. For all other types of layers, select the surface layer thickness, *n* or *k*, as variable to be fitted.
3) In the fresnel modelling tab, click "Add new fresnel analysis". 
4) Adjust the initial guess and lower and upper bounds if needed. By default, the current value of the fitted variable is presented as an initial guess, with lower and upper bounds calculated as initial guess / 4 or initial guess * 2 respectively. 
5) Then choose a suitable angle range to fit using the sliders, unless the initial automatic guess is already satisfactory. A suitable range should be covering say 20-60 % of the lower part of the SPR minimum. In the config.toml file one can tune how many points the automatic guess should include above and below the minimum value of the SPR dip.
6) It is recommended (at least initially) to simultaneously fit an intensity offset and the prism *k* value, but this may also be disabled using the two checkboxes. 
7) Finally, press run calculations and wait for the result to show up in the bottom of the settings field and variables being updated in the sensor table.
#### Batch analysis 

For modelling of several replicates with the same sensor layer structure and materials, the batch analysis button is both convenient and time saving. It requires an example sensor and example analysis that has already been run and which parameters will be copied over. There are then two main options to choose from: 1) the example sensor is used directly as a template for new sensor instances for each replicate,  or 2) individual sensor backgrounds are selected as templates for each selected measurement file and a new layer is added according to the surface layer of the analysis example. For option 2), one may also choose between adding the new layer directly to the background sensor instance,  or instead making a new copy for each.

### Exclusion height determination

In the exclusion height determination tab the non-interacting height probe method can be used to determine the exclusion height of a probing particle for a swollen layer in solution. The method is based on the following peer-reviewed papers:

Schoch, R. L. and Lim, R. Y. H. (2013). Non-Interacting Molecules as Innate Structural Probes in Surface Plasmon Resonance.
_Langmuir_, _29(12)_, _4068–4076_. 
https://doi.org/10.1021/la3049289

Emilsson, G., Schoch, R. L., Oertle, P., Xiong, K., Lim, R. Y. H., and Dahlin, A. B. (2017). 
Surface plasmon resonance methodology for monitoring polymerization kinetics and morphology changes of brushes—evaluated with poly(N-isopropylacrylamide). 
_Applied Surface Science_, _396_, _384–392_. 
https://doi.org/10.1016/j.apsusc.2016.10.165

The non-interacting probe method have the following requirements: 
   * The probing particle does not interact with the layer of interest (it should be excluded from entering into it and not stick to it)
   * The layer is sufficiently thin compared to the sensor decay length such that it is able to give a large enough response to the injected probe (use longer wavelengths for thicker films).
   * Swollen material layers with non-zero *k* values are seemingly unsuitable for the non-interacting probe method
#### Experimental procedure

Measure the SPR/TIR response from a surface containing your layer of interest while making 2-3 repeated injections of the probing particles for 5-10 minutes each. High contrast in the response with the probe compared to running buffer is required for accurate determination of the exclusion height. For protein or polymer probes, 10-20 g/L is typically used depending on the expected layer thickness to get sufficient contrast. Verify that the baseline response returns to the same level after rinsing out each probe injection (normally within < 10-20 millidegrees deviation, approximately).

#### Modelling procedure 

1) Response quantification & Fresnel model tabs: It is necessary to first fresnel model a liquid scan from the same height probing measurement file. Use a scan immediately before the first probe injection starts for the modelling. The right scan can be selected by highlighting this part of the sensorgram in the response quantification tab with the hover lock unselected (be careful to move your cursor up and down without accidentally selecting the wrong part of the trace). Add a sensor layer corresponding to your swollen layer of interest, with *n* value of its dry bulk state (the thickness and *n* value doesn't actually matter, but make sure *k* is 0, non-negative *k* for swollen layers will likely not work with this method). Make sure to include the offset and prism k correction fits here, as this is will improve the result from the height exclusion height algorithm. Note that the obtained fitted result will not be physically accurate as it would correspond to a 0 % hydrated state (if that is what you set *n* to be), so its value should be ignored for now.
2) Exclusion height determination tab: Click add new exclusion height analysis. A prompt will pop up asking for the required background fresnel object before proceeding.
3) Check that the height bounds are reasonable for your expected swollen layer. By default, they are calculated according to the 0 % hydrated state from the fresnel background and 6 times this value.
4) Change to "Choose injection points" under the SPR sensorgram and click the SPR data trace before and after each probe injection (so 2 points per injection). These points will be used to plot the SPR angle vs TIR angle traces during probe injections, which may help to verify if the probe is truly non-interacting (linear relationship with low degree of hysteresis means non-interaction).
5) Switch to "Choose buffer points". At this point it may help to click the legend of the TIR angle trace and injection point markers to hide them. A stable range of scans without probe present just before and after each injection should be selected, i.e. a total of 4 selected points per injection. 
6) Next, switch to "Choose probe points". Choose a suitable stable range on top of the probe injection, i.e. a total of 2 points per injection.
7) Once all points are selected, click "Initialize model". For each range of previously selected buffer and probe points all scans within it will be averaged into a single average scan. Scroll down to verify that the SPR vs TIR angle and averaged buffer and probe plots look OK for each injection step. In some cases errors may appear due to how the points were selected, then try clearing the selected points and make a new attempt (try only clicking one of the traces and avoid any markers). 
8) Finally, click "Start calculations". The exclusion height will then be calculated for each "step" in the response: the buffer -> injected probe (step up), and, probe -> buffer rinse (step down). Thus, 2 exclusion heights are calculated for each probe injection, and all of them are also averaged into a single value with a standard deviation. Once the calculations are finished, new plots of possible fitted pairs of thickness and refractive index for both the buffer and probe averaged scans are presented for each step. The exclusion height is found where these curves graphically intersect (this is automatically detected, but good to verify it worked correctly if the values seem odd). If the exclusion height values differ significantly between different steps, there could be problems with the selected points (try again with a new set of points). Problems may also occur if the probe interacts with something on the sample over time, partly adsorbs to the surface, or needs longer time to rinse properly from the flow cell (shift buffer range to further after probe rinsing). Sometimes no intersection occurs for a data set no matter which points are selected, then one has to retry the experiment, and if this still doesn't work deeper investigation are needed, alternatively the swollen layer or probe may not be suitable for the non-interacting probe method. Note that the calculations may take several minutes. If they take way too long, the "Resolution" setting can be lowered to gain some speed if needed (at a loss of accuracy). While generally not needed, the fitting may be further improved across all points of the thickness/refractive index pairs by checking the two "Refitting" options (again at the expense of slightly longer computation times). Remember to press "Initialize model" again before rerunning calculations if any settings has been changed since.


### Result summary (Planned feature, WIP)

Here you can select various results and export and plot them in different ways.

### The dual-wavelength method (Planned feature, WIP) 
The dual-wavelength method can be used to determine the extension of swollen layers with unknown refractive index, 
based on the following requirements: 
   * The refractive index increment for the layer material for each wavelength is known.
   * The thickness of the measured layer is _much smaller_ than the decay length for each wavelength.
   * The sensitivity factors for each wavelength of the instrument is known (easily determined).

It is based on the peer-reviewed paper by:

Rupert, D. L. M., et al. (2016). Dual-Wavelength Surface Plasmon Resonance for Determining the Size and Concentration of Sub-Populations of Extracellular Vesicles. 
_Analytical Chemistry_, _88(20)_, 9980–9988. 
https://pubs.acs.org/doi/full/10.1021/acs.analchem.6b01860


