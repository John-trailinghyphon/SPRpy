# SPRpy: GUI analysis methods for MP-SPR measurements

This program can be used to perform data analysis on multi-parameter surface plasmon resonance (MP-SPR) measurements
acquired using [Bionavis SPR instruments](https://www.bionavis.com/en/technology/why-choose-mpspr/). Specifically, fresnel modelling and exclusion height analysis of full angular scans is currently available.

Apart from launching SPRpy and selecting files in the file dialog windows, SPRpy is designed to be a fully interactive graphical user interface (GUI) running inside a web browser, thus requiring no programming knowledge to operate (of course, the core functions of the code may be adapted for your own programming scripts). 

Fresnel calculations are based on MATLAB implementations of the [transfer-matrix-method](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics)) 
by [Andreas Dahlin](https://www.adahlin.com/matlab-programs). The GUI elements are built using [Plotly Dash](https://dash.plotly.com/).

To cite this work, navigate to the "Cite this repository" menu on the [SPRpy github repository](https://github.com/John-trailinghyphon/SPRpy), or click this Zenodo badge:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13479400.svg)](https://doi.org/10.5281/zenodo.13479400)

## Installation

### Recommended: install with uv

[uv](https://docs.astral.sh/uv/) installs SPRpy in its own isolated environment **and automatically downloads a compatible Python version** (SPRpy works with Python 3.11-3.12), so Python does not have to be installed separately and other Python projects on the computer are not affected.

1) Install uv (once). Windows (PowerShell):
```powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"```
Linux/Mac:
```curl -LsSf https://astral.sh/uv/install.sh | sh```

2) Install SPRpy:
```uv tool install SPRpy```

3) Start SPRpy by typing ```sprpy``` in a terminal (or use the desktop shortcut, see below). The following commands are also installed:
- ```sprpy-convert``` converts .spr2 files to .csv (select files in a dialog, or give them directly: ```sprpy-convert "C:/data/measurement.spr2"```)
- ```sprpy-xcal``` creates the X-calibration file of an instrument (see below)
- ```sprpy-folder``` opens the SPRpy folder (with config.toml, X-calibration files and the "SPRpy sessions" folder)
- ```SPRpy-desktop-shortcut``` (Windows) adds two desktop shortcuts: "SPRpy" starts SPRpy and "SPRpy folder" opens the SPRpy folder

To update SPRpy: ```uv tool upgrade SPRpy```. To install a specific version (e.g. to open older sessions): ```uv tool install SPRpy==X.Y.Z --force```, or run an older version once without installing it: ```uvx --from SPRpy==X.Y.Z sprpy```.

NOTE: Updating replaces the files in the SPRpy folder, including config.toml. Keep a backup copy of your config.toml (with a different filename) and copy your settings back after updating. New settings that are missing in an older config.toml are taken from config.toml.defaults automatically.

### Alternative: install with pip

SPRpy is written in Python 3.11 and also works with 3.12. It is not yet compatible with python > 3.13 (issue with tkinter for python 3.13.0, should be fixed in 3.13.1), and has not been tested on earlier versions of python. Python 3 releases can be found [here](https://www.python.org/downloads/). NOTE: It is recommended to check the box during installation that adds python.exe to the PATH environment variable (or see manual instructions [here](https://datatofish.com/add-python-to-windows-path/)) to allow running SPRpy scripts by double-clicking them in your file explorer. Alternatively, SPRpy can be set up as any python project in your favourite IDE.

SPRpy is available on [PyPI](https://pypi.org/project/SPRpy/) and can be installed using pip (after installing python):

Windows (in cmd or Powershell):
```python -m pip install SPRpy```
(the source code is also available to clone from the [SPRpy github repository](https://github.com/John-trailinghyphon/SPRpy))

Linux/Mac (always remove "python -m"):
```pip install SPRpy```

The commands ```sprpy```, ```sprpy-convert```, ```sprpy-xcal```, ```sprpy-folder``` and ```SPRpy-desktop-shortcut``` described above are installed with pip as well. If the shortcut does not appear, for Windows you can usually find the installation folder by pasting the following in the explorer window (also note to make sure the python version e.g. \Python311\ or \Python312\, matches the one you installed) ```%USERPROFILE%\AppData\Local\Programs\Python\Python311\Lib\site-packages\SPRpy```

To update to a newer version of SPRpy (overwriting the previously installed version), run the following command:
```python -m pip --no-cache-dir install --upgrade SPRpy```

To install an additional copy of a specific version of SPRpy, run the following command:
```python -m pip install --target ADD_FOLDER_PATH --ignore-installed SPRpy==X.Y.Z```

(change ADD_FOLDER_PATH to desired folder and change X.Y.Z to desired version). This may sometimes be necessary to properly open older SPRpy sessions no longer compatible with the latest release.

### Development version (from the github repository)

With uv, the repository can be run directly; uv creates the environment with the right Python version from pyproject.toml:
```
uv run sprpy
uv run pytest
```

Exclusion height calculations take a few seconds per analysis with the vectorised fresnel model and the adaptive height search (they ran in parallel processes and took minutes in SPRpy <= 1.2).

## Running SPRpy

A text configuration file can be found as "config.toml" that contain some default initial settings that can be tuned. Especially make sure to check that the single 'instrument_TIR_sensitivity' value and several values below [instrument_SPR_sensitivity] matches your specific instrument (when downloading it initially matches measured values of the SPR from A. Dahlin's lab, other users should switch to the highlighted defaults from BioNavis or other determined values).

Before running SPRpy, you need to convert your MP-SPR measurement files (.spr2) to a specific .csv format. This can be achieved by running two separate scripts (simply double-click):
1) "SPRpy_X_cal.py", a script which generally only needs to be run once to convert the stepper motor values to angles for a particular Bionavis instrument (depending on its setup). Requires a full range scan at highest or medium angular resolution (slow 1/1 or medium 1/5 scans), along with its .spr2 file and corresponding exported .dto files from the Bionavis Viewer for each instrument wavelength. The script produces a .csv file that is used by the second script in 2). You can update the default calibration .csv that is loaded in config.toml at the entry 'default_poly_file'.
2) "SPRpy_spr2_to_csv.py", a script that is used to convert measurements (.spr2) to a specific .csv format used by SPRpy. You will be prompted to select a .spr2 measurement file to convert (and X_cal.csv, file unless the script finds the default). One .csv file will be created for each wavelength in the same folder as the original file with the filename of the original + channel and wavelength information (NOTE! The appended part of the file name must not be changed, it is used by SPRpy). The conversion is vectorised and takes a few seconds even for long measurements (one .spr2 file of ~170 MB converts in ~7 s instead of ~90 s). Files can also be given on the command line, e.g. `python SPRpy_spr2_to_csv.py "C:/data/measurement.spr2"`. The intensity calibration now interpolates the calibration curve at each scan angle; set CALIBRATION_MODE = 'original' at the top of the script to reproduce .csv files from earlier SPRpy versions exactly.

To run SPRpy, type ```sprpy``` in a terminal, use the "SPRpy" desktop shortcut, or double-click "SPRpy.py" in the SPRpy folder (or run it inside a python interpreter). SPRpy starts with a menu window with three choices:
1) **Open measurement**: select a converted .csv file, or a .spr2 file directly. A .spr2 file is converted first (one .csv file per channel, next to the .spr2 file), after which you choose which channel to open.
2) **Reopen session**: select the session file (.pickle) of a previous session.
3) **Convert .spr2 files only**: convert one or several .spr2 files without starting the analysis. The menu stays open, so you can continue with 1) or 2) afterwards.

Closing the menu (the cross in the corner, or Esc) quits SPRpy. Cancelling a file dialog simply returns to the menu. The menu window then shows the loading progress, and the GUI opens in your web browser automatically when it is ready. Set show_start_menu = false in config.toml to use the console prompts and file dialogs instead (and open_browser = false to not open the browser automatically). The conversion scripts in step 1 and 2 above can likewise be started with ```sprpy-xcal``` and ```sprpy-convert```.

Without the start menu (show_start_menu = false), SPRpy will first prompt you in the console if you wish to load a previous session or start a new one. All sessions are initially created and stored in a subfolder as ...\\SPRpy\\SPRpy sessions\\SESSION EXAMPLE FOLDER. By default, each new session folder is generated with a name containing the date and time of its creation (thus giving it a unique name), but it can also be renamed to whatever you want inside the GUI while SPRpy is running. One can also rename or move the session folder using the file explorer when SPRpy is **not** running. However, its content structure or .pickle file names must not be changed! If you choose to load a previous session, you will be prompted to select a previous session.pickle file from a session folder. If you choose to start a new session, you will instead be prompted to select an initial SPRpy converted .csv measurement data file to load. NOTE! If you open the converted .csv files in a 3rd party program (like excel), it is recommended to **not** save them as the default .csv option as this may break the formatting (if this happens, rerun the SPRpy_spr2_to_csv.py conversion script for that measurement). Additional measurement files can be added later in the GUI workflow.

Next, the GUI will be initiated to a local IP address (by default http://127.0.0.1:8050/). Simply open this address in your browser to access the GUI. It is recommended to add this address as a bookmark for easy access. If you wish ro run multiple instances of SPRpy simultaneously, you can increment the host number in the config.toml file (e.g. http://127.0.0.2:8050/) and add the new IP address in your browser window before running SPRpy again. NOTE: It is a bad idea to open the same session file in two simultaneously running instances of SPRpy...

## Usage

### Workflow with sessions, sensors and analysis instances 

SPRpy is designed to perform the data analysis of multiple measurements within a working *session*. A session starts upon running SPRpy, with the optional choice to reload a previous session. Each session instance contains all the needed parameters to perform the analysis, along with any obtained results, and the session is automatically saved every time a change occurs. Thus, apart from when calculations are still running, SPRpy can be aborted at will and later be relaunched where a previous session may be reloaded to pick up where you left off. 

Within each session, any number of *sensor* instances can be added. A sensor instance contains information corresponding to the type of sensor chip that was used in a SPR measurement via the parameters describing its optical layer structure (layer thicknesses and wavelength dependent refractive index values). By default, an initial gold sensor will be instanced when starting a new SPRpy session. Additional sensor instances may be quick-added as gold, SiO2, platinum or palladium, with refractive index values matching the wavelength of the currently loaded measurement file. However, any custom layer structure may be built by interacting with the layers of the sensor table in the GUI. The grey button can be used to bring up a temporary table of refractive index values (this can even be customized in the config.toml file for your own materials). Refractive index values for various materials can be found at  [refractiveindex.info](https://refractiveindex.info/). 

Separate to the sensor instances are *analysis* instances for each type of analysis method (currently fresnel modelling (FM) and exclusion height determination (EH)). These keep track of selected additional model specfic parameters and the results. When a new analysis instance is added, it draws its data and optical parameters from the currently loaded measurement file and the currently selected sensor. However, when selecting previous analysis instances with mismatching data paths to the currently loaded measurement file (indicated by green data trace instead of blue), rerunning calculations will pull data from the initial path (this can fail if the folders or file along the path has been moved or renamed since).

Sharing sensor or analysis instances between different sessions via their .pickle files is currently not supported.

### Session name and session log

A log field is available to write comments that are saved to the session. The session name can also be renamed here. SPRpy also writes the important steps of the analysis to the log automatically (with time stamps), so that the analysis can be traced afterwards: loaded measurement files (number of scans, scan speed and TIR range), added, copied, edited, renamed and removed sensors, TIR/SPR fitting settings, drift corrections, surface coverage results, fresnel modelling results (fitted value, offset, prism k, angle range and data file), batch analysis results, exclusion height model settings and results, and exported files. The messages are also printed in the console window.

### File and sensor controls

Next you can find the file and sensor controls. The file controls are used to load new data files and define optical parameters for the sensor that is used. NOTE! In general, all calculations and changes that are executed will pull the data from the currently loaded measurement file and parameter values from the currently loaded sensor in the sensor table (except during batch processing). 

When adding new sensors the refractive index values will be updated according to the wavelength from the filename of the current measurement file and the refractive index of the bulk medium is always calculated based on the TIR angle of the currently loaded measurement data. Copying a sensor will update the wavelength and channel name based on the currently loaded measurement file, but will **NOT** update the refractive index values if a different wavelength measurement was loaded (i.e. only from copy sensors with the same wavelength). 

New layers can be added to the sensor table with the "Add layer" button and layers can also be removed using the crosses next to the label column. The sensor table values can be directly edited, however, note that "Save edited values" must be clicked to commit any user edited values in the sensor table before they take effect (including changing number of layers). Clicking outside of the table cells while editing a value will abort the editing. Instead, select a neighboring cell when finished typing by pressing enter/tab/arrow keys or clicking to make the table accept what was typed. For fresnel fitting, the main variable to be fitted is selected using the sensor table by highlighting a value (red tint) and clicking on the green button "Select variable to fit" (no need for clicking the red save edited values button afterwards). The sensor may also be renamed to what you wish by clicking "Rename sensor", however it will always have a unique identifier SX, where X will be number of created sensors for this session. It is recommended to include some short unique identifier based on the type of layer that can be easily associated with a measurement file (like Au1, Glass2, PEG3 etc.).

### Response quantification

The first tab amoung the analysis options shows two figures for the loaded measurement file. The first left figure shows the angular trace for the last scan of the measurement file by default. The currently presented trace of this figure corresponds to the data that is used for fresnel modelling. Below it is a button that allows for loading more angular scan traces to the figure from other .csv files. Additionally, a second button can be used to plot theoretical fresnel model traces based on the values present in the currently selected sensor data table.

The right plot shows the full measurement sensorgram, including the SPR, TIR and bulk corrected angle traces. When hovering the cursor over particular data points in the sensorgram, the left figure updates with the corresponding angular trace (unless the "Stop mouse hover updates" toggle is activated). This may be used to select exactly which angular trace to use for fresnel modelling from the loaded measurement. Additionally, clicking a datapoint in the sensorgram will create a new offset in the Y-axis at this timepoint. Clicking the legend of any of the data traces will hide it (this goes for all figures in general). The bulk corrected trace is calculated according to the formula presented at the bottom of the page, where each parameter may be adjusted (defaults can be changed in config.toml). The toggle "Show TIR/SPR fitting parameters" can be clicked to show two additional plot windows that are updated when hovering over points in the sensorgram figure, along with additional parameters for tuning the TIR and SPR angle fitting. NOTE! Tuning these fitting parameters will be applied and saved to the current session as a whole (i.e. also affecting further modelling runs). The default loaded TIR and SPR fitting parameters for newly created sessions can be changed in the file "config.toml".

Two methods are available for each angle (selected in the fitting parameters panel):
- TIR angle, *Savitzky-Golay derivative* (default for new sessions): the derivative of the reflectivity is calculated with a Savitzky-Golay filter whose width is set in degrees (so it adapts automatically to the scan speed), and a parabola is fitted to the top of the derivative peak. This gives a considerably less noisy and more robust TIR angle than the *legacy* moving average method of SPRpy <= 1.2.
- SPR angle, *adaptive dip window* (default for new sessions): a cubic is fitted to the points below a set fraction (the "dip level", default 0.2) of the dip depth, so the fit window adapts to the dip width of each wavelength (broad at 670 nm, narrow at 980 nm). The *legacy* method fits a fixed number of points around the minimum. If no complete dip is found (e.g. Pd or Pt films in air), the lowest point is used instead.

Sessions created with older SPRpy versions keep the legacy methods so that previous results are reproduced, but can be switched to the new ones. When loading a measurement, the TIR range (liquid or air, see config.toml) is chosen from where the TIR edge actually is, so single scans in liquid are handled correctly. Spectra where an angle can not be determined get the value 0 in the sensorgram (instead of being left out), so that they can still be selected by hovering.

*SPR drift correction:* a linear drift of the SPR angle can be removed using the controls below the sensorgram. Select the box select tool in the sensorgram toolbar and drag a box over a stable baseline (e.g. buffer before the injections), or type the start and end time, then click "Correct drift". The fitted slope is removed relative to the first time point (the TIR angle is not corrected), the baseline is shaded in the figure and a note is added to the session log. "Reset" removes the correction. The corrected trace is also used when exporting the sensorgram as .csv.

*Surface coverage:* the surface coverage (ng/cm²) of a thin adsorbed layer can be calculated from the SPR shift between two points of the sensorgram, shown at the bottom of the page. If the layer thickness d of the bulk correction parameters is 0, the thin film approximation m_S = δ·Δθ_SPR / (2·(∂n/∂C)·S_SPR) is used, otherwise m_S = d·Δθ_SPR / ((∂n/∂C)·S_SPR·(1 − exp(−2d/δ))) with the estimated layer thickness d (the two agree for thin layers). The equation that is used is highlighted and marked, and stated with the result. Turn on "Select surface coverage points" below the sensorgram and click two points on the SPR or bulk corrected trace (e.g. on the baseline before an injection and after rinsing). While the toggle is on, clicks select points instead of setting the zero point of the sensorgram; a third click starts a new pair. The shift is taken from the trace of the first point (including the drift correction and the current bulk correction parameters), the decay length δ and the SPR sensitivity S_SPR are taken from the bulk correction parameters, and the refractive index increment ∂n/∂C (cm³/g, about 0.18 for proteins) is given next to the equation (the default can be changed with dn_dC in config.toml). The result is shown below the equation and written to the session log.

The analysis methods for the TIR and SPR angles, the vectorised Fresnel calculation and the drift correction are adapted from MATLAB scripts by Jesper Medin (Chalmers).

Further reading for the bulk correction method: 

Accurate Correction of the “Bulk Response” in Surface Plasmon Resonance Sensing Provides New Insights on Interactions Involving Lysozyme and Poly(ethylene glycol)
Justas Svirelis, John Andersson, Anna Stradner, and Andreas Dahlin
ACS Sensors 2022 7 (4), 1175-1182
https://doi.org/10.1021/acssensors.2c00273

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
6) It is recommended (at least initially) to simultaneously fit an intensity offset and the prism *k* value, but this may also be disabled using the two checkboxes. By default, the refractive index of the bulk medium is calculated from the TIR angle of the measurement every time the model is run. Uncheck "Bulk RI from TIR angle?" to instead use the bulk refractive index entered in the sensor table (e.g. if the TIR angle can not be determined reliably). 
7) Finally, press run calculations and wait for the result to show up in the bottom of the settings field and variables being updated in the sensor table.
#### Layer series

For repeated samples that are measured after each added layer (e.g. clean SiO2 → APS → SMCC → PEG, dry or in liquid), the "Layer series" button fits all samples, steps and channels (wavelengths) at once:

1) Prepare a base sensor in the sensor table for the first step (e.g. "Add new sensor" → Glass for Ti/Au/SiO2), and select the variable to fit in the first step (e.g. SiO2|d).
2) Click "Layer series" and choose the measurement files (all channels of all samples and steps). The sample and step of each file are read from the file name, `<date>-<time>_<sample>_<step>-L<channel>_<wavelength>nm.csv` (as named by the instrument and the converter, e.g. 20260116-135311_SiO2_1_APS-L3_670nm.csv is sample SiO2_1, step APS), and the steps are ordered by measurement time. The detected samples, steps and channels are listed, including missing measurements.
3) Fill in the layer added in each following step (name, n, k and an initial thickness guess). Layers listed in [default_sensor_values] in config.toml (e.g. PEG or PHEAA) get the n and k of each wavelength from there, other layers use the given values.
4) Choose the fit options (reflectivity offset, prism k, bulk RI from the TIR angle) and click "Run layer series".

For each sample and channel, the first step is fitted with a copy of the base sensor with the optical constants of the channel's wavelength, and each following step adds the layer of that step and fits its thickness, with all earlier layers fixed at the results of the same sample and channel (a missing step is skipped, and the next step builds on the previous one). One sensor is added to the session for each sample and channel (e.g. "SiO2_1 L3 670nm"), and each step adds its layer to it (as when adding layers manually), so it ends up with the full layer stack and all fitted thicknesses. Each fit is added as a fresnel analysis (e.g. "SiO2_1 APS L3 670nm"), so it can be selected and inspected like any other analysis (note that the analyses of earlier steps then refer to the sensor with all layers: their fitted results and plots are kept, but running them again uses the full layer stack). The results are shown as the mean ± standard deviation over the samples for each step and channel, and as a table of all fits; fits that end at a bound or describe the dip poorly (rms above 0.01) are flagged. "Export results (.csv)" saves all fits (value, offset, prism k, bulk RI, rms, file) and the summary. For the example data (2 samples × 4 steps × 4 channels) the 32 fits take a few seconds.

#### Batch analysis 

For modelling of several replicates with the same sensor layer structure and materials at one wavelength, the batch analysis button can be used (the layer series above is more automatic, and handles all wavelengths and successive layers). It requires an example sensor and example analysis that has already been run and which parameters will be copied over. There are then two main options to choose from: 1) the example sensor is used directly as a template for new sensor instances for each replicate,  or 2) individual sensor backgrounds are selected as templates for each selected measurement file and a new layer is added according to the surface layer of the analysis example. For option 2), one may also choose between adding the new layer directly to the background sensor instance,  or instead making a new copy for each.

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
8) Finally, click "Start calculations". The exclusion height will then be calculated for each "step" in the response: the buffer -> injected probe (step up), and, probe -> buffer rinse (step down). Thus, 2 exclusion heights are calculated for each probe injection, and all of them are also averaged into a single value with a standard deviation. Once the calculations are finished, new plots of possible fitted pairs of thickness and refractive index for both the buffer and probe averaged scans are presented for each step. The exclusion height is found where these curves graphically intersect (this is automatically detected, but good to verify it worked correctly if the values seem odd). If the exclusion height values differ significantly between different steps, there could be problems with the selected points (try again with a new set of points). Problems may also occur if the probe interacts with something on the sample over time, partly adsorbs to the surface, or needs longer time to rinse properly from the flow cell (shift buffer range to further after probe rinsing). Sometimes no intersection occurs for a data set no matter which points are selected, then one has to retry the experiment, and if this still doesn't work deeper investigation are needed, alternatively the swollen layer or probe may not be suitable for the non-interacting probe method. The calculations take a few seconds (with the adaptive height search, the resolution only sets the step of the first sweep and does not affect the accuracy). While generally not needed, the fitting may be further improved across all points of the thickness/refractive index pairs by checking the two "Refitting" options (again at the expense of slightly longer computation times). Remember to press "Initialize model" again before rerunning calculations if any settings has been changed since.

With "Height search: Adaptive" (default), the search is done in two passes, as in SPR_WetHeight.m: a coarse sweep of 31 heights between the height bounds finds where the buffer and probe curves cross, after which the curves are calculated with "Resolution" points (30-100 recommended) in a narrow window around that crossing (the coarse step containing it plus half a step on each side), and the crossing is located within 0.001 nm with a root finder. The d-n pair plots show the curves of this fine pass, so they are well resolved around the crossing also for a wide height range: the default upper bound is 300 nm (or 12 times the background layer thickness if larger), which also covers thick brushes. Heights below 0.5 nm and heights where a fitted refractive index lies at the bounds are skipped. "Fixed grid" calculates all "Resolution" height steps between the bounds and interpolates the crossing between them, as in SPRpy <= 1.2 (use a narrow height range with this option, since a coarse grid makes the interpolated crossing less accurate).

The exclusion height depends strongly on the bulk refractive index change between buffer and probe (calculated from the TIR angle shift), but not on the absolute bulk refractive index. Typically, an error of 1 millidegree in the TIR shift changes the exclusion height by about 1 nm, which is shown together with the results. By default ("Bulk RI change from: TIR edge overlay"), the TIR shift is found by overlaying the TIR edges of the buffer and probe spectra, which gives the same result regardless of the TIR angle settings (method, smoothing). The alternative "Separate TIR angles" uses the difference between the TIR angles of buffer and probe (as in SPRpy <= 1.2), which can differ by a few millidegrees between TIR methods since the TIR edge changes shape slightly with the probe. The TIR shifts are shown with the results and written to the session log.


### All wavelengths (channels)

A .spr2 measurement contains several channels (e.g. L1 670 nm, L2 980 nm, L3 670 nm and L4 785 nm), which are converted to one .csv file each (name-L1_670nm.csv etc.) and share the same scan times. Fresnel modelling and exclusion height determination made for one channel can be repeated for all channels with the "All wavelengths" panels below the fresnel modelling and exclusion height settings (the channels of the measurement are listed as soon as an analysis is selected):

1) Fresnel modelling: first fit the model for one channel as usual. "Fit all wavelengths" then fits the same scan of the other selected channels, with the same layers, fitted variable and fit options. The angle range of each channel covers the same part of its SPR dip as the range of the first channel: the distances from the SPR minimum are scaled with the width of each dip (which differs a lot between wavelengths), with at least half as many data points. The thickness bounds are widened (up to 10 times the first result), and the table flags channels whose fit ends at a bound or needs a large intensity offset: the model then does not describe the dip well, typically because the optical constants for that wavelength in config.toml (gold and prism in particular) do not match the sensor and instrument, and the fit compensates with a wrong layer thickness, which also affects the exclusion heights. Layer thicknesses are kept, while n and k of each layer are taken from config.toml for the wavelength of each channel (layers that are not listed in config.toml, e.g. a custom layer, keep their values; add a column for the layer to [default_sensor_values] in config.toml to give it a refractive index for each wavelength, e.g. from ellipsometry). A sensor and a fresnel analysis are added for each channel (and reused when fitting again), and the results of all channels are shown in a table, together with notes on the optical constants used. Make sure that the sensor of the first channel uses the same set of values as config.toml (e.g. the prism refractive index), otherwise the channels are modelled with different sets of optical constants.
2) Exclusion height determination: select the points and settings for one channel as usual (its fresnel background must have been fitted for all wavelengths). "Run all wavelengths" runs the analysis for all selected channels with the same points and settings, using the fresnel background of each channel. The table shows the exclusion height of each channel, how much it changes per millidegree of TIR shift (lower is more robust), and the TIR shifts. The results are compared per measurement spot and written to the log. Each channel analysis can be selected afterwards to inspect its curves.

The channels are grouped by measurement spot, set with channel_spots in config.toml (default [[1, 2], [3, 4]]: L1 and L2 measure one spot of the sample, L3 and L4 another; change it for instruments with another layout). Channels of the same spot should ideally agree, but the wavelengths have different sensitivity: 670 nm is generally most accurate for thin layers, and the longer wavelengths for thick layers, so some difference between wavelengths is not unreasonable. Differences between spots can indicate an inhomogeneous sample. The results are therefore compared per spot rather than averaged over all channels. The TIR shift of a probe injection only depends on the bulk solution, so it should be (nearly) the same in all channels, which is a useful check. Note also that results at other wavelengths depend on the optical constants for those wavelengths in config.toml (gold and prism in particular).

### Result summary

The results from fresnel modelling and exclusion height determination are presented in two tables in the left-most column and can be exported into a .csv format using the button at the top.

The right barplot groups and plots fresnel model results based on what sensor object was used as background and all its layers that were fitted.
NOTE: For session files < v1.0.0 only the label for the latest fitted layer is available, but the value will correctly represent earlier layers (you can verify with the analysis column in the table to the left).

### Automated tests

The folder "tests" contains an automated test suite (pytest) that checks critical calculations: that the vectorised Fresnel model reproduces the original implementation, that the TIR and SPR angle methods find the correct angles of modelled spectra, and an end-to-end run on the example measurement (conversion, loading and sensorgram). Run it from the repository folder before publishing a new version:

```
uv run pytest                    # with uv (installs pytest automatically)
python -m pip install pytest     # or with pip...
python -m pytest                 # ...all tests (~40 s)
python -m pytest -m "not slow"   # skip the example data tests
```

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


