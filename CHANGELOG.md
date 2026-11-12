# Changelog

## v2.0.1

### Features
- Fits match the SPR dip angle ("Match dip angle?" option, default for new fresnel and exclusion height analyses, the layer series and the batch analysis): the model dip is placed at the measured SPR angle, and the other fitted values (offset, prism k) describe the rest of the dip. Before, fitting the reflectivity over the angle range could leave the model dip several millidegrees from the measured one depending on the range (about +9 mdeg on average with the default range in the SiO2/APS/SMCC/PEG test data); now the result hardly depends on the angle range. The fit result and the log show the difference between the model dip and the measured SPR angle. Analyses from older sessions keep the reflectivity-only fit until the option is checked.

### Fixes
- The package no longer includes converted example data (.csv files made from the example .spr2 files), which made it too large for PyPI.
- README: what to do when the sprpy command is not found, uv cannot hardlink files (e.g. with OneDrive), or an older version was installed right after a release.
- Corrected the name of Jesper Medin in the documentation and source comments.

## v2.0.0

### Features
- Much faster calculations: the Fresnel transfer-matrix model is now vectorised over all angles (identical results, ~90 times faster). This speeds up fresnel modelling, batch analysis and exclusion height determination.
- Much faster conversion of .spr2 files (SPRpy_spr2_to_csv.py), with a more accurate intensity calibration (interpolation of the calibration curve). Files can also be given on the command line.
- New TIR angle method (Savitzky-Golay derivative with a parabola fit of the derivative peak) and SPR angle method (cubic fit of an adaptive dip window), adapted from MATLAB scripts by Jesper Medin. Considerably less noisy TIR angles. Selectable in the TIR/SPR fitting parameters panel; old sessions keep the legacy methods by default.
- The TIR range (liquid or air) is chosen from where the TIR edge is in the data when loading a measurement, instead of from the number of scans (single scans in liquid now work).
- SPR drift correction in the response quantification tab (linear fit to a selected baseline).
- Checkbox in fresnel modelling to keep the bulk refractive index of the sensor table instead of calculating it from the TIR angle.
- The TIR and SPR fit plots now show the actual fit windows, and are calculated on demand (lower memory use for long measurements).
- Restyled GUI: top bar, collapsible illustrations, cards for session log, measurement and sensor controls and analysis tabs.
- Missing entries in config.toml (e.g. after upgrading) are taken from config.toml.defaults and listed at startup.
- Empty columns/rows in .csv files (e.g. from Excel) are ignored when loading data.
- Automated test suite (folder "tests", run with pytest).
- Installation with uv (`uv tool install SPRpy`), which also installs a compatible Python version automatically. New commands: `sprpy` (start SPRpy), `sprpy-convert` (.spr2 to .csv), `sprpy-xcal` (X-calibration) and `sprpy-folder` (open the SPRpy folder). The desktop shortcut now starts SPRpy, with a second shortcut to the SPRpy folder.
- Start menu with an SPR curve (calculated with the Fresnel model), sharp at any Windows display scaling and always on top: open a measurement (.csv, or .spr2 which is converted first), reopen a session, or only convert .spr2 files. The window then shows the loading progress, and the GUI opens in the browser automatically. Can be turned off in config.toml (show_start_menu, open_browser).
- Many more automatic session log messages (data loading, sensors, fit settings, fresnel and exclusion height results, batch results, exports); the log updates and scrolls automatically.

- Exclusion height: the bulk refractive index change between buffer and probe is now determined by overlaying their TIR edges ("Bulk RI change from" option). The exclusion height depends strongly on this change (about 1 nm per millidegree of TIR shift), and the TIR edge changes shape slightly with the probe, so separately determined TIR angles gave results that depended on the TIR method and its settings (e.g. 7.4 nm with the legacy TIR method and 4.2 nm with the Savitzky-Golay method for the example data, compared with 6.8 nm for both with the overlay). The old behaviour is available as "Separate TIR angles", and is kept for analyses from older sessions until they are initialized again.
- Exclusion height results show the TIR shift of each step and how much the exclusion height changes per millidegree error in the TIR shift.

- All wavelengths: fresnel modelling and exclusion height determination can be repeated for all channels (wavelengths) of a measurement in the "All wavelengths" panels of the fresnel modelling and exclusion height tabs. Each channel gets its own sensor (optical constants of its wavelength from config.toml), fresnel analysis and exclusion height analysis, and the results are compared per measurement spot (channels measuring the same spot are set with channel_spots in config.toml, default L1+L2 and L3+L4).

- Exclusion height: two-pass adaptive height search (default, "Height search" option), as in SPR_WetHeight.m: a coarse sweep finds the crossing of the buffer and probe curves, the curves are then calculated with "Resolution" points in a narrow window around it, and the crossing is located within 0.001 nm with a root finder. The default height range is now up to 300 nm (covering thick brushes).
- All wavelengths: the angle range of each channel's fresnel background is scaled with the width of its SPR dip, the thickness bounds are wider than in the first channel (the 980 nm fit could end at the bound before), and fits at a bound or with a large offset are flagged.
- Exclusion height calculations run in the SPRpy process instead of parallel processes, which only added overhead with the vectorised fresnel model: a full calculation on the example data takes under a second with the adaptive search, and about 9 s with the fixed grid (about 45 s before).

- Surface coverage in the response quantification tab: select two points in the sensorgram ("Select surface coverage points" toggle) to calculate the surface coverage m_S in ng/cm², with the thin film approximation m_S = δ·Δθ_SPR / (2·(∂n/∂C)·S_SPR) for d = 0, or m_S = d·Δθ_SPR / ((∂n/∂C)·S_SPR·(1 − exp(−2d/δ))) for an estimated layer thickness d. ∂n/∂C is given next to the equations at the bottom of the page, and d, δ and S_SPR are taken from the bulk correction parameters.

- PHEAA added to the default sensor values in config.toml, with its refractive index at each wavelength from spectroscopic ellipsometry (1.464, 1.457, 1.454 and 1.449 at 670, 785, 850 and 980 nm), used when fresnel modelling all wavelengths.

- Layer series (fresnel modelling tab): fits repeated samples measured after each added layer (e.g. clean SiO2, APS, SMCC, PEG) for all channels at once. Samples and steps are read from the file names, each step builds on the fitted layers of the same sample and channel with the optical constants of each wavelength, one sensor per sample and channel gets a layer added for each step, every fit is added to the session as a fresnel analysis, and the results (mean ± std over samples per step and channel, all fits, flagged fits) can be exported to .csv.

### Fixes
- Spectra where no TIR or SPR angle is found get the value 0 instead of NaN, so they can still be selected in the sensorgram (e.g. Pd or Pt films).
- Automatic angle ranges no longer fail when the SPR minimum is close to the end of the scan (Pd films), and ignore low reflectivity below the TIR angle.
- Points in exclusion height determination can be selected in any order within each pair.
- Exclusion height analysis objects are saved as EH... instead of FM... by "save all".
- Batch analysis file dialog opens in the folder of the current measurement file.
- k values at 850 nm for PEG, protein and DNA in the default sensor values were set to their n values, now 0.
- Compatibility fixes when loading old sessions (missing fitted layer attributes, fitted_result = None).
- Removed an unused Windows-only import.
- SPRpy_X_cal.py works with the new .spr2 converter.
- Cancelling a file dialog no longer gives errors: at startup SPRpy closes with a message (or returns to the start menu), and in the GUI the action is skipped (load data, add data trace, batch files, save/export).
- The desktop shortcut script was registered twice (as console and GUI script).


## v1.2.1

### Features

### Fixes
- Updated default thickness for SiO2 layers to be 20 nm

## v1.2.0

### Features
- Can now tune initial guess and lower and upper bounds for the fitting of refractive index during exclusion height determination

### Fixes
- Fixed rare problem with previous default x0, lower/upper bounds for RI fitting during exclusion height determination of extremely hydrated layers 
- Increased the default upper height limit in exclusion height determination to a factor of 12 times the lower limit (instead of 6).
- Fixed problem with selecting different exclusion height objects that had not yet run any calculations
- Fixed problem with removing exclusion height objects

## v1.1.9

### Features
- Add data trace now also defaults to current data file's folder location

## v1.1.8

### Features
- Adding Ti as optional adhesion material in default values
- Ti is now default adhesion layer for SiO2 films

## v1.1.7

### Fixes
- Updated default TIR range to be broader and extended to lower angles
- Added missing dependency for version 1.1.6 (yanked).

## v1.1.6

### Features
- When loading new data files the explorer window now defaults to the previous data file's location.

## v1.1.5

### Fixes
- Fixed default k value to 5.42 for Au 850 nm
- Fixed dependency issue with pandas > 3.0, limiting to version 2.3.1

## v1.1.4

### Fixes
- Fixed default k value to 0 for 850 nm wavelength sensors
- Tweaked default k values for Au 
- Fixed crash in SPRpy_X_cal.py, now correctly reads "," as decimal points 

## v1.1.3

### Fixes
- Fixed slow loading times when selecting points in exclusion height determination
- When loading new measurement data the wavelength dependent bulk correction parameters are now also updated

## v1.1.2

### Fixes
- Now correctly updates sensorgram graph after applying new TIR and SPR fitting parameters

## v1.1.1

### Fixes
- Now correctly handles selecting TIR fitting parameters according to scanspeed when loading new measurements in one session.

## v1.1.0

### Features
- Optical parameters for sensor tables are now loaded from values entered in "[default_sensor_values]" in config.toml
- Can now change TIR and SPR angle fitting parameters and refit them using inputs and a button in a collapsable hidden menu in the Response quantification tab. NOTE: This will change the TIR and SPR fitting for hte whole session.
- Default TIR and SPR fitting parameters added to config.toml
- Displays the TIR and SPR angle fits in the reflectivity window when hovering over data points in the SPR sensorgram
- A custom folder location for the "SPRpy sessions" folder and for the data prompt folder can now be set in config.toml 

### Fixes
- Now includes time point information when TIR or SPR angle is not found.
- Angle range markers in Fresnel modelling figure now updates correctly

## v1.0.1

### Fixes
- Add data trace and fresnel calculation now correctly updates the reflectivity figure

## v1.0.0

### Features
- Result summary tab with barplot for fresnel data and data export as .csv now available 
- Included .csv data saving option for all figures and changed to filedialog window for choosing filenames
- Included bulk correction to sensorgram in quantification tab

### Fixes
- Included pywin32 as dependency for Windows platforms (required for command "SPRpy-desktop-shortcut")
- Limited dependency for dash to v2.18.2 as graph changes in dash >= 3.0.0 cause CSS problems with indefinitely growing plots
- Limited dependency for dash-bootstrap-components to v1.7.2 as last compatible release to dash v2.18.2
- Fixed sensor remove logic similarly to removal of fresnel objects
- Added compatibility fix for loading older session versions. NOTE: Will still cause erroneous results in result summary tab if multiple layers were fitted for one sensor object
- Clarified text for batch fresnel modelling options

## v0.2.3

### Fixes
- Fixed build problem with python 3.12

## v0.2.2

### Fixes
- Added support for python 3.12 (still problem with python 3.13.0 due to tkinter issue until patch 3.13.1)
- Added support for latest kaleido version (v0.4.1)

## v0.2.1

### Features
- Added support for using fewer than the maximum number of available logical cores in parallel computing steps. Change "max_logical_cores" in config.toml

### Fixes

- Fixed X-calibration not accounting for scanspeed
- Fixed Y-calibration error for calibration scans with reflectivity above 1.0 
- Moved default X-cal filename to config.toml
- Generalised TIR_determination() to work for any scanspeed (different filtering and fitting for scanspeeds below or above 5)

## v0.2.0

### Features

- Generalised `SPRpy_X_cal.py` and `SPRpy_spr2_to_csv.py` to work with Bionavis instruments with any wavelength setup.
- Added support for 850 nm lasers
- Added batch analysis processing for fresnel modelling
- Added offset and prism extinction fit settings to fresnel modelling and exclusion height determination
- Added automatic angle range detection based on SPR minimum to fresnel modelling (can be tuned in `config.toml`)

### Fixes

- Updated README with accurate documentation
- Updated algorithm in `SPRpy_spr2_to_csv.py` producing more accurate absolute angles (5 mdeg variation between measurements)
- `SPRpy` is now compatible with `numpy 2.0`
- New refractive index values for Cr, Pt
- Changed fresnel algorithm from iterating reflectivity offsets to fitting a reflectivity offset for the whole range

## v0.1.3

### Features

-

### Fixes

- Fixed bug with angle range slider max values
- Fixed bug with renaming sessions

## v0.1.1 

- First working release of `SPRpy`