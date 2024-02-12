# SPRpy - Analysis of MP-SPR measurements with python and dash

This program can be used to perform fresnel modelling and exclusion height analysis of multi-parameter surface plasmon resonance (MP-SPR) measurements
acquired using [Bionavis SPR instruments](https://www.bionavis.com/en/technology/why-choose-mpspr/). 

Fresnel calculations are based on MATLAB implementations of the [transfer-matrix-method](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics)) 
by [Andreas Dahlin](https://www.adahlin.com/matlab-programs). The GUI elements are built using [Dash](https://dash.plotly.com/).

## Installation

Version 0.1 of SPRpy is written in Python 3.11. As of 2024-01-12 it is not yet compatible with python 3.12 (dependency on package [bottleneck](https://pypi.org/project/Bottleneck/) <= python 3.11), and has not yet been tested on earlier versions (in that case clone it from the [SPRpy github repository](https://github.com/John-trailinghyphon/SPRpy) instead). Python 3 releases can be found [here](https://www.python.org/downloads/). NOTE: It is recommended to check the box during installation that adds python.exe to the PATH environment variable (or see manual instructions [here](https://datatofish.com/add-python-to-windows-path/)) to avoid the need to set up your own environment.

SPRpy is available on [PyPI](https://pypi.org/project/SPRpy/) and can be installed using pip:

Windows:
```python -m pip install SPRpy```

Linux/Mac (always remove "python -m"):
```pip install SPRpy```

To add a shortcut to the SPRpy folder to the desktop after installation, run the following command in the command prompt (windows only):
```SPRpy-desktop-shortcut```

To install a newer version of SPRpy`and uninstall the previously installed version, run the following command:
```python -m pip --no-cache-dir install --upgrade SPRpy```

To install a specific version (change X.Y.Z) of SPRpy, run the following command:
```python -m pip install SPRpy==X.Y.Z```

## Usage

### Running SPRpy

A text configuration file can be found "config.toml". The path in "default_data_folder" can be set to a folder of your choice where you will be initially prompted when loading new data.

Before running SPRpy, you need to convert your MP-SPR measurement files (.spr2) to a particular .csv format. This can be done using the "spr2_to_csv.py" script found in the SPRpy folder. To run it, simply double-click the script or run it from a python interpreter. You will be prompted to select a .spr2 file to convert. The converted file will be saved in the same folder as the original file with a similar name but with the laser wavelength and .csv extension.

To run SPRpy, run "SPRpy.py" from the SPRpy folder (double-click) or inside a python interpreter.

SPRpy will first prompt you if you wish to load a previous session or start a new one. All sessions are by default stored in a folder under "SPRpy\SPRpy sessions". The session folder can also be moved and renamed when SPRpy is not running, but its contents must not be changed. If you choose to load a previous session, you will be prompted to select a session file. If you choose to start a new session, you will instead be prompted to select a .csv data file to load.

Next, the GUI will be initiated to a local IP address (by default http://127.0.0.1:8050/). Simply open this address in your browser to access the GUI. It is recommended to add this address as a bookmark for easy access. If you wish ro run multiple instances of SPRpy, you can increment the host number in the config file and add the new IP address in your browser window. NOTE: It is a bad idea to open the same session in two simultaneously running instances...

### Session name and session log

A log field is available to write comments that are saved to the session or rename the session folder. The session name will also be used as a title in the GUI.

### File and sensor controls

Next you can find the file and sensor controls. The file controls are used to load new data files and define optical parameters for the sensor that is used. All calculations that are performed will pull the values from the currently loaded sensor in the sensor table. Note that "Save edited values" must be clicked to commit any user edited values in the table before they take effect.

### Fresnel modelling of non-swollen layers

The basic fresnel modelling approach for multilayered thin films. Each layer is assumed to have a homogenous refractive 
index throughout its thickness.

#### Experimental

1) Measure the angular reflectivity trace (including the TIR region) of the cleaned sensor (or use a previous representative one, depending on required accuracy).

2) Add the layer of interest or perform desired surface treatment and measure the sample again. 

#### Modelling (Text WIP)

1) Load the measurement file corresponding to the cleaned sensor. Select the desired wavelength and corresponding 
literature values for the refractive index and thickness of each material layer for the sensor (defaults are loaded). 

2) Choose which layer and property should be fitted to match a model background to the reference measurement. For empty 
metal sensors, choose the imaginary part of the refractive index (extinction coefficient, _k_) for the metal layer (the 
plasmonic metal layer thickness should not be altered). For all other types of layers, select the surface layer 
thickness (or refractive index) as variable to be fitted. Provide its refractive index as a constant (or vice versa if 
thickness is selected) and give an initial guess as well as a lower and upper bound of the variable thickness (or 
refractive index) for the fitting algorithm. 

3) Run the fitting once to see the plot. Use the interactivity in the plot to adapt the selected range of values to fit 
against (typically around the reflectivity minimum) until the fit looks good.

4) If needed, add a new surface layer to the model on top of the previous and repeat the procedure. 

### Exclusion height determination (Text WIP)

The non-interacting height probe method can be utilized in MP-SPR measurements to determine the exclusion height of a 
probing particle for a swollen layer in solution. The method is based on the following peer-reviewed papers:

Schoch, R. L. and Lim, R. Y. H. (2013). Non-Interacting Molecules as Innate Structural Probes in Surface Plasmon Resonance.
_Langmuir_, _29(12)_, _4068–4076_. 
https://doi.org/10.1021/la3049289

Emilsson, G., Schoch, R. L., Oertle, P., Xiong, K., Lim, R. Y. H., and Dahlin, A. B. (2017). 
Surface plasmon resonance methodology for monitoring polymerization kinetics and morphology changes of brushes—evaluated with poly(N-isopropylacrylamide). 
_Applied Surface Science_, _396_, _384–392_. 
https://doi.org/10.1016/j.apsusc.2016.10.165

The non-interacting probe method have the following requirements: 
   * The probing particle does not interact with the layer of interest (it should be excluded from entering into it and not stick to it)
   * The layer is sufficiently thinner than the sensor decay length (enough to give a significant response to the injected probe)

Below follows the principle of the technique.

#### Experimental (Text WIP)

1) Measure a reference sample (clean sensor) in air or the solvent of interest and acquire the optical parameters that 
will be used as a background for the analysis. For best accuracy, this should ideally have been treated as similarly 
to your sample as possible without containing the layer(s) of interest.

2) Take a dry scan of the sensor with the polymer layer. The height of 
the dry layer can later be modelled and used as a lower bound when determining the exclusion height. 

3) Measure the SPR/TIR response from a surface containing your layer of interest while making 2-3 repeated injections of 
the probing particles for 5-10 minutes each. High contrast between the response with and without the probe is required 
for accurate determination of the exclusion height. For protein or polymer probes, 10-20 g/L is typically used.
Verify that the baseline returns to the same level prior to probe injections (within ~ 10 millidegrees).

#### Modelling (Text WIP)

1) Select a previous fresnel modelled background to obtain fit parameters that match the measurement well in liquid. It is recommended to fresnel model a liquid scan (it will not be physically correct) immediately before the first probe injection and select this as background.

2) The program assumes the background was modelled with the layer refractive index in air will suggest the resulting layer thickness as a lower bound and six times this value as upper bound when later determining the exclusion height. Adjust bounds manually if necessary.

3) Select injection points around a range from the SPR sensorgram covering the probe injection of interest. Plot the SPR angle vs TIR angle traces 
during probe injections to verify that the probe is non-interacting (linear relationship with low degree of hysteresis 
means non-interaction).

4) Select suitable ranges before, during and after probe injections to generate the averaged angular reflectivity traces used for
fresnel model fitting. Check that the fits of your model and parameters for all pairs of solvent/probe are sufficient in
the region of interest (typically around the reflectivity minimum). If the fit is poor, try the following:
    * Try a different initial guess of the swollen layer refractive index (typically n = 1.33 - 1.5 in water)
    * Broaden or narrow the range of the region of interest
    * Tweak the offset in intensity
    * Select new ranges around the probe injections (as stable as possible).
    * Use a different reference sample (or optimize your reference fitting)

5) Select a reasonable height range for the layer of interest in its solvent (this you can iterate on if needed).
Run the calculations of the exclusion heights for each probe injection. The calculations are run twice for each injection,
once based on the response before the probe injection and once for after the probe injection has been rinsed. In each case,
the calculations will yield a modelled swollen layer height for every refractive index value within the provided range, 
with and without the presence of the probe. The exclusion height is then found where these graphically intersect.

If the two exclusion height values differ significantly between each other for each probe injection, the probe likely 
interacts with something on the sample over time, partly adsorbs to the surface, or needs longer time to rinse properly 
from the flow cell (shift solvent range to further after probe rinsing). 

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


