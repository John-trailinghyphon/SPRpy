# Modelling of MP-SPR measurements with python

This program can be used to perform fresnel modelling of multi-parameter surface plasmon resonance (MP-SPR) measurements
performed using [Bionavis SPR instruments](https://www.bionavis.com/en/technology/why-choose-mpspr/). 

It is based on MATLAB implementations of the [transfer-matrix-method](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics)) 
by [Andreas Dahlin](https://www.adahlin.com/matlab-programs). 

## Performing fresnel modelling of non-swollen layers

The basic fresnel modelling approach for multilayered thin films. Each layer is assumed to have a homogenous refractive 
index throughout its thickness.

### Experimental

1) Measure the angular reflectivity trace (including the TIR region) of the cleaned sensor (or a representative one).

2) Attach the layer of interest and measure the sample again.

### Modelling

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

## The non-interacting probe method for height determination

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

### Experimental

1) Measure a reference sample (clean sensor) in air or the solvent of interest and acquire the optical parameters that 
will be used as a background for the analysis. For best accuracy, this should ideally have been treated as similarly 
to your sample as possible without containing the layer(s) of interest.

2) Take a dry scan of the sensor with the polymer layer. The height of 
the dry layer can later be modelled and used as a lower bound when determining the exclusion height. 

3) Measure the SPR/TIR response from a surface containing your layer of interest while making 2-3 repeated injections of 
the probing particles for 5-10 minutes each. High contrast between the response with and without the probe is required 
for accurate determination of the exclusion height. For protein or polymer probes, 10-20 g/L is typically used.
Verify that the baseline returns to the same level prior to probe injections (within ~ 10 millidegrees).

### Modelling

1) Use fresnel modelling to obtain model parameters that match the reference measurement.

2) The height of the layer in air can be modelled and used as a lower bound when later determining the exclusion height. 
You can skip this part if you instead want to set a lower bound manually.

3) Select a range from the SPR sensorgram covering the probe injection of interest. Plot the SPR angle vs TIR angle traces 
during probe injections to verify that the probe is non-interacting (linear relationship with low degree of hysteresis 
indicates non-interaction).

4) Select suitable ranges before, during and after probe injections to generate the averaged angular reflectivity traces used for
fresnel model fitting. Check that the fits of your model and parameters for all pairs of solvent/probe are sufficient in
the region of interest (typically around the reflectivity minimum). If the fit is poor, try the following:
    * Try a different initial guess of the swollen layer refractive index (typically n = 1.33 - 1.5 in water)
    * Broaden or narrow the range of the region of interest
    * Tweak the offset in intensity
    * Select new ranges around the probe injections (as stable as possible).
    * Use a different reference sample (or optimize your reference fitting)

5) Select a reasonable refractive index range for the polymer layer in its solvent (this you can iterate on if needed).
Run the calculations of the exclusion heights for each probe injection. The calculations are run twice for each injection,
once based on the response before the probe injection and once for after the probe injection has been rinsed. In each case,
the calculations will yield a modelled swollen layer height for every refractive index value within the provided range, 
with and without the presence of the probe. The exclusion height is then found where these intersect.

If the two exclusion height values differ significantly between each other for each probe injection, the probe likely 
interacts with something on the sample over time, partly adsorbs to the surface, or needs longer time to rinse properly 
from the flow cell (shift solvent range to further after probe rinsing). 

## (TODO) The dual-wavelength method
The dual-wavelength method can be used to determine the extension of swollen layers with unknown refractive index, 
based on the following requirements: 
   * The refractive index increment for the layer material for each wavelength is known.
   * The thickness of the measured layer is _much smaller_ than the decay length for each wavelength.
   * The sensitivity factors for each wavelength of the instrument is known (easily determined).

It is based on the peer-reviewed paper by:

Rupert, D. L. M., et al. (2016). Dual-Wavelength Surface Plasmon Resonance for Determining the Size and Concentration of Sub-Populations of Extracellular Vesicles. 
_Analytical Chemistry_, _88(20)_, 9980–9988. 
https://pubs.acs.org/doi/full/10.1021/acs.analchem.6b01860

TODO: Update this section when implemented.
