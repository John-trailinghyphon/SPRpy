# Changelog

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