# Changelog

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