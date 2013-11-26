This is a prototype of WCS for JWST instruments.

Data is in FITS files but reference files are in JSON format.
The WCS information is in the FITS headers and in the reference files as follows:

Images:

- basic WCS in in primary header
- distortion is in distortion_image.json

Long slits:

- basic WSC is in primary header
- distortion is in distortion_image.json
- spectral models are in spec_wcs.json

IFUs:

- basic WCS of primary slit is in primary header
- basic WCS for all other slits is relative to the primary slit and is in wcsc_regions.json
- regions definitions on the detector in terms of polygons are in regions_miri.json
- distortions for each region are in distortion_regions.json
- spectral models for each region are in spec_wcs.json

