# pyHIIexplorerV2

 pyHIIexplorer is a package for extracting the integrated spectra of HII regions from integral field spectroscopy (IFS) datacubes. This version is based on [HIIexplorer](http://www.astroscu.unam.mx/~sfsanchez/HII_explorer/index.html) perl edition.
 
 ## Quick Description
 The detection of HII regions performed by pyHIIexplorer es bases on two assumptions: 1) \HII\ regions
 have strong emission lines that are clearly above the continuum emission and the average ionized gas emission across each galaxy. 2) the typical size of \HII\ regions is about a few hundreds of parsecs, which corresponds to a usual projected size of a few arcsec at the distance of our galaxies. These assumptions will define clumpy structures with a high Ha emission line contrast in comparison to the continuum.
 
### Input Parameters
1) Emission map (like Ha emission line map)
2) a minimum flux intensity threshold for the peak emission
3) a minimum relative fraction with respect to the peak emission to keep aggregating nearby pixels
4) a maximum distance to the ionized region peak emission to stop the aggregation
5) a minimum absolute intensity threshold in the adjacent pixels to continue the aggregation

For further information, check [Espinosa-Ponce, C. et al 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.1622E) and [Sánchez, S. F. et al 2012](https://ui.adsabs.harvard.edu/abs/2012A%26A...546A...2S)

## Requeriments 
* python 3.x
* os package
* csv package
* numpy package
* pandas package
* astropy package
* argparse package

If you do not get any error or warning message you have suscessfully installed pyHII_explorer. There are two example python scripts to run the segregation and extraction subscripts [run_pyHIIexp.py](run_pyHIIexp.py) and [run_pyHIIexp_mp.py](run_pyHIIexp_mp.py). If you find any problem, please send an email to C. Espinosa Ponce (cespinosa@astro.unam.mx).
