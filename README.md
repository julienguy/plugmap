python scripts to evaluate fiber losses in eBOSS
```
usage: plugmap.py [-h] -p PLUGMAP --holes HOLES [-c CONFIG] [--oha OHA] -o
                  OUTFILE

optional arguments:
  -h, --help            show this help message and exit
  -p PLUGMAP, --plugmap PLUGMAP
                        path to plPlugMapP-XXXX.par (default: None)
  --holes HOLES         path to plateHoles-00XXXX.par (default: None)
  -c CONFIG, --config CONFIG
                        path to configuration file (default: None)
  --oha OHA             observed hour angle, overwrite config (default: None)
  -o OUTFILE, --outfile OUTFILE
                        path to output ASCII file (default: None)
```

Example :
```
./plugmap.py -p $PLATELIST_DIR/plates/0072XX/007285/plPlugMapP-7285.par --holes $PLATELIST_DIR/plates/0072XX/007285/plateHoles-007285.par -o myPlugMapP-7285.list --oha -12.966267 --config eboss-r.conf
```
with

eboss-r.conf :
```
OBSERVED_GUIDE_WAVE 5000. # poorly known
OBSERVED_QSO_WAVE 7450. # to study r1/r2
OBSERVED_STAR_WAVE 7467. # to study r1/r2
OBSERVED_LRG_WAVE 7498. # to study r1/r2
OBSERVED_ELG_WAVE 7498. # to study r1/r2
FIBER_LOSS_SIGMA_ARCSEC 0.87 # valid for a seeing of 1.5 arcsec FWHM
```


