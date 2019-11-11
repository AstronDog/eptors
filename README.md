# eptors.py

EPTA Pulsar Timing Outlier Rejection Scheme(eptors.py), is a python Script for eliminating incorrect and biased TOAs.

## Requirements ##

* Python 2.7
* [numpy](http://numpy.scipy.org)
* [scipy](http://numpy.scipy.org)
* [matplotlib](http://matplotlib.org), for plotting only
* [pandas](https://pandas.pydata.org)
* [statsmodels](https://www.statsmodels.org)
* [tempo2](http://tempo2.sourceforge.net)
* [sklearn](https://scikit-learn.org)

## Usage ##

Before running this code, please creating the tim file with command like:

```
pat -TP(Or FT) -C snr -f  "tempo2 IPTA" -A FDM -s TEMPLATE OBSERVATION > OBSERVATION.tim
```
and flag the parameters your want to fit within tempo2 in your par file.


To use this code, simply type

```
./eptors.py -e yourepe.par -t yourtoa.tim
```
and the code will work automatically to remove the outliers. More details and  annotations can be found in the code file.

## Contact ##

* [_Jun Wang_](mailto:jun.wang.ucas@gmail.com)