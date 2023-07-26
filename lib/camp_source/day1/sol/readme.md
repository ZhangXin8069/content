# Solution for Exercise 1

### `c2pt_read.py`
* Read iog files inside `Data`
* Select pion two-point data, only real part
* Save to `c2pt.npy` file (numpy array file)


### `c2pt_meff.py`
* Load `c2pt.npy` file
* Calculate the effective mass
  * log
  * cosh


### `lsqfit_2pt.py`
* Load `c2pt.npy` file
* Jack-knife resample
* Averge forward / backward propagation
* `lsqfit`
* Plot the data and best fit for comparision
