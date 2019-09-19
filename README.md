## Installation
For example, if you want to install it to `~/MyRepo`, do the following:
```
$ cd ~/MyRepo
$ git clone git@github.com:freemao/schnablelab.git
```
Then add `export PYTHONPATH=~/MyRepo:$PYTHONPATH` to your `~/.bashrc` file and run `$ source ~/.bashrc` to update your `PYTHONPATH`.

### Dependencies
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [Scipy](http://www.scipy.org)
* [scikit-learn](https://scikit-learn.org/stable/)

#### dependencies for some functions
* [Pillow](https://pillow.readthedocs.io/en/stable/)
* [opencv](https://pypi.org/project/opencv-python/)

All the dependencies can be installed through `easy_install`, `pip` or `conda install`.
