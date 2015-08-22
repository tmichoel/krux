# Download and installation instructions #

## Subversion ##

The recommended way to install kruX and keep up to date with any future versions is via the subversion system. If you are unfamiliar with this, see [here](http://subversion.apache.org/packages.html) for a list of relevant software clients. Detailed instructions for how to checkout kruX can be found under the [Source](https://code.google.com/p/krux/source/checkout) tab.

If you are unable to install kruX via subversion, you can [browse](https://code.google.com/p/krux/source/browse/trunk) the source code repository and download the files for your version of interest one-by-one (click on a file name, then click "view raw file", and save).

After checking out kruX via subversion, there will be a new directory called "kruX", with subdirectories

```
kruX/html/
kruX/matlab/
kruX/r/
kruX/python/
kruX/test/
```

The "html" directory contains documentation files, the "test" directory contains test data and scripts and the three others contain their respective kruX implementations.

## Matlab ##

To use kruX, it needs to be in the Matlab search path. This can be accomplished by executing:

```
>> addpath('<installation dir>/kruX/matlab/');
```

on the Matlab prompt, where <installation dir> is the location where kruX was checked out. See [here](http://www.mathworks.co.uk/help/matlab/ref/addpath.html) for more help on setting the search path. You can now change into the "kruX/test" directory and go [step by step](http://www.mathworks.co.uk/help/matlab/matlab_prog/run-sections-of-programs.html) through the test script `kruX_tutorial_matlab.m`, the output of which is described on this [tutorial page](https://krux.googlecode.com/svn/trunk/html/kruX_tutorial_matlab.html).


## R ##

Use the [source](http://stat.ethz.ch/R-manual/R-devel/library/base/html/source.html) command to load kruX, like

```
> source('<installation dir>/kruX/r/kruX.R')
```

The example script "test\_kruX.R" in the "kruX/test" directory illustrates how to run kruX.

## Python ##

Append kruX to your python path, like

```
>>> import sys
>>> sys.path.append('<installation dir>/kruX/python/')
>>> import kruX
```

The example script "test\_kruX.py" in the "kruX/test" directory illustrates how to run kruX.