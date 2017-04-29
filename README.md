# EvoDormReview


## Re-running the analyses and re-generating the figures

Files in the data folder ending with ".zip" will have to be unzipped before the code can be run.


The code accepts the following arguments from the user.

**Flags**

### Order of operations

If you want to regenerate the data (warning, this is very computationally intensive and will take awhile to complete) start with step 1, otherwise run step 2 for the figures.

1) For example, if you want to generate the data, run the following commands (this will take a very long time).

  python runSimulations.py -d -3 -N 1000 -u  0.0001 -G 100 -g 10000 -r 100

  python runSimulations.py -d -5 -N 100 -M 100 -u 0.001 -R 0.01 -G 50 -g 1000 -r 100

	python runSimulations.py -d -6 -N 100 -M 100 -u 0.001 -G 100 -g 10000 -r 100

2) If you want to generate figure 4, run the following commands. Figure 1 is conceptual and contains no data.

  python makeFigs.py -f -2

  python makeFigs.py -f -3

  python makeFigs.py -f -4

  python makeFigs.py -f -5

  python makeFigs.py -f -6

## Dependencies

Python version 2.7.13 is used.

The following Python modules/versions are used in this analysis.

+ numpy 1.10.4

+ matplotlib 1.5.1

+ pandas 0.18.1

+ scipy 0.17.1

## The MIT License (MIT)

Copyright (c) 2017  William R. Shoemaker and Jay T. Lennon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
