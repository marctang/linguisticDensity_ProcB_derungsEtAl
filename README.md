# Linguistic Density and Environmental Variables

This repository contains the R code used for computing all major results shown and discussed in the paper:

"Environmental factors drive language density more in food-producing than in hunter-gatherer populations" by Derungs et al. (currently in review)

## Prerequisites
In order to run these scripts locally the following two folder have to be downloaded:

- The R scripts require input data. All required input data can be downloaded from [this link]( https://www.dropbox.com/s/77i75beotws60bj/input.rar?dl=0). The downloaded input.rar has to be unpacked and stored in the same directory as the R code (there should be a folder called "input/").

- The scripts store all results in an output folder. An the folder, including all results from a previous execution of the code, can be downloaded from [this link]( https://www.dropbox.com/s/p0326myjj002dlm/output.rar?dl=0). The downloaded output.rar has to be unpacked and stored in the same directory as the R code and the input folder.


## Content
This repo contains six R scripts with the following purpose:
- 001_regularGrid.R: Computes the regular grid for three spatial resolutions
- 002_langPrep.R: Prepares the language data such that it can be used in the follow up analysis
- 003_langAggreg.R: Generates language counts for the grid points from the input language information
- 004_random_NS.R: Simulates 500 random distributions of points, which result in the expected counts per regular grid point
- 005_addingPredictors.R: Retrieves environmental information for all grid points
- 006_prediction.R: Predicts LD peaks and troughs with environmental information, using Random Forest

## Authors

* [**Curdin Derungs**](https://github.com/curdon) - *Implementation, Analysis*
* [**Balthasar Bickel**](http://www.comparativelinguistics.uzh.ch/en/bickel.html) - *Research Design*

