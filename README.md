# Createive Component, STAT 6990
**Project:** A Spatial Model of Home School Enrollment Trends with Missing and Censored Data

**Terms** Spring 2024 & Fall 2024

**Code Contact:** Nick Grener, grener@gmail.com      

## Introduction  
 A recent article in The Washington Post describes a widespread increase in home schooling over the last decade in America, and furthermore indicates a 54% increase of home school enrollment in the state of Ohio (my home state) over the time frame from 2017 through 2022. The Washington Post’s analysis was unable to identify any
 important drivers of this increase, instead concluding that, “Home schooling’s surging popularity crosses every measurable line of politics, geography and demographics.” For my creative component
 project (STAT 6990), I wanted to investigate this phenomenon at the level of local school districts. This project provided me with an opportunity to learn about modeling spatial areal data, missing
 data mechanisms, Bayesian inference (including MCMC techniques), and spatio-temporal residual autocorrelation, via the analysis of an interesting real-world data set. For my analysis, I developed and evaluated a succession of linear spatial models in an effort to identify which predictors were associated with the recent trends in home school enrollment. In contrast to the newspaper’s findings, I discovered that- at least in the case of Ohio- some variables were significantly
 associated with the increase in homeschooling at the district level.

### Required R Packages
maps
tidycensus
units
tidyverse
sf
tigris
mapview
tmap
leafsync
ggspatial
remotes
readxl
ggplot2
car
MASS
Matrix
fdrtool
mvtnorm
truncnorm
spdep
here
epinet

## Parameters
The code needs to be run in the following order:
1. DataPrep.R
2. Analysis.R

## Acknowledgments
Nathan Wikle (University of Iowa) was an extremely helpful mentor throughout the course of this project. 