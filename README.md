# Minimal model for human atrial action potentials

This repository contains the code for a minimal model of human atrial cardiomyocytes. Using the [Hodgkin/Huxley formalism](http://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) for creating models of cellular electrical activity using a system of ordinary differential equations, this code provides a realistic simulation of activity in a typical human atrial cell. 

### Notes:

* The code as it appears here is not a direct implementation of the original minimal model (MM) presented by Bueno-Orvio et al. in their 2008 [paper](http://www.sciencedirect.com/science/article/pii/S0022519308001690). Some slight modifications to the original equations were made that allowed the model to be more closely fit to data from a patient.
* All of the different FK4V_5p0_... files are versions of the original code set up to be used with different versions of the genetic algorithm (not included here) that we used to fit the model to the patient data.
* Any code with "cable" in the name includes a series of single cell models coupled together with a partial differential reaction-diffusion equation.
* I didn't write main.cpp...
