# spatial_rcs
Code associated with "Inference of malaria reproduction numbers in three elimination settings by combining  temporal data and distance metrics" by Routledge et al. https://www.nature.com/articles/s41598-021-93238-0

## Running the tutorial code
The file tutorial.R contains the code required to set up the algorithm, including TensorFlow and Tensorflow-probability. This file sources the function_genRC.R code containing the algorithm used to generate estimates of individual reproduction numbers as well as the posteriors for the parameters delta, epsilon and alpha used in the model.
This file also loads Res.csv, which is a simulated line list, converted  in the tutorial into the appropriate matrices for inputting into the model.



