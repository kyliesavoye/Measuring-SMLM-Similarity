# Measuring-SMLM-Similarity

We propose a novel approach for comparing SMLM-derived marked point-clouds, integrating three distinct similarity measures:

1. Analysis of molecular coordinates
2. Examination of mark values
3. Assessment of mark distribution

This method represents the similarity between two marked point-clouds as a coordinate in 3D space. The distance of this coordinate from the origin correlates with the overall similarity of the clouds. We validated our approach using simulated data, demonstrating its ability to consistently identify dissimilarities across various point-clouds. Additionally, we applied the method to experimental data (available in folder "Experimental Data", specifically GP-marked point patterns from SMLM imaging using the polarity-sensitive membrane probe di-4-ANEPPDHQ. This framework is versatile, extending beyond SMLM to any modality generating marked point-patterns, and is applicable to other environmentally-sensitive probes.

In this repository, you will find the code to generate all the figures on our paper "Measuring the Similarity of SMLM-derived marked point-clouds". All the code can be found in the Jupyter Notebook, where sections 1 to 4 can be ran independently of each other. The folders "Method Schematic Figures", "Simulated Data Examples Figure", " Testing Method Figures" and "Experimental Data Figure" contain the figures resulting form the Jupyter Notebook. 
