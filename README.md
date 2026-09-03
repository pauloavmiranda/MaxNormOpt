# Supplementary materials of the paper "A Comparative Analysis of Max-Norm Optimization Frameworks in Image Segmentation" published in 39th Conference on Graphics, Patterns and Images (SIBGRAPI):

## Abstract

Graph-based frameworks by combinatorial optimization can handle image segmentation as a graph partition problem subject to soft and hard constraints. In this work, we conduct a comparative analysis of three different algorithms proposed for the optimization of max-norm objective functions in image segmentation: Max-Norm Optimization by Removal of Atoms (MNORA), Oriented Image Foresting Transform (OIFT), and Oriented Relative Fuzzy Connectedness (ORFC). OIFT and ORFC are efficient methods specially designed for image segmentation on directed graphs by optimum oriented cuts. Following the ideas of MNORA, we also theoretically extend OIFT and ORFC to handle a more general class of problems satisfying a condition analogous to the submodularity condition encountered in max-flow/min-cut optimization. Quantitative experiments are performed on natural images at different resolutions on graphs with asymmetric weights, with source code available.

## Authors

- Justino F. Alvarez
- F&aacute;bio A.M. Cappabianco
- Paulo A.V. Miranda

## Source code

Our source code is available in the folder **"MaxNormOpt"**.
The code was implemented in C/C++ language, compiled with gcc 13.3.0, and tested on a Linux operating system (Ubuntu 24.04.4 LTS 64-bit), running on an Intel® Core™ i5-10210U CPU @ 1.60GHz × 8 machine. 
The code natively only supports images in the PPM format (_Portable Pixel Map_).
The [ImageMagick](https://imagemagick.org/) command-line tools are required to convert images from different formats.



## Contact

If you have any doubts, questions or suggestions to improve this code, please contact me at:
**pmiranda@ime.usp.br**
