# Supplementary materials of the paper "A Comparative Analysis of Max-Norm Optimization Frameworks in Image Segmentation" published in 39th Conference on Graphics, Patterns and Images (SIBGRAPI):

## Abstract

Graph-based frameworks by combinatorial optimization can handle image segmentation as a graph partition problem subject to soft and hard constraints. In this work, we conduct a comparative analysis of three different algorithms proposed for the optimization of max-norm objective functions in image segmentation: Max-Norm Optimization by Removal of Atoms (MNORA), Oriented Image Foresting Transform (OIFT), and Oriented Relative Fuzzy Connectedness (ORFC). OIFT and ORFC are efficient methods specially designed for image segmentation on directed graphs by optimum oriented cuts. Following the ideas of MNORA, we also theoretically extend OIFT and ORFC to handle a more general class of problems satisfying a condition analogous to the submodularity condition encountered in max-flow/min-cut optimization. Quantitative experiments are performed on natural images at different resolutions on graphs with asymmetric weights, with source code available.

## Authors

- Justino F. Alvarez
- F&aacute;bio A.M. Cappabianco
- Paulo A.V. Miranda

## Source code

Our source code is available in the repository folder **"MaxNormOpt"**.
The code was implemented in C/C++ language, compiled with gcc 13.3.0, and tested on a Linux operating system (Ubuntu 24.04.4 LTS 64-bit), running on an Intel® Core™ i5-10210U CPU @ 1.60GHz × 8 machine. 
The code natively only supports images in the PPM format (_Portable Pixel Map_).
The [ImageMagick](https://imagemagick.org/) command-line tools are required to convert images from different formats.

The public image dataset with ground truth used in the comparative evaluation is available in subfolder **"MaxNormOpt/dat"**.
It consists of 240 natural images of household objects taken with a mobile phone, divided into 3 different resolutions (300&times;400, 600&times;800, and 1200&times;1600 pixels), aiming to measure the execution times and memory consumption of the different max-norm optimization algorithms evaluated.
All commands executed in the experiments to generate the results are available in the subfolder **"MaxNormOpt/exp"**.

To compile the program, enter the repository folder **"MaxNormOpt"** and type **"make"**.
To segment an image, you must run the **"OptMaxNorm"** executable, which will be located inside the repository folder after compilation.

### usage:

```
OptMaxNorm <image> <file_seeds> [polarity] [gtruth]
image ............. the color input image (.ppm, .png, .jpg),
file_seeds ........ the file with seeds (.txt),
polarity .......... boundary polarity in the interval [-1, 1] (default = 0.5),
gtruth ............ ground truth image in the PGM format.
```

As output, the program generates the label images of the resulting segmentations (files **"label_mnora.pgm"**, **"label_oift-1.pgm"**, **"label_oift-2.pgm"**, and **"label_orfc.pgm"**).


### Program execution examples:

The following command computes the segmentations for the first low-resolution image.

```
./OptMaxNorm dat/low/fig01.png dat/low/seeds01.txt 1 dat/low/fig01_gt.pgm
```

The next command calculates the results for the third medium-resolution image.

```
./OptMaxNorm dat/med/fig03.png dat/med/seeds03.txt 1 dat/med/fig03_gt.pgm
```

All the above commands were taken from the files available in the subfolder **"MaxNormOpt/exp"**.
To run the complete experiment for a given dataset resolution (low, medium, or high resolution), simply run the corresponding script **"./exp/run_exp_low.txt"**, **"./exp/run_exp_med.txt"**, or **"./exp/run_exp_high.txt"** (don't forget to grant execution permission to these files).
The experimental output report, as a spreadsheet, will be produced in CSV format in file **"MaxNormOpt/report.csv"**, to be located in the repository folder.


## Other versions

To calculate the methods on a small custom graph, the executable **"OptMaxNorm_graph"** can be used.
For example, consider the following command:

```
./OptMaxNorm_graph dat/graph/graph02.txt
```

A version for use on Region Adjacency Graphs (RAGs) of superpixels is also made available via the **"OptMaxNorm_rag"** executable.


## Contact

If you have any doubts, questions or suggestions to improve this code, please contact me at:
**pmiranda@ime.usp.br**
