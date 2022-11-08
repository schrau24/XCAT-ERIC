# MRXCAT with ERIC: extra-dimensional respiration with inflow of contrast

#### a MATLAB-based app designed for testing of MR sampling and reconstruction of abdominal dynamic contrast-enhanced MRI

<img src="/utilities/icons/PancreasSim.gif?raw=true" width="700px">

![Main screen](utilities/icons/mainScreen.png)

Required MATLAB toolboxes:
* *Image Processing Toolbox*
* *Parallel Computing Toolbox*

Additional toolboxes:
* [gpuSparse](https://github.com/marcsous/gpuSparse) and [nufft_3d](https://github.com/marcsous/nufft_3d) from Mark Bydder. Note these are provided with the package but may need to be recompiled within MATLAB
* [bart](https://github.com/mrirecon/bart), either on a linux system or using a [Windows-based install](https://bart-doc.readthedocs.io/en/latest/install.html). This code was tested using bart *version 0.4.03*

To start, open MATLAB at the folder containing *MRXCATwERIC.mlapp* and type:
```
MRXCATwERIC
```