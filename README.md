# Filter - a C++ framework for solving non-linear inverse problems  

## Overview  

This framework is devoted to solution of problems which may be stated as follows.  
Given a (non-linear) model of the form  
```math
y = F(x, \theta),
```
(noised) measurements  
```math
\hat{y} \approx F(x_0, \theta)
```
and parameters 
```math 
\theta
``` find a good estimate of
```math x_0.
```
This formulation is applicable in various domains. Examples folder demonstrates application to estimating parameters of polynomial model (a toy example), and to solving geophysical problem of inverting airborne electromagnetics data. The latter one is detailed in several papers listed below.  

## Structure  
TODO  

## License  
This project is licensed under the Apache License 2.0 — see the LICENSE file for details.  

## External links  
[My first paper on airborne electromagnetics](https://www.geotechnologies.ru/publications/2023_AEM_kar_khliu.pdf)  
[Another paper, demonstrating several beautiful inversions](https://www.geotechnologies.ru/publications/2025_11_Induced%20polarization%20effects%20in%20frequency%20domain%20AEM%20data.pdf)  
[Also a geophysical paper of mine](https://www.geotechnologies.ru/publications/2025_9_The%20impact%20of%20AEM%20system%20configuration%20on%20inversion%20results.pdf)  
