# NonConvexPLF_code


This repository is part of a Research Project developed by Salvador Guerrero García and directed by Javier Garcia Gonzalez (PhD from I.C.A.I). 

## Research Goal
The goal is to apply the linearization methods proposed in Nonconvex piecewise linear functions: Advanced formulations and simple modeling tools, (by Joey Huchette and Juan Pablo Vielma) in modeling the non linear relationships presented in Power Storage devices such as lithium bateries.


The paper is implemented by the papers authors in the Julia Lang function: PiecewiseLinearOpt
We implemented the paper developing our own library: myPiecewiseLinearOpt 

## understanding_library.jl

In this file the on-the-shelf library PiecewiseLinearOpt was compare with our implementation myPiecewiseLinearOpt. Different non lnear functions were minimized and maximized.

## hydro.jl

Small Economic dispatch model with a thermal power plant and a hydropower plant. The non linear relationship between the water level of the reservoir and the generated power was linealied using the paper implementation.

## Report.pdf

The results of the case study are shwon.
