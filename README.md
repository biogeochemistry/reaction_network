### Chemical reaction network toolbox ###

Prototype of the reaction network toolbox. This toolbox will give a modeler the possibility to easy specify chemical reactions and parameters in the text file. The main purpose of the toolbox is to solve the coupled nonlinear partial differential equations which represent the transport and reactions processes in the porous aqueous media:
 
![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7B%5Cpartial%20%28%5Cvarepsilon%20C_i%29%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%20%5Cleft%28%5Cvarepsilon%20D_i%20%5Cfrac%7B%5Cpartial%20C_i%7D%7B%5Cpartial%20x%7D%20%2B%20%5Cvarepsilon%20D_%7Bbio%7D%20%5Cfrac%7B%5Cpartial%20C_i%7D%7B%5Cpartial%20x%7D%5Cright%29%20-%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%20%5Cleft%28%5Comega%5Cvarepsilon%20C_i%5Cright%29%20%2B%20%5Cvarepsilon%5Csum%20R%28x%2Ct%2CC_i%2C...%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)


#### The toolbox is still under the development and collaborators are really welcome!

Please, contact me, if you know how to improve or have ideas about some cool features. 
Generally, I expect collaborators to create the pull request with new feature/improvement! 

#### Main development line: ####

- [x] 1D ODE module for advection-diffusion;
- [x] Parallel comp: OpenMP for linear solver;
- [x] 1D ode 6th order approximation of derivatives.;
- [ ] Reaction parsing module: create parsing string mechanism for obtaining dC/Dt for all reactions. User can provide list of reactions with R values as a text file.;
- [ ] 2D ODE extension;
- [x] PDE 1D extension;
- [x] Coupled PDE;
- [ ] OpenCL test for small Ax=b problems;
- [ ] The output of the model should include the concentrations and mass balances of all simulated substances, the magnitudes of all imposed and simulated parameters, and the mass fluxes of all simulated processes.;
- [ ] Submodules;
- [ ] Someday in the future;
- [ ] Python wrappers;

#### Reaction parsing module: 0% ####
- [ ] model could be defined by means of an input file with a system definition (substances, processes and pertinent coefficients), a computational grid, an initial composition, timers, flow fields (a dynamic water balance), dispersion coefficients, loads and meteorological forcing. All input parameters in the model can be specified as constants or temporally and/or spatially varying parameters.;
- [ ] Boost. RegExp;

#### pH module 100% ####
- [x] Levenberg-Marquardt minimization algorithm;


#### PDE 1D extension 100% ####
- [x] PDE solver C++;
- [x] Test solution with Mathematica;
- [x] Write tests for PDE module;


#### Coupled System of PDEs: 100% #####
- [x] update concentrations according to R term;
- [x] Solve PDE for an each species..;
- [x] Parallel solution of each species with OpenMP(2 Species = 2 processors);


#### Submodules: 0% ####
- [ ] PHREEQC;
- [ ] FABM;


#### 2D ODE extension: 58% ####
- [x] 2d nodes North, south, West, East, Center;
- [x] add number to Node?;
- [x] 2d finite difference grid;
- [x] interior and BC nodes vs all nodes in 2d?;
- [x] Boundary conditions;
- [x] formation of 2d grid + tests;
- [x] assigning of the functions to 2D boundary + tests;
- [ ] alpha betta gamma coefficients in 2D (North,West,South,East,Center nodes);
- [x] Populate matrix + tests;
- [x] Populate vector +tests;
- [ ] apply solver + tests;
- [ ] test on Poisson's eq with known solution;


#### Someday in the future: 67% ####
- [ ] specify Boundary conditions for 2D or 3D as a function f(x,y);
- [ ] solve linear systems fully in parallel (Hybrid MP and MPI);
- [ ] auto-calibration (MCMC-DREAM vs MCMC) influenced by LibBi(http://libbi.org/index.html);


#### Python wrappers: 33%#####

- [ ] cython;
- [x] swig: node.hpp:10: Warning 325: Nested struct not currently supported (coordinates ignored);
- [ ] boost.python;


##### The MIT License (MIT)

Copyright (c) 2016 Igor Markelov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
