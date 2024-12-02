[![DOI](https://zenodo.org/badge/402256220.svg)](https://doi.org/10.5281/zenodo.14262568)

# OpenQKDSecurity Version 2.0

OpenQKDSecurity is a software package based in MATLAB that allows users to calculate key rates for quantum key distribution (QKD) protocols using the [Winick et al. framework](https://quantum-journal.org/papers/q-2018-07-26-77/) ([arXiv](https://arxiv.org/abs/1710.05511)). It is extensible, allowing for user-defined protocols to be implemented, and modular, allowing for users to change specific aspects of a protocol. Our software can be used to interface with experimental data, demonstrate the theoretical scalability of protocols in various conditions, and optimize parameters to maximize key rate. Its modular structure helps break down the colossal task of calculating key rates into small areas that require only domain specific knowledge. Therefore, no single person must be an expect in all areas.

We break down a numerical QKD security proof into a preset and 5 differnt modules described bellow:

| Layer       | Difficulty            | Audience                                                                                                                                                             | Description                                                                                                                                                                                                                                                                                            |
|-------------|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Preset      | Basic to Intermediate | Everyone                                                                                                                                                             | The Preset sets up a QKD protocol by defining parameters, options and modules to use. Use cases range from simply picking a predefined preset and tweaking a parameters to selection of modules and technical options.                                                                                 |
| Description | Intermediate          | Users with intermediate knowledge of QKD protocols, including proof techniques.                                                                                      | Using the initial parameters given in the preset, a description module constructs complete description of a QKD protocol. Measurement observables, announcement structure, and the key map are typically defined on this layer.                                                                        |
| Channel     | Intermediate          | Users with intermediate to advanced knowledge of channel simulation or computing key rate from experimental data.                                                    | The channel module produces expectation values or observed frequencies from Alice and Bob's measurements on their exchanged states. These values can be obtained from numerical simulation or imported from an external source. Experimental data should also be loaded in this layer.                 |
| Key Rate    | Advanced              | Users with advanced knowledgeable of numerical proof techniques in QKD.                                                                                              | This layer is responsible for determining the key rate of a class of protocols based on the given parameters. Proof techniques such as decoy state analysis and squashing should be implemented in this layer. This layer is also responsible for processing input to be forwarded to the math solver. |
| Math Solver | Advanced              | Users with advanced knowledge of convex optimization.                                                                                                                | Calculates a strict lower bound on the quantum relative entropy between the key and Eve via convex optimization.                                                                                                                                                                                       |
| Optimizer   | Advanced              | Users with knowledge of general optimization techniques with the intent to minimize the number of calls to the objective function needed to find the global maximum. | Routines to perform unstructured optimization on selected parameters, using the key rate as the objective. This layer is given a wrapped version of the protocol with options modified to reduce precision but decrease calculation time.                                                              |

Once parameters, modules, and options are packaged into a preset file, running a numerical key rate calculation is as simple as,
```matlab
% Qubit BB84 prepare and measure. We use schmidt decomposition to reduce
% the dimension of Alice's state for enchanced speed and stability.
qkdInput = BasicBB84Alice2DPreset();

%run the QKDSolver with this input
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("BasicBB84Alice2DResults.mat","results","qkdInput");

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results)
```

Further documentation and examples may be found in the included user guide pdf.

## Installation

Before installing the software, ensure you have the latest version of MATLAB installed on your machine.
Our software requires *at least version 2020b* for full functionality, but installing the latest version is preferable.

Our software has the following dependencies for its default settings:

- [CVX](https://cvxr.com/cvx/download/) v2.2, a library for solving convex optimization problems in MATLAB.
- [QETLAB](https://github.com/nathanieljohnston/QETLAB) *above* v0.9, a MATLAB toolbox for operations involving quantum channels and entanglement. Note that you cannot use the version from their website as it has a bugs associated with implementing Choi matrices. *You must use their latest copy on Github*. At the time of writing, this means downloading their code with the green "code" button and *not* the v0.9 release.
- [ZGNQKD](https://www.math.uwaterloo.ca/~hwolkowi/henry/reports/ZGNQKDmainsolverUSEDforPUBLCNJuly31/) solver (optional), an alternative to our Frank-Wolfe solver for quantum relative entropy. Currently, it only supports equality constraints. Download the zip and the "Solver" folder and sub-folders to your path.
- [MOSEK](https://www.mosek.com/) (optional), a more advanced semidefinite programming (SDP) solver than the default (SDPT3) used in CVX. Note that the MOSEK solver can be downloaded together with CVX, but requires a separate license to use. See [this page](https://cvxr.com/cvx/doc/mosek.html) for more information.

Please refer to the documentation of each of these software packages for installation instructions.

To install the openQKDSecurity package, download the repository on GitHub as a zip file and unzip to a preferred directory. This can be done on the main page of the repository by pressing the green “Code” button at the top of the page and then selecting "Download ZIP".

Next, ensure that our software and the above dependencies are on your MATLAB path.To place a folder on your path, navigate to it in MATLAB, then right-click and select "Add to Path\textgreater Selected folder and Subfolders". Make sure you do this for OpenQKDSecurity, QETLAB and CVX. We also recommend you run "cvx_setup" to check if CVX is properly configured and which solvers are available.

Before you close MATLAB, go to "Home\>Environment\>Set Path" and click "save" at the bottom. If you don't, then you have to add everything to your path each time you restart MATLAB.

> [!IMPORTANT]
> We strongly encourage you to run "testInstall.m" at this point to check for basic installation issues.

## Current Status
Currently in this version of the software, we provide the following example protocols:
* Qubit BB84 with no loss, where Alice is 4-dimensional
* Qubit BB84 with no loss, where Alice has been reduced to two dimensions (this reduction is used in all of the following protocols)
* Qubit BB84 with loss
* Qubit BB84 with a finite number of signals sent
* Weak coherent pulse based BB84 with decoy-state analysis

We also provide two math solver modules:
* FW2StepSolver, which uses a Frank-Wolfe iteration algorithm followed by a linearization and conversion to dual in order to calculate a reliable lower bound on the quantum relative entropy. This is based on work done by our group published at https://arxiv.org/abs/1710.05511v2. 
* FRGNSolver (requires the ZGNQKD solver), which uses facial reduction and a Gauss-Newton method to simplify and solve the SDP to arrive at quantum relative entropy. This solver was written by members of Henry Wolkowicz's group at the University of Waterloo and is published at https://arxiv.org/abs/2104.03847v2. Note that this solver can only handle strict equality constraints at the moment.

In the future, we plan to add more protocols. We also welcome custom-defined protocols and solver modules.

## Contributing


If you would like to contribute, please contact Professor Lütkenhaus at lutkenhaus.office@uwaterloo.ca.

## License
OpenQKDSecurity is licensed under the MIT license. See the LICENSE file for details.

## Acknowledgements
This work has been performed at the Institute for Quantum Computing (IQC), which is supported by Innovation, Science and Economic Development (ISED) Canada. This research has been supported by NSERC Collaborative Research and Development (CRD) Program and Discovery Grants Program. Financial support for this project has been partially provided by Huawei Technologies Canada Co., Ltd.
