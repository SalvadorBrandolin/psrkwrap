Julia 1.8.3 and Python 3.10

Julia: 
Clapeyron v0.5.7
IJulia v1.23.3
PackageCompiler v2.1.5
PyCall v1.95.1
PyPlot v2.11.0

Python
julia==0.6.0
numpy==1.21.1
scipy==1.10.0

The [pyjulia](https://github.com/JuliaPy/pyjulia) library is used to make calls 
to Julia and evaluate the PRSK eos functions easily and quickly, without 
dealing directly with Julia and memorizing all the functions. Note that using 
the Python wrapper introduces an overhead that can make the execution up to 10 
or 100 times slower.

The idea of this mini-library is to perform the flashes for my PhD with ease 
and convenience.

As of December 27, 2022, Clapeyron does not have the following functions 
implemented:

Ideal gas heat capacity DIPPR equation 107
Mathias Coppeman's ${\alpha}$ function.
Defining DIPPR107 is complicated because it requires defining the ideal
Helmholtz free energy function. Thus, it is easier to take the DIPPR function 
and adjust the ReidIdeal cubic polynomial and use that as the cp function.

On the other hand, for Mathias Coppeman's function, there is no other option.

This function is implemented in the 
[mathiascopeman.jl](/psrk/julia/mathiascopeman.jl) file, imitating the other 
alpha functions in Clapeyron. It has been tested.


The installation is not very complex but has its quirks.

INSTALLATION

Easy way  
In your python interpreter:

```python
from psrk import install

install()
```

Manual way:
1) Install Julia. This is a bit cumbersome; you need to download the Julia 
folder and add it to the .bashrc path. For example:
    export PATH="$PATH:~/julia-1.8.3/bin"
or the desired path for Julia...

2) From the terminal, run:
```bash
    julia
```
Julia should run if it is installed correctly. Once there, execute:
```Julia
import Pkg; Pkg.add("Clapeyron")
import Pkg; Pkg.add("PyCall")
```  
This installs the Clapeyron library and PyCall, which is necessary for calling Julia from Python.

3) Create a Python virtual environment to work comfortably and execute from the terminal:
```bash
    pip install julia
```
Then, in Python:
```python
    import julia
    julia.install()
```
This installs PyJulia.

4) 
Now, calling Julia from Python won't work due to incompatibilities. 
The solution is to create a system image where Julia is compiled with 
Clapeyron and PyCall. See the documentation at: 
https://pyjulia.readthedocs.io/en/latest/sysimage.html

Execute from the library's root:
python3 -m julia.sysimage psrk/sys.so

This command takes some time to run. It creates a compiled version of Julia and 
saves it in a sys.so file (system image). 

I will create a sys.so in the psrk folder within this library so that it can 
be used directly. However, it may be necessary to create one on each machine 
because the compiled version may be different, even if the

 To check if it works and learn a bit more, refer to the file calling_julia.ipynb.

    Install the library in the virtual environment:
        pip install .

    Optional:
    you can run the tests to see if all works well.

    pip install -r requirements_dev.txt

    tox -r

    Ejecutar desde bash:  
    python3 -m julia.sysimage sys.so