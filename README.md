RGIsearch
=========

RGIsearch is a program designed to find invariants of renormalization group equations of any physical theory. An elaborate explanation on the purpose of these invariants and the methodology behind the program can be found in my [thesis](http://www.ru.nl/publish/pages/517373/thesis_rob_verheyen.pdf)

The results of the application of this code to the renormalization group equations to the MSSM will soon be published.

RGIsearch currently implements the beta-functions of the SM and the MSSM, but any set of beta-functions can be implemented. 

#Compilation
For first-time compilation, running 

``` 
make
```

will generate all binary files and link them to an executable. This can then be run with

``` 
./RGIsearch 
```

If any Beta-functions, settings or any of the sources are changed, ```make``` should update the binaries. If for some reason it does not, or it is otherwise required, one can do

```
make clean
make
```

#Program Settings
The program settings can be found in `settings.h`.

```const int MAX_TERM```

This parameter controls the amount of different parameters allowed in a monomial term. For instance, if it is set to 2, terms like xy^2 can occur, but xyz is forbidden. 
Increasing this parameter can have huge influence of the required computation times and memory.

```const int FILTER_THRESHOLD```

This parameter controls the filtering mechanism. While searching through the dimensionalities defined by the first two parameters, the algorithm might encounter invariants that are 
products of previously found invariants. A filtering algorithm is implemented to avoid these products. If this parameter is set to a higher value, the filtering algorithm will perform 
a more elaborate attempt to find all possible products of previously found invariants to compare against any newly found invariants. The default value should almost always be sufficient.

```bool INCLUDE_TWO_LOOP```

Set to true to search for invariants at two-loop level. Set to false to search for invariants at one-loop level

```bool FACTORIZE```

Set to true to include the algorithm that searches for factorized polynomial invariants. Refer to the above sources for details on this algorithm. 
The code does not support factorization for two-loop level searches.

```bool REPORT```

Set to true to make the program report its current activities.

#Beta-Functions
The beta-fcuntions can be found in `equations.cpp`.

To define the beta-functions of a theory, the physical parameters of the theory first need to be defined. This is done with:

```
Param newPar("parName", parSize);
```

ParSize is the size of the parameter (single number, 2x2 matrix, 3x3 matrix etc..).

If a parameter is complex, its daggered counterpart must be defined. This is done with:

```
Param newParDagger = dagger(newPar);
```

Several methods are included in the Param class to allow for simplifications such as those required for going from the MSSM to the pMSSM. They include:
1. `botRight()` Makes everything except for the bottom right component equal to zero.
2. `diag()` Makes all offdiagonal terms equal to zero.
3. `botRightDiag()` Makes all offdiagonal terms equal to zero, and the all components except for the bottom right degenerate.


After defining the physical parameters of the theory, one can define its beta-functions. These are stored as a vector of the objects BetaFunc (`vector<BetaFunc> nameOfBetaFuncs`). 
The beta-function of a parameter can then be defined as:

```
BetaFunc bNewPar(newPar);
bNewPar = <polynomial> 
```

The polynomial can be constructed using regular arithmetic involving the physical parameters of the theory. Irrational numbers (a/b) can be represented as `ir(a,b)`. 
A trace function is available as `Tr()`. For the complex conjugate parameters, their beta-functions must be included as:

```
BetaFunc bNewParDagger = Conjugate(bNewPar);
```

Finally, the algorithm can be initiated by calling:

```
findInvariants(1loopBetaFuncs, 2loopBetaFuncs)
```

If only one-loop beta-functions are available, it is possible to pass an empty vector as 2loopBetaFuncs as long as `bool INCLUDE_TWO_LOOP` is set to `false`. 




