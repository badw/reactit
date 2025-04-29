[![Tests](https://github.com/badw/reactit/actions/workflows/tests.yml/badge.svg)](https://github.com/badw/reactit/actions/workflows/tests.yml)




### `reactit` -  **React**ion **It**erator <img src="./static/reactit.png" width="300" align="right">

`reactit` is a code for generating all possible reactions between species in an efficient manner.


#### Installation 

```
git clone https://github.com/badw/reactit.git 
cd reactit 
pip install . 
```

#### Examples 

an example Jupyter Notebook can be found here: 


`./examples/reactit.ipynb`


simple usage: 

```
from reactit import ReactionGenerator

rg = ReactionGenerator( compounds = {0:'CO2',1:'H2O',2:'H2',3:'CO'} )
reactions = rg.iterate(max_length=4)

reactions
```

> ['1 CO2 + 1 H2 = 1 H2O + 1 CO']

support for [https://github.com/bjodah/chempy.git](chempy) `Equilibrium` and [https://github.com/materialsproject/pymatgen.git](pymatgen) `BalancedReaction` objects are included: 

```
chempy_reactions = rg.to_chempy()
chempy_reactions[0]
```
> Equilibrium(CO2 + H2 ↔ CO + H2O)

```
pymatgen_reactions = rg.to_pymatgen()
pymatgen_reactions[0]
```
> BalancedReaction(H2 + CO2 -> H2O + CO)



#### 
Todo: 
- [ ] ability to pass `SMILES` strings to form `SMILES` equations 
- [ ] multiprocessing for particularly large reactions 
- [ ] save to file 
- [ ] `CLI`
