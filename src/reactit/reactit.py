import numpy as np 
import itertools as it 
import math 
from sympy import symbols, Eq, solve    
import functools 
from fractions import Fraction 
import re 
from collections import defaultdict 
import warnings
import tqdm_pathos


class ReactionGenerator:
    
    def __init__(
            self,
            compounds
            ):
        self.nc = len(compounds)
        if isinstance(compounds,list):
            self.compounds = {i:c for i,c in enumerate(compounds)}
        else:
            self.compounds = compounds          

    @staticmethod
    def convert_number_strings(reactions):
        """
        converts a list of reaction strings from self.numeric_reaction_filter back into a list of numbers
        """
        for reaction in reactions:
            r,p = reaction.strip().split('],')
            r = tuple(int(x) for x in r.split('[[')[1].split(',') if x)
            p = tuple(int(x) for x in p.split('[')[1].split(']]')[0].split(',') if x)
            yield ((r,p))

    @staticmethod
    def numeric_reaction_filter(reaction:list)->str:
        """
        filters numeric reactions (e.g. [[1,2,3],[4]]) based on:
        1. if reactants = products 
        2. if any of the reactants are in the products (this is to avoid symmetry reactants 

        returns a string of that reaction which is later converted back to a list by "self.convert_number_strings"
        """
        reactants = sorted(reaction[0])
        products = sorted(reaction[1])
        if not reactants == products:
            if not any([x in reactants for x in products]):
                return (str(
                    sorted(
                        tuple(
                            (reactants, products)
                        )
                    )
                )
                )
            
    @staticmethod
    def string_reaction_filter(r:list)->list:
        """
        takes a reaction list string i.e. [['N2','O2','H2'],['HNO3']] and screens it for discrepancies between the elemental compositions of reactants and products. 

        impossible reactions are not returned.
        """
        re,pr = r
        c = []
        for i in re:
            rs = [x for x in i]
            for j in rs:
                try:
                    int(j)
                except Exception:
                    c.append(j)
    
        dr = set(dict.fromkeys(c))
    
        c = []
        for i in pr:
            ps = [x for x in i]
            for j in ps:
                try:
                    int(j)
                except Exception:
                    c.append(j)
    
        dp = set(dict.fromkeys(c))
    
        if dr == dp:
            return(r)

    def enumerate_combinations(self,max_length=4):#->tuple:
        """
        enumerates possible combinations of self.compounds as a numeric entity given a max length (minimum reaction length of 3)
        returns a list of lists
        """
        reactions = []
        for reaction_length in range(3,max_length+1):
            sizing = [
                x for x in it.combinations_with_replacement(
                    np.arange(
                        1, reaction_length
                    ),
                    2
                )
                if np.sum(x) == reaction_length
            ]    


            for i,size in enumerate(sizing):
                reactants = it.combinations(
                    [
                        x for x in range(self.nc)
                    ],
                    size[0]
                )    
    

                products = it.combinations(
                    [
                        x for x in range(self.nc)
                    ],
                    size[1]
                )
                combined = tuple(
                    it.product(reactants, products)
                    )
                filtered = []
                for v in map(lambda x: self.numeric_reaction_filter(x), combined):
                    if v:
                        filtered.append(v)
                filtered = list(set(filtered))
                reactions.extend(filtered)

        reactions = self.convert_number_strings(reactions)

        return(reactions)

    def convert_to_string(self,numeric_reactions:tuple)->list:
        """
        takes a numeric reaction i.e. [[0,1,2],[3]] and converts the numbers into a string given by self.compounds

        i.e.
        self.compounds = {0:'N2',1:'O2',2:'H2',3:'HNO3'} 

        [[0,1,2],[3]] -> [['N2','O2','H2'],['HNO3']]

        returns the converted list
        """
        converted = []
        for r in numeric_reactions:
            d = [[self.compounds[i] for i in r[0]],[self.compounds[i] for i in r[1]]]
            ds = self.string_reaction_filter(d)
            if ds:
                converted.append(ds)
        return(converted)
    
    @staticmethod
    def parse_molecule(formula:str)->dict: 
        """
        parses a molecule string i.e. 'H2O' into a dictionary broken down into elemental counts 
        i.e. {'H':2,'O':1}
        """
        # Regular expression to match elements and their counts 
        pattern = r'([A-Z][a-z]?)(\d*)' 
        matches = re.findall(pattern, formula) 
         
        atom_count = defaultdict(int) 
     
        for element, count in matches: 
            if count == '': 
                count = 1  # Default count is 1 if not specified 
            else: 
                count = int(count)  # Convert count to integer 
             
            atom_count[element] += count 
     
        return dict(atom_count)  
       
    @staticmethod
    def find_coefficients(fractional_coefficients:list)->int:
        """
        finds the common denominator and an integer factor to multiply reactions by
        returns an integer multiplyer to get new coefficients 
        """
        denoms = [Fraction(x).denominator for x in fractional_coefficients]
        return (
            functools.reduce(
                lambda a, b: a*b//math.gcd(a, b), denoms
            )
        )    

    def balance_reaction(self,reactants:dict, products:dict)->str:
        """
        balances a reaction given a dictionary of self.parse_molecule reactants and products. 
        by default it allows for undetermined systems and applies coefficients (if you desire something more agnostic go check out chempy.balance_reaction)
        returns a reaction string 
        """
        # Create a list of all elements
        reactant_keys = list(reactants.keys())

        product_keys = list(products.keys())

        elements = set()
        for i in list(reactants.values())+list(products.values()):
            for k in i.keys():
                elements.add(k)
        elements = sorted(elements)  

        # Set up symbolic coefficients
        reactant_coeff = [
                symbols(f'r{i}') for i in range(len(reactants))
        ]
        product_coeff = [
                symbols(f'p{i}') for i in range(len(products))
        ]
        equations = []    

        # For each element, write a conservation equation
        for element in elements:
            reactant_sum = sum(reactants[species].get(element, 0) * reactant_coeff[i]
                               for i, species in enumerate(reactants))
            product_sum = sum(products[species].get(element, 0) * product_coeff[i]
                              for i, species in enumerate(products))
            equations.append(Eq(reactant_sum, product_sum))    

        # Add an extra equation to fix one variable to 1 (to avoid trivial solution)
        equations.append(Eq(list(reactant_coeff)[0], 1))    

        # Solve the system
        solution = solve(equations, reactant_coeff + product_coeff, dict=True)    

        if not solution: ## this can be improved...
            return None  # No solution
        else:
            try:
                common_factor = self.find_coefficients(solution[0].values())
                new_coeffs = [x*common_factor for x in solution[0].values()]
                for i,(k,v) in enumerate(solution[0].items()):
                    solution[0][k] = new_coeffs[i]
            except TypeError:
                return None 
            reactant_string = ''
            product_string = ''
            rs = 0 ; ps = 0
            for symbol,coeff in solution[0].items():
                try:
                    if float(coeff) <= 0:
                        return(None)
                except TypeError:
                    return(None)
                
                rp,num = list(str(symbol))
                if rp == 'r':
                    if rs == 0:
                        reactant_string += f'{coeff} {reactant_keys[int(num)]}'
                    else:
                        reactant_string += f' + {coeff} {reactant_keys[int(num)]}'
                    rs += 1
                else:
                    if ps == 0:
                        product_string += f'{coeff} {product_keys[int(num)]}'
                    else:
                        product_string += f' + {coeff} {product_keys[int(num)]}'
                    ps += 1
            equation_string = reactant_string + ' = ' + product_string

        return(equation_string) 
    
    def balance_function(self,iterable:list):
            """
            balance_function checks the balance of a given list of strings, this is here for future multiprocessing
            """
            reactants = {
                molecule:self.parse_molecule(molecule) for molecule in iterable[0]
                }
            products = {
                molecule:self.parse_molecule(molecule) for molecule in iterable[1]
                }
            balanced = self.balance_reaction(reactants,products)
            if balanced:
                return(balanced)
    
    def iterate(
            self,
            max_length:int = 4,
            n_cpus:int = 4,
            **tqdm_pathos_kws:dict,
            ) -> list:
        """
        iterates possible reactions up to a maximum length (DEFAULT: 4) i.e. CO + H2O = CO2 + H2
        returns a list of reactions that are also accessible via self.reactions
        default 
        """

        numeric_reactions = self.enumerate_combinations(max_length = int(max_length)) #this is a generator
        strings = self.convert_to_string(numeric_reactions) # this is a list...
        #screened = tqdm_pathos.map(self.mp_function,strings)

        screened = tqdm_pathos.map(self.balance_function,strings,n_cpus=n_cpus,**tqdm_pathos_kws)
        balanced = []
        for reaction in screened:
            if reaction:
                balanced.append(reaction)
        #serial
        #balanced = []
        #for reaction in strings:
        #    screen = self.balance_function(reaction)
        #    if screen:
        #        balanced.append(screen)
        #balanced = [x for x in balanced if x]
        self.reactions = balanced 
        return(balanced)
    
    @staticmethod
    def get_reactants_products(reaction:str) -> list:
            """
            generates a list of dictionaries of reactants and products from a reaction string.
            i.e. '2 H2 + 1 O2 = 2 H2O' :
            [
            {'H2':2,'O2':1}
            {'H2O':2}
            ]
            """
            r,p = reaction.split('=')
            rs = r.split('+')
            rdict = {}
            for i in rs:
                rsplits = [x for x in i.split(' ') if x] 
                rdict[rsplits[-1]] = int(rsplits[0])
            
            ps = p.split('+')
            pdict = {}
            for i in ps:
                psplits = [x for x in i.split(' ') if x]
                pdict[psplits[-1]] = int(psplits[0])    
    
            return([rdict,pdict])
    
    def as_dict(self) -> dict:
        """
        takes the list of reactions and converts it to a dictionary with broken down reactants and products.

        {
        <reaction num.>: {
        'reaction_string: <reaction string>,
        'reactants': <reactants dict>,
        'products': <products dict>
        }
        .
        .
        .
        }
        """
        _dict = {i:{'reaction_string':r} for i,r in enumerate(self.reactions)}
        for i,reaction in _dict.items():
            r,p = self.get_reactants_products(reaction['reaction_string'])
            _dict[i]['reactants'] = r 
            _dict[i]['products'] = p 
        return(_dict)
    
    def to_chempy(self) -> list:
        """
        outputs a list of chempy Equilibrium objects
        REQUIRES chempy to be installed
        """
        try:
            from chempy import Equilibrium 
            return(
                [
                    Equilibrium.from_string(eq) for eq in self.reactions
                    ]
                    )
        except ImportError:
            warnings.warn("chempy not installed - use 'pip install chempy'")
    
    def to_pymatgen(self) -> list:
        """
        outputs a list of pymatgen BalancedReaction objects
        REQUIRES pymatgen to be installed
        """
        try:
            from pymatgen.analysis.reaction_calculator import BalancedReaction 
            return(
                [
                    BalancedReaction.from_dict(r) for r in self.as_dict().values()
                ]
            )
        except ImportError:
            warnings.warn("pymatgen not installed - use 'pip install pymatgen'")

    def to_networkx_graph(self):
        """
        converts the reactions list to a networkx.Graph object 
        """
        try:
            import networkx as nx
            G = nx.Graph()

            for reaction,breakdown in self.as_dict().items():
                G.add_node(reaction) #or just an integer breakdown['reaction_string']
                [
                    G.add_node(molecule) for molecule in list(breakdown['reactants'])+list(breakdown['products'])
                    ]
                [
                    G.add_edge(reactant,reaction) for reactant in list(breakdown['reactants'])
                    ]
                [
                    G.add_edge(reaction,product) for product in list(breakdown['products'])
                    ]
            return(G)
        except ImportError:
            warnings.warn("networkx not installed - use 'pip install networkx'")
