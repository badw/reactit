import numpy as np 
import itertools as it 
import tqdm 
import gzip 
import os 
import math 
import warnings 
from sympy import symbols, Eq, solve    
import functools 
import math 
from fractions import Fraction 
import re 
from collections import defaultdict 
import tqdm_pathos 



class ReactionsListGenerator:
    '''
    `ReactionsDictionaryGenerator` is a class that generates a dictionary of values:
    ([0], [1, 2, 3])
([0], [1, 2, 4])
([0], [1, 2, 5])
([0], [1, 2, 6])
([0], [1, 2, 7])
based on a given set of placemarker values
'''
    
    def __init__(
            self,
            no_compounds=None,
            path='.'
            ):
        self.nc = no_compounds
        self.path = path

    def convert_strings(self,reactions):
        for reaction in reactions:
            r,p = reaction.strip().split('],')
            r = tuple(int(x) for x in r.split('[[')[1].split(',') if x)
            p = tuple(int(x) for x in p.split('[')[1].split(']]')[0].split(',') if x)
            yield ((r,p))

    @staticmethod    
    def reaction_filter(reaction):
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

    def iterate(
            self,
            max_length=4,
    ):
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
                        x for x in range(self.nc+1)
                    ],
                    size[0]
                )    
    

                products = it.combinations(
                    [
                        x for x in range(self.nc+1)
                    ],
                    size[1]
                )
                combined = tuple(
                    it.product(reactants, products)
                    )
                filtered = []
                for v in map(lambda x: ReactionsListGenerator.reaction_filter(x), combined):
                    if v:
                        filtered.append(v)
                filtered = list(set(filtered))
                reactions.extend(filtered)

        reactions = self.convert_strings(reactions)
        return(reactions)

class MappingtoReaction:
    ''' a class that takes a premade reactions reference dictionary and generates a list of reactions with it that are then further balanced and filtered'''
    
    def __init__(self,reactions,compounds):
        self.reactions = reactions
        if isinstance(compounds,list):
            self.compounds = {i:c for i,c in enumerate(compounds)}
        else:
            self.compounds = compounds            
            
    def remove_indices(self,reaction_indexes):
        ''''''
        @staticmethod
        def filterfunc(iterable,valid_nums):
            spec = list(it.chain(*iterable))
            if not [x for x in spec if x not in valid_nums]:
                return(iterable)
                    

        filtered = list(
            filter(
                lambda x: x is not None, map(
                    lambda x: filterfunc(x,self.compounds.keys()),self.reactions
                    )
                    )
        )
        return(filtered)
    
    def convert_to_string(self,approved_list):
    
        def _screen_string(r):
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
    
        converted = []
        for r in approved_list:
            d = [[self.compounds[i] for i in r[0]],[self.compounds[i] for i in r[1]]]
            ds = _screen_string(d)
            if ds:
                converted.append(ds)
        return(converted)
    
    @staticmethod
    def parse_molecule(formula): 
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
    def find_int(_lst):
        denoms = [Fraction(x).denominator for x in _lst]
        return (
            functools.reduce(
                lambda a, b: a*b//math.gcd(a, b), denoms
            )
        )    

    def balance_reaction(self,reactants, products):
        '''this function allows for undetermined systems
        '''
        # Create a list of all elements
        product_keys = list(products.keys())
        reactant_keys = list(reactants.keys())
        elements = set()
        for species in (reactant_keys + product_keys):
            elements.update(species)
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
                common_factor = self.find_int(solution[0].values())
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
                    rp == 'p'
                    if ps == 0:
                        product_string += f'{coeff} {product_keys[int(num)]}'
                    else:
                        product_string += f' + {coeff} {product_keys[int(num)]}'
                    ps += 1
            equation_string = reactant_string + ' = ' + product_string
        return(equation_string) 
    
    def balance_function(self,iterable):
            '''balance_function checks the balance of a given list of strings, this is here for future multiprocessing
            '''
            reactants = {
                molecule:self.parse_molecule(molecule) for molecule in iterable[0]
                }
            products = {
                molecule:self.parse_molecule(molecule) for molecule in iterable[1]
                }
            balanced = self.balance_reaction(reactants,products)
            if balanced:
                return(balanced)
    
    def run_all(self):
        warnings.filterwarnings('ignore')
        approved = self.remove_indices(self.reactions)
        strings = self.convert_to_string(approved)
        #screened = tqdm_pathos.map(self.mp_function,strings)
        screened = []
        for reaction in tqdm.tqdm(strings):
            screen = self.balance_function(reaction)
            if screen:
                screened.append(screen)
        screened = [x for x in screened if x]
        self.screened = screened 
        return(screened)
    
    @staticmethod
    def get_reactants_products(reaction):
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
    
    def as_dict(self):
        _dict = {i:{'reaction_string':r} for i,r in enumerate(self.screened)}
        for i,reaction in _dict.items():
            r,p = self.get_reactants_products(reaction['reaction_string'])
            _dict[i]['reactants'] = r 
            _dict[i]['products'] = p 
        return(_dict)