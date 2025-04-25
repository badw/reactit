import numpy as np 
import itertools as it 
import tqdm 
import gzip 
import os 
import math 
from chempy import balance_stoichiometry, Equilibrium
import warnings 


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
        for reaction_length in range(2,max_length+1):
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
        indexes =  tuple(self.compounds)
        approved = []
        for r in reaction_indexes:
            if not [x for x in it.chain(*r) if x not in indexes]:
                approved.append(r)
        return(approved)
    
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
    
    def convert_to_equation(self,converted_strings):
        converted = []
        for r in tqdm.tqdm(converted_strings):
            re,pr = r
            try:
                converted.append(balance_stoichiometry(list(re),list(pr),underdetermined=None))
            except Exception:
                try:
                    converted.append(balance_stoichiometry(list(re),list(pr)))
                except Exception:
                    pass                            
        return(converted) 
    
    def screen_converted(self,converted_reactions): # this should be multiprocessed - > perhaps a mp.Pool? as it only needs the big list
        def _convert_ord_to_dict(r):
            re,pr =  r
            try:
                reacs = {k:int(re[k]) for k in re}
                prods = {k:int(pr[k]) for k in pr}
            except Exception:
                warnings.warn('\n error with {}'.format(r))
            return(reacs,prods)
    
        screened = []
        for r in converted_reactions:
            try:
                re,pr = _convert_ord_to_dict(r)
                try:
                    screened.append(Equilibrium(re,pr))
                except Exception:
                    pass
            except Exception :
                pass
        return(screened)    
    
    def run_all(self):
        warnings.filterwarnings('ignore')
        print('orig =',len(self.reactions),end='...')
        approved = self.remove_indices(self.reactions)
        print(' approved = ',len(approved),end='...')
        strings = self.convert_to_string(approved)
        print(' prescreening = ',len(strings))
        equations = self.convert_to_equation(strings)
        screened = self.screen_converted(equations)
        print(' final = ',len(screened),end='...')
        return(screened)