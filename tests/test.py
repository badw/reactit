from re import L
import pytest 
from reactit import ReactionGenerator

compounds = {0:'CO2',1:'CO',3:'H2O',4:'H2'}

def test_convert_number_strings():

    string_reactions = [
        '[[0,1],[2,3]]',
    '[[4,5,6],[7]]']

    real_reactions = [((0, 1), (2, 3)), ((4, 5, 6), (7,))]
    
    rg = ReactionGenerator(compounds = compounds)
    converted_string_reactions = list(
        rg.convert_number_strings(string_reactions)
    )
    assert converted_string_reactions == real_reactions

def test_numeric_reaction_filter():

    reaction1 = [[0,1],[2,3]] #normal
    reaction2 = [[0,1],[0,1]] #same
    reaction3 = [[0,1,2],[2,1,0]] # mirror 

    rg = ReactionGenerator(compounds = compounds)
    reaction1_filter = rg.numeric_reaction_filter(reaction1)
    reaction2_filter = rg.numeric_reaction_filter(reaction2)
    reaction3_filter = rg.numeric_reaction_filter(reaction3)

    assert reaction1_filter == str(reaction1)
    assert reaction2_filter is None
    assert reaction3_filter is None

def test_string_reaction_filter():

    reaction1 = [['N2','O2','H2'],['HNO3']] # should be fine 
    reaction2 = [['N2','H2'],['CH4']] # no C in reactants, not N in products

    rg = ReactionGenerator(compounds = compounds)
    reaction1_filter = rg.string_reaction_filter(reaction1)
    reaction2_filter = rg.string_reaction_filter(reaction2)

    assert reaction1_filter == reaction1
    assert reaction2_filter is None 

def test_enumerate_combinations():

    real_max_three = [((0,), (1, 2)),
                      ((0,), (1, 3)),
                      ((0,), (2, 3)),
                      ((0, 1), (2,)),
                      ((0, 1), (3,)),
                      ((0, 2), (1,)),
                      ((0, 2), (3,)),
                      ((0, 3), (1,)),
                      ((0, 3), (2,)),
                      ((1,), (2, 3)),
                      ((1, 2), (3,)),
                      ((1, 3), (2,))]
    
    real_max_four = [((0,), (1, 2)),
                     ((0,), (1, 2, 3)),
                     ((0,), (1, 3)),
                     ((0,), (2, 3)),
                     ((0, 1), (2,)),
                     ((0, 1), (2, 3)),
                     ((0, 1), (3,)),
                     ((0, 1, 2), (3,)),
                     ((0, 1, 3), (2,)),
                     ((0, 2), (1,)),
                     ((0, 2), (1, 3)),
                     ((0, 2), (3,)),
                     ((0, 2, 3), (1,)),
                     ((0, 3), (1,)),
                     ((0, 3), (1, 2)),
                     ((0, 3), (2,)),
                     ((1,), (2, 3)),
                     ((1, 2), (3,)),
                     ((1, 3), (2,))]

    rg = ReactionGenerator(compounds=compounds)
    
    max_two = sorted(
        list(
            rg.enumerate_combinations(max_length = 2)
            )
            )
    max_three = sorted(
        list(
            rg.enumerate_combinations(max_length=3)
        )
    )
    max_four = sorted(
        list(
            rg.enumerate_combinations(max_length=4)
        )
    )

    assert max_two == []
    assert max_three == real_max_three
    assert max_four == real_max_four

    