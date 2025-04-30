from reactit import ReactionGenerator


def test_class_initialisation():
    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    compounds_list = ['CO2','CO','H2O','H2']

    #list
    rg = ReactionGenerator(compounds=compounds_list)
    assert rg.compounds == compounds
    #dict
    rg = ReactionGenerator(compounds=compounds)
    assert rg.compounds ==compounds

def test_convert_number_strings():

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

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

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

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

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    reaction1 = [['N2','O2','H2'],['HNO3']] # should be fine 
    reaction2 = [['N2','H2'],['CH4']] # no C in reactants, not N in products

    rg = ReactionGenerator(compounds = compounds)
    reaction1_filter = rg.string_reaction_filter(reaction1)
    reaction2_filter = rg.string_reaction_filter(reaction2)

    assert reaction1_filter == reaction1
    assert reaction2_filter is None 

def test_enumerate_combinations():

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

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

def test_convert_to_string():
    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    rg = ReactionGenerator(compounds = compounds)

    result = rg.convert_to_string([((1, 2), (0, 3))])

    assert result == [[['CO', 'H2O'], ['CO2', 'H2']]]

def test_parse_molecule():
    
    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    rg = ReactionGenerator(compounds = compounds)

    molecule = rg.parse_molecule('CO2')
    
    assert molecule == {'C': 1, 'O': 2}

def test_find_coefficients():

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    test_coefficients = [1,0.5,1] # H2 + 0.5 O2 = H2O

    rg = ReactionGenerator(compounds = compounds)

    coeffs = rg.find_coefficients(fractional_coefficients = test_coefficients)

    assert coeffs == 2

def test_balance_reaction():

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    reactants = {'H2':{'H':2},'O2':{'O':2}}
    products = {'H2O':{'H':2,'O':1}}

    rg = ReactionGenerator(compounds = compounds)

    balanced_reaction = rg.balance_reaction(reactants,products)

    assert balanced_reaction == '2 H2 + 1 O2 = 2 H2O'


def test_balance_function():
    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2'}

    reaction_list = [[['H2O'], ['H2', 'O2']]]

    rg = ReactionGenerator(compounds = compounds) 

    balanced_reaction = rg.balance_function(reaction_list[0])

    assert balanced_reaction == '2 H2O = 2 H2 + 1 O2'


def test_iterate():

    compounds = {0:'CO2',1:'CO',2:'H2O',3:'H2',4:'NO2',5:'N2',6:'O2'}

    rg = ReactionGenerator(compounds = compounds) 

    three_reaction = sorted(rg.iterate(3))
    four_reaction = sorted(rg.iterate(4))

    assert three_reaction == ['2 CO2 = 2 CO + 1 O2', '2 H2O = 2 H2 + 1 O2', '2 NO2 = 1 N2 + 2 O2']

    assert four_reaction == ['1 CO2 = 1 CO + 1 H2O',
                             '2 CO2 = 2 CO + 1 NO2',
                             '2 CO2 = 2 CO + 1 O2',
                             '2 H2O = 2 H2 + 1 NO2',
                             '2 H2O = 2 H2 + 1 O2',
                             '2 NO2 = 1 N2 + 2 O2']
    

def test_get_reactants_products():

    compounds = {0:'O2',1:'H2O',2:'H2'}

    rg = ReactionGenerator(compounds = compounds)

    _dict = rg.get_reactants_products('2 H2O = 2 H2 + 1 O2')

    assert _dict == [{'H2O': 2}, {'H2': 2, 'O2': 1}]


def test_as_dict():
    compounds = {0:'O2',1:'H2O',2:'H2'}

    rg = ReactionGenerator(compounds = compounds)

    rg.iterate(4)

    _dict = rg.as_dict()

    assert _dict == {0: {'reaction_string': '1 O2 + 2 H2 = 2 H2O',
                         'reactants': {'O2': 1, 'H2': 2},
                         'products': {'H2O': 2}}}








    
    


    