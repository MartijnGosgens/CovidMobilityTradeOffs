# Parameter choices of obtained divisions.
mobres = [0.15, 0.25, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0, 8.0, 16.0, 32.0]
init2admobres = {
    'concentrated': [
        200,
        3000,
        4000,
        8000,
        12000,
        30000,
        100000,
        200000,
        800000
    ],
    'evenlydistributed': [
        3000,
        4000,
        8000,
        12000,
        20000,
        30000,
        40000,
        60000,
        100000,
        200000,
        400000
    ],
    'historical0310': [
        350,
        373,
        519,
        1488,
        2519,
        4263,
        7216,
        12213,
        822413
    ],
    'historical0421': [
        30,
        40,
        60,
        100,
        200,
        400,
        600,
        1200,
        1600,
        2000,
        3000,
        4000,
        8000
    ]
}
from .comparison import division_path
from .mobility_seir import load_initialization
from .divisions import mobility_regions,adaptive_mobility_regions

def compute_mobility_regions():
    res2mobility_regions = {}
    for res in mobres:
        res2mobility_regions[res] = mobility_regions(res=res).relabel_regions()
        res2mobility_regions[res].to_csv(division_path('Mobility regions/'+str(res).replace('.','_')+'.csv'))

def compute_adaptive_mobility_regions():
    init2res2adaptive_mobility_regions = {}
    for init,ress in init2admobres.items():
        init2res2adaptive_mobility_regions[init] = {}
        df = load_initialization(init+'.csv')
        for res in ress:
            admobreg = adaptive_mobility_regions(init_df=df,res=res).relabel_regions()
            init2res2adaptive_mobility_regions[init][res] = admobreg
            admobreg.to_csv(division_path('{}/Adaptive mobility regions/{}.csv'.format(init,res)))