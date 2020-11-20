init2start2fn = {
    'concentrated': {
        'Adaptive mobility regions res=30000': 'concentrated/Adaptive mobility regions/30000',
        'Adaptive mobility regions res=800000': 'concentrated/Adaptive mobility regions/800000',
        'Mobility regions res=1.0': 'Mobility regions/1_0',
        'Mobility regions res=16.0': 'Mobility regions/16_0',
        'Security regions': 'Benchmark divisions/Security regions',
        'Provinces': 'Benchmark divisions/Provinces'
    },
    'evenlydistributed': {
        'Adaptive mobility regions res=200000': 'evenlydistributed/Adaptive mobility regions/200000'
    },
    'historical0310': {
        'Adaptive mobility regions res=12213': 'historical0310/Adaptive mobility regions/12213',
        'Mobility regions res=1.0': 'Mobility regions/1_0',
        'Mobility regions res=16.0': 'Mobility regions/16_0',
        'Security regions': 'Benchmark divisions/Security regions'
    },
    'historical0421': {
        'Adaptive mobility regions res=3000': 'historical0421/Adaptive mobility regions/3000',
        'Mobility regions res=1.0': 'Mobility regions/1_0',
        'Mobility regions res=16.0': 'Mobility regions/16_0',
        'Security regions': 'Benchmark divisions/Security regions'
    }
}
horizon = 30
from time import time
from .comparison import division_path
from .mobility_seir import load_initialization
from .division import Division
from .divisions import ObjectiveScorer,Louvain
from .mobility_seir import  MobilitySEIR
from .mobility import mobility

def optimize_init(init):
    # Obtain initialization
    df = load_initialization(init+'.csv')

    # Find \gamma^*
    min_restr = Division.FromCSV(division_path("Benchmark divisions/Min restrictions.csv"))
    max_restr = Division.FromCSV(division_path("Benchmark divisions/Max restrictions.csv"))
    min_seir, max_seir = (MobilitySEIR(horizon=30, init_df=df, division=div) for div in [min_restr, max_restr])
    min_seir.simulate_all(), max_seir.simulate_all()
    min_restr_infections, max_restr_infections = (seir.cumulative_infections() for seir in [min_seir, max_seir])
    gamma_star = int(round(horizon * mobility.total_weight / (min_restr_infections - max_restr_infections)))
    tradeoffs = [gamma_star,2*gamma_star]

    for save_fn,start_fn in init2start2fn[init].items():
        for tradeoff in tradeoffs:
            div = Division.FromCSV(division_path(start_fn + '.csv'))
            start = div.copy().nextlevel()
            scorer = ObjectiveScorer(df, init_hc=start, tradeoff=tradeoff)
            opt = Louvain(scorer, iterative=True)
            start_time = time()
            opt.optimize(verbose=False)
            end_time = time()
            Division(
                opt.scorer.candidate().relabel_regions(),
                previouslevel=div.copy()
            ).flatClustering().to_csv(division_path('{}/Optimized divisions/{} optimized for tradeoff {}.csv'.format(init,save_fn,tradeoff)))
            print(
                "{fn} resulted in a division with {r} regions and score of {value} and took {t}s".format(
                    fn=save_fn, value=opt.scorer.current_value, t=end_time - start_time,
                    r=len(opt.scorer.HC.regions)))