from .mobility import mobility
from .constants import CoronaConstants
import os
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from .division import Division

'''
    Usage:
        from COVID_tools.comparison import perform_comparison,gather_divisions,plot_comparison
        name = "concentrated"
        init_file = "COVID_tools/initializations/concentrated_after10days_deterministicsimulation.csv"
        mobility_performances,optimized_performances,benchmark_performances = perform_comparison(name,init_file,**gather_divisions(name))
        plot_comparison(mobility_performances,optimized_performances,benchmark_performances)
    Comments:
        different divisions can be chosen by supplying a name2division dict to
        init_seir_models().
'''
def division_path(division_name):
    return os.path.join(
        os.path.dirname(__file__),
        'divisions',
        division_name
    )


def results_path(name):
    return os.path.join(
        os.path.dirname(__file__),
        'results',
        name
    )


# Select divisions
division_names = [
    "Min restrictions",
    "Max restrictions",
    "Provinces",
    "Security regions",
]

name2filename = {
    name: name+'.csv'
    for name in division_names
}
'''
    Returns a list of all the adaptive mobility resolution values for which there is a
    division in divisions/[name]/Adaptive mobility regions/.
'''
def res_in_folder(folder=None,name=None):
    import os
    if folder==None:
        folder = division_path(name+'/Adaptive mobility regions/')
    ress = [
        float(f.split('.')[0].replace('_','.')) if '_' in f else int(f.split('.')[0])
        for f in os.listdir(folder)
    ]
    ress.sort()
    return ress

def get_divisions(name2filename=name2filename):
    return {
        name: Division.FromCSV(division_path(filename))
        for name,filename in name2filename.items()
    }

def init_seir_models(init_df,horizon=30,name2division=None, mobility=mobility, constants = CoronaConstants):
    from .mobility_seir import MobilitySEIR
    if name2division == None:
        name2division = get_divisions()
    return {
        name: MobilitySEIR(constants = constants, horizon=horizon, init_df=init_df, division=division, mobility=mobility)
        for name,division in name2division.items()
    }

def perform_simulations(seirs):
    for seir in seirs:
        seir.simulate_all()

def get_performances(name2seir,t=30,total_pop=None,total_mob=None):
    if total_pop is None:
        from .mezuro_preprocessing import gemeente_shapes
        total_pop = gemeente_shapes.inhabitant.sum()
    if total_mob is None:
        from .mobility import mobility
        total_mob = mobility.total_weight
    return {
        name: {
            "Mobility": 100*mobility.intra_region_mobility(seir.division)/total_mob,
            "Infections": 100000*seir.cumulative_infections(t=t)/total_pop
        }
        for name, seir in name2seir.items()
    }

def get_intervals(performances):
    from collections import OrderedDict
    perf = OrderedDict(sorted(performances.items(), key=lambda kv: kv[1]['Infections']))

    min_gamma_A_betterthan_B = {
        (A,B): (
            perf[B]['Mobility']-perf[A]['Mobility']
        ) / (
            perf[B]['Infections']-perf[A]['Infections']
        )
        for A,B in it.combinations(perf.keys(),2)
    }

    A2B2thresshold = {}
    for pair,thresshold in min_gamma_A_betterthan_B.items():
        if not pair[0] in A2B2thresshold:
            A2B2thresshold[pair[0]] = {}
        A2B2thresshold[pair[0]][pair[1]] = thresshold

    # Figure out the starting value of the interval where A is best
    intervalstart2division = {}
    min_gamma = float('inf')
    for A,B2thresshold in A2B2thresshold.items():
        highest = max(B2thresshold.values())
        if highest > min_gamma:
            continue
        intervalstart2division[highest] = A
        min_gamma = highest
    return intervalstart2division

'''
    Plots a line for the performances of a heuristic for various res values.
'''
def plot_heuristic(res2performance,ax,marker='s',color='b',label=None,x_offset=0,y_offset=0,res_labels=[]):
    from collections import OrderedDict
    perf = OrderedDict(sorted(res2performance.items(), key=lambda kv: kv[1]['Mobility']))
    x = [v['Mobility'] for v in perf.values()]
    y = [v['Infections'] for v in perf.values()]

    ax.plot(x,y,marker=marker,linestyle=':',label=label,color=color)

    for res in res_labels:
        p = perf[res]
        ax.text(p['Mobility']+x_offset,p['Infections']+y_offset,res,horizontalalignment='center')

'''
    Gives the names of the divisions for an initialization with name name.
'''
def gather_divisions(name):
    import os
    mobility_ress = res_in_folder(folder=division_path("Mobility regions/"))
    adaptive_mobility_ress = res_in_folder(name=name)
    benchmarks = [
        "Min restrictions",
        "Max restrictions",
        "Provinces",
        "Security regions"
    ]
    custom_benchmarks = []
    custom_benchmark_dir = division_path(name+"/Benchmark divisions/")
    optimized_dir = division_path(name+"/Optimized divisions/")
    start_and_tradeoff = []
    if os.path.isdir(custom_benchmark_dir):
        custom_benchmarks = [
            f.split(".")[0]
            for f in os.listdir(custom_benchmark_dir)
        ]
    if os.path.isdir(optimized_dir):
        start_and_tradeoff = [
            f[:-4].split(" optimized for tradeoff ")
            for f in os.listdir(optimized_dir)
        ]
    return {
        "mobility_ress": mobility_ress,
        "adaptive_mobility_ress": adaptive_mobility_ress,
        "start_and_tradeoff": start_and_tradeoff,
        "benchmarks": benchmarks,
        "custom_benchmarks": custom_benchmarks,
    }

'''
    Obtains the results for the given benchmark and heuristic divisions.
'''
def perform_comparison(
        name, init_file,
        mobility_ress=[],
        adaptive_mobility_ress=[],
        start_and_tradeoff=[],
        benchmarks=[],
        custom_benchmarks=[],
        horizon=30,
        save_results=False,
        constants=CoronaConstants):
    # Load df
    init_df = pd.read_csv(init_file,index_col='name')

    # Obtain file names
    mobility_res2file = {
        res: 'Mobility regions/'+str(res).replace('.','_')+'.csv'
        for res in mobility_ress
    }
    adaptive_mobility_res2file = {
        res: name+'/Adaptive mobility regions/'+str(res).replace('.','_')+'.csv'
        for res in adaptive_mobility_ress
    }
    start_and_tradeoff2file = {
        (start,tradeoff): name+'/Optimized divisions/{} optimized for tradeoff {}.csv'.format(start,tradeoff)
        for start,tradeoff in start_and_tradeoff
    }
    benchmark2file = {
        bench: 'Benchmark divisions/'+bench+'.csv'
        for bench in benchmarks
    }
    custombenchmark2file = {
        bench: name+'/Benchmark divisions/'+bench+'.csv'
        for bench in custom_benchmarks
    }
    # Obtain divisions
    res2mobility = get_divisions(mobility_res2file)
    adaptive_res2mobility = get_divisions(adaptive_mobility_res2file)
    start_tradeoff2optimized = get_divisions(start_and_tradeoff2file)
    benchmarks = get_divisions(benchmark2file)
    benchmarks.update(get_divisions(custombenchmark2file))
    # Initialize seir models
    mobility_seirs = init_seir_models(name2division=res2mobility,init_df=init_df,constants=constants)
    adaptive_mobility_seirs = init_seir_models(name2division=adaptive_res2mobility,init_df=init_df,constants=constants)
    optimized_seirs = init_seir_models(name2division=start_tradeoff2optimized,init_df=init_df,constants=constants)
    benchmark_seirs = init_seir_models(name2division=benchmarks,init_df=init_df,constants=constants)
    # Run
    perform_simulations(mobility_seirs.values())
    perform_simulations(adaptive_mobility_seirs.values())
    perform_simulations(optimized_seirs.values())
    perform_simulations(benchmark_seirs.values())
    # Obtain performances
    mobility_performances = get_performances(mobility_seirs,t=horizon)
    adaptive_mobility_performances = get_performances(adaptive_mobility_seirs,t=horizon)
    optimized_performances = get_performances(optimized_seirs,t=horizon)
    benchmark_performances = get_performances(benchmark_seirs,t=horizon)
    if save_results:
        import json
        file_name = name
        # Make sure that we don't override results from different constants
        if len(constants.changed) > 0:
            for var,val in constants.changed.items():
                file_name += '_{}_{}'.format(var,val)
        with open(results_path("{}.json".format(file_name)), 'w') as outfile:
            json.dump({
                "mobility_performances": mobility_performances,
                "adaptive_mobility_performances": adaptive_mobility_performances,
                "optimized_performances": {
                    "{} optimized for {}".format(start,tradeoff): perf
                    for (start,tradeoff),perf in optimized_performances.items()
                },
                "benchmark_performances": benchmark_performances
            }, outfile, indent=4)
    return (mobility_performances,adaptive_mobility_performances,optimized_performances,benchmark_performances)

'''
    Places markers for each of the benchmark_performances and
    draws lines for the heuristics.
'''
def plot_comparison(mobility_performances,adaptive_mobility_performances,
                    optimized_performances,benchmark_performances,
                    label_x_offsets={},label_y_offsets={},ax=None):
    ax = plot_division_performances(benchmark_performances,exp_fit=False,ax=ax,label_x_offsets=label_x_offsets,label_y_offsets=label_y_offsets)
    plot_heuristic(mobility_performances,ax=ax,label='Mobility regions',color='g',x_offset=-0.2*10**8)
    plot_heuristic(adaptive_mobility_performances,ax=ax,marker='v',label='Adaptive mobility regions',color='r',x_offset=+0.2*10**8)

    colors = ['c','m','y','brown','orange']
    tradeoffs = {tradeoff for _,tradeoff in optimized_performances.keys()}
    tradeoff2c = {tradeoff: colors[i] for i,tradeoff in enumerate(tradeoffs)}

    mob_prefix = "Mobility regions res="
    admob_prefix = "Adaptive mobility regions res="
    for (start,tradeoff),perf in optimized_performances.items():
        perf_start = {}
        if mob_prefix in start:
            perf_start = mobility_performances[float(start[len(mob_prefix):])]
        elif admob_prefix in start:
            perf_start = adaptive_mobility_performances[float(start[len(admob_prefix):])]
        else:
            perf_start = benchmark_performances[start]
        ax.annotate("", xy=(perf["Mobility"], perf["Infections"]), xytext=(perf_start["Mobility"], perf_start["Infections"]),
            arrowprops=dict(arrowstyle="-|>",color=tradeoff2c[tradeoff],mutation_scale=15,linewidth=2))
        #ax.text(perf["Mobility"],perf["Infections"],str(tradeoff),horizontalalignment='left')
    # Make legend with entries for the heuristics and add entries
    # for the tradeoffs
    leg = ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    for tradeoff,c in tradeoff2c.items():
        handles.append(Line2D([0], [0], color=c))
        labels.append(r"Optimized for $\gamma={}$".format(tradeoff))
    ax.legend(handles,labels)
    ax.set_xlabel("Percentage mobility allowed",fontsize='large')
    ax.set_ylabel("Infections per 100.000",fontsize='large')
    return ax

def plot_division_performances(performances,ax=None, exp_fit=False,label_y_offsets={},label_x_offsets={}):
    if ax==None:
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
    results_df = pd.DataFrame.from_dict(performances,orient='index')
    results_df.plot.scatter(x='Mobility',y='Infections',ax=ax,marker='x',s=75)

    y_offset = (results_df['Infections'].max()-results_df['Infections'].min())/70

    for name, point in results_df.iterrows():
        y = point['Infections']+y_offset
        x = point['Mobility']
        if name in label_y_offsets:
            y += label_y_offsets[name]*y_offset
        if name in label_x_offsets:
            x += label_x_offsets[name]
        ax.text(x, y, name,horizontalalignment='center')

    if exp_fit:
        ax.plot(*exp_fit_xy(performances))
    return ax

def fit_exponent(performance):
    return np.log(
        performance['Min restrictions']['Infections']/performance['Max restrictions']['Infections']
    ) /(performance['Min restrictions']['Mobility'] - performance['Max restrictions']['Mobility'])

def exp_fit(performance):
    exponent=fit_exponent(performance)
    min_infections = performance['Max restrictions']['Infections']
    return lambda mobility: min_infections*np.exp(exponent*mobility)

def exp_fit_xy(performance):
    from collections import OrderedDict
    ordered = OrderedDict(sorted(performance.items(), key=lambda kv: kv[1]['Mobility']))
    fit = exp_fit(performance)
    xs = [v['Mobility'] for v in ordered.values()]
    ys = [fit(x) for x in xs]
    return xs,ys
