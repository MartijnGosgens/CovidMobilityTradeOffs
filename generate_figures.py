import matplotlib.pyplot as plt

from .mobility_seir import load_initialization,MobilitySEIR,init_path
from . import visualizations as vis
from . import comparison as comp
import geopandas as gpd
import numpy as np
from . import mezuro_preprocessing as mzp
from .mobility import mobility
from .division import Division
'''
    plot_all() generates all the figures that are used in the paper.
'''

def plot_initialization(ax,init_file,title,vmax=10):
    df = load_initialization(init_file)
    # Compute fraction infections / fraction inhabitants
    infected=df['infected_tested']+df['infected_nottested']
    sum_infective = np.sum(infected)
    population = np.sum(np.sum(df[MobilitySEIR.compartments]))
    fraction_infected = infected/sum_infective
    fraction_inhabitants = sum(df[c] for c in MobilitySEIR.compartments) / population
    df['normalized'] = fraction_infected/fraction_inhabitants
    # add geometry
    df = gpd.GeoDataFrame(df.merge(mzp.gemeente_shapes.reset_index(), left_on='name', right_on ='name'))
    # Plot
    ax.axis('off')
    ax.set_title('%s' % title,fontsize=18)
    df.plot(column = 'normalized',ax=ax,cmap='OrRd',legend=False,vmax=vmax,vmin=0)

# FIGURE 2: concentration mapping
def concentration_mapping():
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(20,10),sharey=True)
    fig.patch.set_facecolor('white')
    vmax = 10

    plot_initialization(ax1, "evenlydistributed.csv", "Evenly distributed", vmax=vmax)
    plot_initialization(ax2, "historical0310.csv", "March 10th, 2020", vmax=vmax)
    plot_initialization(ax3, "concentrated.csv", "Concentrated", vmax=vmax)

    sm = plt.cm.ScalarMappable(cmap='OrRd', norm=plt.Normalize(vmin=0, vmax=vmax))
    sm._A = []
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    cb = fig.colorbar(sm, cax=cbar_ax,orientation='horizontal')
    cb.set_label(label="Active infections density / population density", size=18)
    cb.ax.tick_params(labelsize=18)
    plt.tight_layout()

    plt.savefig('initializations'+'.pdf',bbox_inches='tight', transparent=True)

# FIGURE 3: Plot of the concentration values of the RIVM data over time
def historical_concentrations(load_results=False):
    import matplotlib.dates as mdates
    from datetime import datetime,timedelta
    import pandas as pd
    from .mobility import mobility
    from .rivm_loader import rivm
    from .mobility_seir import MobilitySEIR
    date = datetime(2020, 3, 1)
    concentrations = []
    dates = []
    if load_results:
        df = pd.read_csv(comp.results_path("historical_concentration.csv"),index_col='date')
        dates = [datetime.strptime(d,'%Y-%m-%d') for d in df.index]
        concentrations = list(df['concentration'].values)
    else:
        while date<datetime(2020,11,1):
            try:
                mmdd = date.strftime('%m%d')
                concentrations.append(
                    MobilitySEIR(init_df=rivm.SEI2R2_init(mmdd), mobility=mobility).initial_concentration()
                )
                dates.append(date)
                date = date + timedelta(days=1)
            except:
                break
        # Save results
        df = pd.DataFrame(index=dates)
        df.index.name = 'date'
        df['concentration'] = concentrations
        df.to_csv(comp.results_path("historical_concentration.csv"))

    # Get concentrations of evenly distributed and concentrated
    c_even,c_conc = (
        MobilitySEIR(
            init_df=load_initialization(file),
            mobility=mobility
        ).initial_concentration()
        for file in ["evenlydistributed.csv",
                     "concentrated.csv"]
    )
    months = mdates.MonthLocator()  # every month
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(dates, concentrations, c='b')
    # lines for evenly distributed and concentrated
    ax.plot([dates[0], dates[-1]], [c_even]*2, linestyle='--', c='orange')
    ax.plot([dates[0], dates[-1]], [c_conc]*2, linestyle='--', c='g')
    ax.annotate(xy=(dates[7], c_even+0.01), s="Evenly distributed", c='orange')
    ax.annotate(xy=(dates[7], c_conc-0.02), s="Concentrated", c='g')
    ax.annotate(xy=(dates[7], 0.75), s="Historical (RIVM)", c='b')
    ax.xaxis.set_major_locator(months)
    ax.set_xlabel("Date")
    ax.set_ylabel("Concentration")
    for j,y, date, s in [(17, 0.5, datetime(2020, 3, 10), "March 10"), (35, 0.4, datetime(2020, 4, 21), "April 21")]:
        i = dates.index(date)
        ax.annotate(s, xy=(date, concentrations[i]), xytext=(dates[j], y),
                    arrowprops=dict(arrowstyle="->", mutation_scale=10, linewidth=1))
    plt.savefig('historical_concentrations.pdf', bbox_inches="tight", transparent=True)

# FIGURE 4: Mobility loss: provinces vs mobility regions
def provinces_vs_mobility_regions():
    provinces = Division.FromCSV(comp.division_path("Benchmark divisions/provinces.csv"))
    mobilityreg = Division.FromCSV(comp.division_path("Mobility regions/2_0.csv"))

    fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(20,14),sharey=True)
    fig.patch.set_facecolor('white')

    ax1.axis('off')
    vis.plot_mobility_loss(mobility.G,provinces,ax=ax1,show_legend=False)

    ax2.axis('off')
    vis.plot_mobility_loss(mobility.G,mobilityreg,ax=ax2,show_legend=False)

    # Add marker for Almere
    almere = mobility.G.nodes['Almere']['geometry'].centroid
    ax1.annotate("Almere", xytext=(almere.x-80000, almere.y+30000), xy=(almere.x, almere.y),fontsize=18,
        arrowprops=dict(arrowstyle="->",mutation_scale=20,linewidth=3))

    sm = plt.cm.ScalarMappable(cmap='RdYlGn_r',norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    cb = fig.colorbar(sm, cax=cbar_ax,orientation='horizontal')
    cb.set_label(label= "Fraction mobility lost",size=18)
    cb.ax.tick_params(labelsize=18)

    plt.tight_layout()

    plt.savefig('blocked_mobility'+'.pdf',bbox_inches='tight', transparent=True)

# FIGURE 5: Mobility regions vs adaptive mobility regions (show distribution of active infections?)
def mobility_regions_vs_adaptive_mobility_regions():
    div_amob = Division.FromCSV(comp.division_path("concentrated/Adaptive mobility regions/30000.csv")).relabel_regions()
    div_mob = Division.FromCSV(comp.division_path("Mobility regions/2_0.csv")).relabel_regions()

    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(20,14))
    ax1.axis('off')
    ax2.axis('off')

    uden = mzp.gemeente_shapes.loc['Uden','geometry'].centroid
    ax1.annotate("Uden", xytext=(uden.x-110000, uden.y+110000), xy=(uden.x, uden.y),fontsize=18,
        arrowprops=dict(arrowstyle="->",mutation_scale=20,linewidth=3))

    vmax = 40
    init = load_initialization("concentrated.csv")
    shapes = mzp.gemeente_shapes.copy()
    shapes["cases"] = 10**5 * (init["infected_tested"]+init["infected_nottested"]) / sum(init[c] for c in init.columns)
    for div,ax in [(div_mob,ax1),(div_amob,ax2)]:
        shapes.plot(column="cases",ax=ax,cmap = 'OrRd',legend=False,vmax=vmax,vmin=0)
        vis.plot_division(div,mzp.gemeente_shapes,ax=ax,contours=True)

    sm = plt.cm.ScalarMappable(cmap='OrRd',norm=plt.Normalize(vmin=0, vmax=vmax))
    sm._A = []
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    cb = fig.colorbar(sm, cax=cbar_ax,orientation='horizontal')
    cb.set_label(label= "Active infections per 100.000 inhabitants",size=18)
    cb.ax.tick_params(labelsize=18)
    plt.tight_layout()
    plt.savefig('mobility_vs_adaptive_mobility.pdf',bbox_inches='tight', transparent=True)

from .constants import CoronaConstants
def performances_plot(name,init_file,file_name=None,shown_starts=[],load_performances=False,constants=CoronaConstants,
                      label_x_offsets={},label_y_offsets={}):
    if file_name is None:
        file_name = name+'_comparison'
    fig,ax = plt.subplots(figsize=(8,8))
    if load_performances:
        import json
        results_file_name = name
        # Load file with the right constants
        if len(constants.changed) > 0:
            for var,val in constants.changed.items():
                results_file_name += '_{}_{}'.format(var,val)
        with open(comp.results_path("{}.json".format(results_file_name))) as json_file:
            performances = json.load(json_file)
            performances["optimized_performances"] = {
                tuple(name.split(" optimized for ")): perf
                for name,perf in performances["optimized_performances"].items()
                if name.split(" optimized for ")[0] in shown_starts
            }
            performances["adaptive_mobility_performances"] = {
                float(res): perf
                for res,perf in performances["adaptive_mobility_performances"].items()
            }
            performances["mobility_performances"] = {
                float(res): perf
                for res,perf in performances["mobility_performances"].items()
            }
            comp.plot_comparison(**performances,label_x_offsets=label_x_offsets,label_y_offsets=label_y_offsets,ax=ax)
    else:
        divisions = comp.gather_divisions(name)
        divisions["start_and_tradeoff"] = [
            (start,tradeoff)
            for (start,tradeoff) in divisions["start_and_tradeoff"]
            if start in shown_starts
        ]
        comp.plot_comparison(*comp.perform_comparison(name,init_file,constants=constants,save_results=True,**divisions),
                             label_x_offsets=label_x_offsets,label_y_offsets=label_y_offsets,ax=ax)
    fig.savefig(file_name.replace('.','_')+'.pdf',bbox_inches='tight', transparent=True)



# FIGURE 6a: evenly distributed performance plot
def performances_evenlydistributed(load_performances=True,constants=CoronaConstants):
    performances_plot(name="evenlydistributed",
                      init_file=init_path("evenlydistributed.csv"),
                      shown_starts=['Adaptive mobility regions res=200000'],
                      load_performances=load_performances,
                      constants=constants,
                      label_x_offsets={'Security regions': -8,'Provinces':-5,'Max restrictions': 6,'Min restrictions': -5})
# FIGURE 6b: concentrated performance plot
def performances_concentrated(load_performances=True,constants=CoronaConstants):
    performances_plot(name="concentrated",
                      init_file=init_path("concentrated.csv"),
                      shown_starts=[
                          'Adaptive mobility regions res=800000',
                          'Adaptive mobility regions res=30000',
                          'Mobility regions res=1.0',
                          'Mobility regions res=16.0',
                          'Security regions',
                          'Provinces'
                      ],
                      load_performances=load_performances,
                      constants=constants,
                      label_x_offsets={'Max restrictions': 6,'Min restrictions': -5,'Isolated security region':-6,'Isolated municipality':-8})
# FIGURE 7a: March 10th performance plot
def performances_historical0310(load_performances=True,constants=CoronaConstants):
    performances_plot(name="historical0310",
                      init_file=init_path("historical0310.csv"),
                      shown_starts=[
                          'Adaptive mobility regions res=12213',
                          'Mobility regions res=1.0',
                          'Mobility regions res=16.0',
                          'Security regions'
                      ],
                      load_performances=load_performances,
                      constants=constants,
                      label_x_offsets={'Security regions': -8,'Provinces':-3,'Max restrictions': 6,'Min restrictions': -5,'Isolate Noord-Brabant':-2})
# FIGURE 7b: April 21st performance plot
def performances_historical0421(load_performances=True,constants=CoronaConstants):
    performances_plot(name="historical0421",
                      init_file=init_path("historical0421.csv"),
                      shown_starts=[
                          'Adaptive mobility regions res=3000',
                          'Mobility regions res=1.0',
                          'Mobility regions res=16.0',
                          'Security regions'
                      ],
                      load_performances=load_performances,
                      constants=constants,
                      label_x_offsets={'Security regions': -8,'Provinces':-3,'Max restrictions': 6,'Min restrictions': -5})

def sensitivity_analysis(name,values,symbol,load_performances=True):
    for value in values:
        performances_plot(name="concentrated",
                          init_file=init_path("concentrated.csv"),
                          file_name="concentrated_{}_{}".format(symbol,value).replace('.','_'),
                          shown_starts=[
                              'Adaptive mobility regions res=800000',
                              'Adaptive mobility regions res=30000',
                              'Mobility regions res=1.0',
                              'Mobility regions res=16.0',
                              'Security regions',
                              'Provinces'
                          ],
                          load_performances=load_performances,
                          constants=CoronaConstants(**{name: value}),
                          label_x_offsets={'Max restrictions': 6,'Min restrictions': -5,'Isolated security region':-6,'Isolated municipality':-8})
    
# FIGURE 14a: Robustness with respect to fraction tested people
def sensitivity_analysis_a(a_min=0.01 ,a_max=1, symbol="a",load_performances=True):
    sensitivity_analysis('fraction_tested',[a_min,a_max],symbol,load_performances=load_performances)
# FIGURE 15: Robustness with respect to fraction local contacts
def sensitivity_analysis_p(p_min=0.1 ,p_max=0.9, symbol="p",load_performances=True):
    sensitivity_analysis('fraction_local_contacts',[p_min,p_max],symbol,load_performances=load_performances)
# FIGURE 16: Robustness with respect to effective reproduction number
def sensitivity_analysis_R(R_min=1 ,R_max=2.5, symbol="R",load_performances=True):
    sensitivity_analysis('basic_reproduction_number',[R_min,R_max],symbol,load_performances=load_performances)
# FIGURE 17: Robustness with respect to infectious period
def sensitivity_analysis_w(w_min=3 ,w_max=7.5, symbol="w",load_performances=True):
    sensitivity_analysis('infectious_period',[w_min,w_max],symbol,load_performances=load_performances)
# FIGURE 18: Robustness with respect to latent period
def sensitivity_analysis_v(v_min=3 ,v_max=7.5, symbol="v",load_performances=True):
    sensitivity_analysis('latent_period',[v_min,v_max],symbol,load_performances=load_performances)
    
# FIGURE 19: Effective reproduction number as a function a and p:
    # figure takes long to generate for resolution = 20 (45 minutes), test with resolution=4 (<1 min)
    # resolution = 20 means 21*21 pixels
def effective_reproduction_number(resolution=20):
    from .rivm_loader import rivm

    # No infections initialized.
    init_df = rivm.concentrated_init(0, "Eindhoven")

    grid = []
    for p in range(1,resolution+2):
        row = []
        for alpha in range (1,resolution+2):
            fraction_local_contacts=(p-1)/(resolution)
            fraction_tested=(alpha-1)/(resolution)
            seir = MobilitySEIR(init_df,constants=CoronaConstants(fraction_local_contacts=fraction_local_contacts,fraction_tested=fraction_tested))
            row.append(seir.basic_reproduction_number())
        grid.append(row)
    
    # Plot numpy array:
    fig, ax = plt.subplots(figsize=(8,8))
    im = ax.imshow(grid,cmap='RdYlGn_r')
    ax.set_xticks((np.arange(len(grid))))
    ax.set_yticks((np.arange(len(grid))))
    
    # Manually set ticklabels for resolution of 20:
    if resolution == 20:
        ax.set_xticklabels([0,"","","","",0.25,"","","","",0.5,"","","","",0.75,"","","","",1])
        ax.set_yticklabels([0,"","","","",0.25,"","","","",0.5,"","","","",0.75,"","","","",1])
    if resolution == 4:
        ax.set_xticklabels([0,0.25,0.5,0.75,1])
        ax.set_yticklabels([0,0.25,0.5,0.75,1])

    ax.set_title("Effective reproduction number",fontsize=13)
    fig.colorbar(im,shrink=0.73)

    ax.set_xlabel('fraction tested: a',fontsize=12)
    ax.set_ylabel('fraction local contacts: p',fontsize=12)
    fig.tight_layout()
    fig.savefig("effective_reproduction_number.pdf")


def plot_main(load_performances=True):
    concentration_mapping()
    historical_concentrations()
    provinces_vs_mobility_regions()
    mobility_regions_vs_adaptive_mobility_regions()
    performances_evenlydistributed(load_performances=load_performances)
    performances_concentrated(load_performances=load_performances)
    performances_historical0310(load_performances=load_performances)
    performances_historical0421(load_performances=load_performances)

def plot_supplementary(load_performances=True):
    sensitivity_analysis_a(load_performances=load_performances)
    sensitivity_analysis_p(load_performances=load_performances)
    sensitivity_analysis_R(load_performances=load_performances)
    sensitivity_analysis_w(load_performances=load_performances)
    sensitivity_analysis_v(load_performances=load_performances)
    effective_reproduction_number()

def plot_all(load_performances=True):
    plot_main(load_performances=load_performances)
    plot_supplementary(load_performances=load_performances)