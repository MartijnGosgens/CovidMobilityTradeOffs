
import geopandas as gdp
import matplotlib.pyplot as plt
from matplotlib import cm

from .networks import fraction_traffic_blocked

'''
    Plots the division of the shapes on ax.

    PARAMETERS
    division: The division of the shapes into regions, a dictionary from area to region.
    shapes: The GeoDataFrame that will be plotted.
    ax: The axes on which the map will be plotted. If None, then a new axis will
        be created.
    contours: If True, then only the borders of the regions will be drawn, leaving
        the areas transparent. Else, the areas will be coloured using 20 random colors.
'''
def plot_division(division,shapes,ax=None,contours=False):
    if ax==None:
        fig, ax = plt.subplots(1, 1, figsize=(15,15))
        ax.axis('off')
    shapes = shapes.copy()
    shapes['region'] = [
        division[shape]
        for shape in shapes.index
    ]
    if contours:
        community_shapes = shapes.dissolve(by='region')
        community_shapes.plot(ax=ax,linewidth=2,facecolor="none",edgecolor='black')
    else:
        shapes.plot(column='region',cmap=cm.tab20,ax=ax)
    return ax

'''
    Given a dict from initialization name to concentration value, this plots
    these numbers on a 1d axis with labels at the values.
'''
def plot_concentrations(name2concentration):
    fig,ax = plt.subplots(figsize=(9,1))
    # Remove lines but the x axis
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    # Move x axis to the middle of the figure
    ax.spines["bottom"].set_position(("data", 0.5))
    ax.axes.get_yaxis().set_visible(False)
    plt.xlim(0,1)

    # I don't understand a single line of the code below but it does what it needs to do
    labs, vals = zip(*sorted([(k,v) for k,v in name2concentration.items()], key=lambda t: t[1] ))
    ticks = ["{:.2f}\n{}".format(v,k) for k,v in zip(labs,vals)]
    ax.set_xticks(vals[::1])
    ax.set_xticklabels(ticks[::1])
    #ax.set_xticks(vals[1::1], minor=True)
    #ax.set_xticklabels(ticks[1::2], minor=True, va="bottom")
    ax.tick_params(which="minor", direction="in", pad=-10 )

    # This draws the arrowheads at the end of the axis. Don't ask me how.
    ax.annotate("", xy=(1,0.5), xycoords="axes fraction", xytext=(5,0),textcoords="offset points", ha="center",
                arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0,facecolor='black'))
    ax.annotate("", xy=(0,0.5), xycoords="axes fraction", xytext=(-5,0),textcoords="offset points", ha="center",
                arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0,facecolor='black'))
    return ax

'''
    Plots the loss in mobility resulting from a division of the shapes.

    PARAMETERS
    G: The graph that contains the mobility between the shapes.
    division: The division of the nodes/shapes into regions, a dictionary from area to region.
    shapes: The GeoDataFrame that will be plotted. Not needed if each node of G
        has a 'geometry' attribute.
    ax: The axes on which the map will be plotted. If None, then a new axis will
        be created.
    legend_kwds: optional dictionary containing additional keywords for the legend of the plot.
'''
def plot_mobility_loss(G,division,ax=None,show_legend=True,legend_kwds={},shapes=None):
    if not shapes:
        shapes = gdp.GeoDataFrame.from_dict(G.nodes,orient='index')
    if ax==None:
        fig, ax = plt.subplots(1, 1, figsize=(15,15))
        ax.axis('off')
    shapes = shapes.copy()
    shapes['fraction_blocked'] = [
        fraction_traffic_blocked(shape,G,division)
        for shape in shapes.index
    ]
    if not 'label' in legend_kwds:
        legend_kwds['label'] = "Fraction mobility lost"
    shapes.plot(column='fraction_blocked',ax=ax,legend=show_legend,legend_kwds=legend_kwds,cmap='RdYlGn_r',vmax=1,vmin=0)
    plot_division(division,shapes,ax=ax,contours=True)
