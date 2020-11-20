import networkx as nx

'''
    When running

        from .networks import avg_graph
        G = avg_graph(gemeente_shapes, gemeente2gemeente)

    The resulting networkx.DiGraph has nodes representing gemeentes and
    edges representing the average mobility (frequent+regelmatig+incidenteel) per day.
'''

'''
    Takes a GeoDataFrame selected_shapes (e.g. gemeente_shapes) and a mobility
    DataFrame shape2shape (e.g. gemeente2gemeente)
    and constructs a networkx.DiGraph from it.
    The nodes of this network will correspond to the index column of selected_shapes
    while the weights will correspond to the values in column col of shape2shape
    corresponding to date.
    The node attributes will be copied from the columns of selected_shapes.

    Parameters
    date: the date for which the mobility should be taken
    selected_shapes: a GeoDataFrame that is indexed by the shapes which we would
        like to use as nodes.
    shape2shape: a DataFrame indexed a pair of keys of selected_shapes and a date.
'''
def graph_of_date(date,selected_shapes, shape2shape,col='totaal_aantal_bezoekers'):
    shape2shape_date = shape2shape.groupby('datum').get_group(date)
    G = nx.DiGraph()
    G.add_nodes_from(selected_shapes.index)
    nx.set_node_attributes(G, values=selected_shapes.to_dict('index'))
    G.add_weighted_edges_from([
        (idx[0],idx[1],row[col])
        for idx,row in shape2shape_date.iterrows()
        if idx[0] in selected_shapes.index and idx[1] in selected_shapes.index
    ])
    return G

'''
    Takes a GeoDataFrame selected_shapes (e.g. gemeente_shapes) and a mobility
    DataFrame shape2shape (e.g. gemeente2gemeente)
    and constructs a networkx.DiGraph from it.
    The nodes of this network will correspond to the index column of selected_shapes
    while the weights will correspond to the values in column col of shape2shape
    (averaged over all dates).
    The node attributes will be copied from the columns of selected_shapes.

    Parameters
    selected_shapes: a GeoDataFrame that is indexed by the shapes which we would
        like to use as nodes.
    shape2shape: a DataFrame indexed a pair of keys of selected_shapes and a date.
'''
def avg_graph(selected_shapes, shape2shape,col='totaal_aantal_bezoekers'):
    shape2shape_agg = shape2shape.groupby(['woon','bezoek']).agg({col: 'mean'})
    G = nx.DiGraph()
    G.add_nodes_from(selected_shapes.index)
    nx.set_node_attributes(G, values=selected_shapes.to_dict('index'))
    G.add_weighted_edges_from([
        (idx[0],idx[1],row[col])
        for idx,row in shape2shape_agg.iterrows()
        if idx[0] in selected_shapes.index and idx[1] in selected_shapes.index
    ])
    return G

def traffic_blocked(node,G,division):
    if type(node) in (set, list):
        return sum(traffic_blocked(i,G,division) for i in node)
    else:
        return sum([
            G[node][j]['weight']
            for j in G.neighbors(node)
            if division[node]!=division[j]
        ])
def total_traffic(node,G):
    if type(node) in (set, list):
        return sum(total_traffic(i,G) for i in node)
    else:
        return sum([
            G[node][j]['weight']
            for j in G.neighbors(node)
        ])
def fraction_traffic_blocked(node,G,division):
    total = total_traffic(node,G)
    if total==0:
        return 0
    else:
        return traffic_blocked(node,G,division)/total
