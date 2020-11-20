from enum import Enum
from collections import deque
import math
import networkx as nx
import numpy as np
import itertools as it
from .mobility import mobility
from .division import Division



'''
    init_df should have columns 'infective', 'inhabitant' and 'susceptible'.
    Setting iterative to False will reduce the running time but may result in
    a suboptimal clustering.
'''
def adaptive_mobility_regions(init_df,res=16000,mobility=None,iterative=True,**kwargs):
    if not mobility:
        from .mobility import mobility as mobmob
        mobility = mobmob
    if not 'infective' in init_df.columns:
        init_df['infective'] = init_df['infected_tested']+init_df['infected_nottested']
    G = mobility.G.copy()
    nx.set_node_attributes(G, init_df['infective'].to_dict(), 'infective')
    nx.set_node_attributes(G, init_df['susceptible'].to_dict(), 'susceptible')

    scorer = AdaptiveMobilityRegionsScorer(G=G,res=res,**kwargs)
    opt = Louvain(scorer,iterative=iterative)
    opt.optimize()
    print("Obtained score",opt.scorer.value())
    return opt.scorer.candidate()

'''
    Returns a division with two regions: One being division[isolated_area] and
    one obtained from merging all other regions.
'''
def isolation_division(division,isolated_area):
    return Division({
        i: 0 if division[i]==division[isolated_area] else 1
        for i in division.keys()
    })

'''
    Returns the mobility regions coresponding to a given resolution parameter res.
        mobility_regions(res=res)
    Set iterative=False to accelerate.
'''
def mobility_regions(res=2,mobility=None,iterative=True,**kwargs):
    if not mobility:
        from .mobility import mobility as mobmob
        mobility = mobmob
    scorer=ModularityScorer(mobility.G,res=res,**kwargs)
    opt = Louvain(scorer,iterative=iterative)
    opt.optimize()
    return opt.scorer.candidate()
class LouvainScorer:
    @staticmethod
    def aggregate_graph(G,C,agg_attrs=[],agg_edge_attrs=[]):
        newG = nx.DiGraph()
        newG.add_nodes_from(C.regions.keys())
        attrs_dict = {
            node: {
                **{
                    attr: sum([
                        G.nodes[subnode][attr]
                        for subnode in subnodes
                    ])
                    for attr in agg_attrs
                }, **{
                    attr: sum([
                        G[i][j][attr]
                        for i,j in G.subgraph(subnodes).edges
                    ]) + sum([
                        G.nodes[subnode][attr]
                        for subnode in subnodes
                        if attr in G.nodes[subnode]
                    ])
                    for attr in agg_edge_attrs
                }
            }
            for node,subnodes in C.regions.items()
        }
        nx.set_node_attributes(newG,attrs_dict)
        newEdges={}
        for i,j in G.edges:
            if not (C[i],C[j]) in newEdges:
                newEdges[C[i],C[j]] = 0
            newEdges[C[i],C[j]] += G[i][j]['weight']
        newG.add_weighted_edges_from([
            e+(w,)
            for e,w in newEdges.items()
        ])
        return newG

    def __init__(self, G, C0=None,agg_attrs=[],agg_edge_attrs=[]):
        self.G = G
        self.G_start = self.G.copy()
        if C0 is None:
            C0 = Division(dict(zip(G.nodes,G.nodes)))
        self.HC0 = C0
        self.HC = self.HC0.copy()
        self.performed_actions = []
        self.agg_attrs=agg_attrs
        self.agg_edge_attrs = agg_edge_attrs

        # Total weight of the graph
        self.in_degree = lambda x: (sum(d for _,d in self.G.in_degree(x, weight='weight')) if hasattr(x, '__iter__') else self.G.in_degree(x,weight='weight'))
        self.out_degree = lambda x: (sum(d for _,d in self.G.out_degree(x, weight='weight')) if hasattr(x, '__iter__') else self.G.out_degree(x,weight='weight'))
        self.total = self.in_degree(self.G.nodes)

        self.current_value = self.recalculate_value()

    # Returns the (current) candidate as a flat clustering
    def candidate(self):
        return self.HC.flatClustering()

    # Only called at initialization. Leaving this at 0 should not influence the
    # optimization while it will result in wrong values of current_value and value()
    def recalculate_value(self):
        return 0

    def value(self):
        return self.current_value

    def relabel_value(self,newlabel,node):
        pass
    def action_value(self,action):
        # Only do something for RelabelNode
        if type(action)==RelabelNode:
            return self.relabel_value(action.newlabel, action.node)
        return self.current_value

    # This is called right before the Aggregate action is performed. Can be
    # implemented to keep track of variables
    def will_aggregate(self):
        pass
    def will_recluster(self):
        pass
    def will_relabel(self,newlabel,node):
        pass
    def will_flatten(self):
        pass

    # The nodes who need to be revisited by MoveNodes after node has been moved to newlabel
    def revisit_nodes_after_relabel(self,relabel_action):
        node = relabel_action.node
        label = relabel_action.newlabel
        return set(self.G.predecessors(node)).union(self.G.successors(node))-self.HC.regions[label]

    # Updates clustering by performing the action
    def perform_action(self,action):
        if type(action) == Aggregate:
            # Aggregate the graph
            self.G = self.aggregate_graph(self.G,self.HC,agg_attrs=self.agg_attrs,agg_edge_attrs=self.agg_edge_attrs)
            self.will_aggregate()
        elif type(action) == Flatten:
            # Use the flat version of G
            self.G = self.G_start
            self.will_flatten()
        elif type(action) == RelabelNode:
            # Use the flat version of G
            self.current_value = self.relabel_value(action.newlabel,action.node)
            self.will_relabel(action.newlabel,action.node)
        elif type(action) != RelabelNode:
            print(action,"of type",type(action),"is not allowed for a",type(self))
        self.HC = action.perform(self.HC)
        self.performed_actions.append(action)

class ProductFormScorer(LouvainScorer):
    def __init__(self, G, node_in_weights, node_out_weights, C0=None):
        newG = G.copy()
        nx.set_node_attributes(newG,node_in_weights,'in_weight')
        nx.set_node_attributes(newG,node_out_weights,'out_weight')
        agg_attrs = ['in_weight','out_weight']
        super().__init__(newG,C0=C0,agg_attrs=agg_attrs)

    def initialize_aggregates(self):
        self.community_aggregates = {
            l: {
                'weight': self.weight(c,c,symmetric=False),
                **{
                    attr: sum(
                        self.G.nodes[i][attr]
                        for i in c
                    )
                    for attr in ['in_weight','out_weight']
                }
            }
            for l,c in self.HC.regions.items()
        }
    def community_contribution(self,weight,in_weight,out_weight):
        return weight - in_weight * out_weight
    def recalculate_value(self):
        # Recompute the aggregates
        self.initialize_aggregates()

        return sum(
            self.community_contribution(**aggs)
            for aggs in self.community_aggregates.values()
        )

    def relabel_value(self,newlabel,node,verbose=False):
        oldlabel = self.HC[node]
        newlabel_contribution_before = self.community_contribution(
            **self.community_aggregates[newlabel]
        ) if newlabel in self.community_aggregates else 0

        # Subtract previous contributions of newlabel and oldlabel and add updated ones
        return self.current_value - newlabel_contribution_before - self.community_contribution(
            **self.community_aggregates[oldlabel]
        ) + self.community_contribution(
            **self.update_aggregates_before_leave(oldlabel,node)
        ) + self.community_contribution(
            **self.update_aggregates_before_enter(newlabel,node)
        )

    # If withoutself and HC[i]==label, then we compute the penalty between {i} and label-{i}
    def node_community_penalty(self,i,label,withoutself=False):
        selfpenalty = 0
        if withoutself and self.HC[i]==label:
            selfpenalty = -2*self.G.nodes[i]['in_weight']*self.G.nodes[i]['out_weight']
        return selfpenalty + (
            self.community_aggregates[label]['in_weight'] * self.G.nodes[i]['out_weight']
            + self.community_aggregates[label]['out_weight'] * self.G.nodes[i]['in_weight']
        )

    def node_community_contribution(self,i,label):
        weight = sum([
            self.G[i][j]['weight'] for j in self.G.successors(i)
        ]) + sum([
            self.G[j][i]['weight'] for j in self.G.predecessors(i)
        ])
        return weight - self.node_community_penalty(i,label)

    def pair_contribution(self,i,j):
        return (
            self.G[i][j]['weight'] if j in self.G[i] else 0
        ) + (
            self.G[j][i]['weight'] if i in self.G[j] else 0
        ) - (
            self.G.nodes[i]['in_weight'] * self.G.nodes[j]['out_weight']
            + self.G.nodes[i]['out_weight'] * self.G.nodes[j]['in_weight']
        )

    def will_relabel(self,newlabel,node):
        self.community_aggregates[newlabel] = self.update_aggregates_before_enter(newlabel,node)
        oldlabel = self.HC[node]
        if len(self.HC.regions[oldlabel])==1:
            del self.community_aggregates[oldlabel]
        else:
            self.community_aggregates[oldlabel] = self.update_aggregates_before_leave(oldlabel,node)

    # Assuming node is not yet in HC[newlabel]
    def update_aggregates_before_enter(self,newlabel,node):
        if not newlabel in self.HC.regions.keys():
            #print(node,'wants to move to the new label',newlabel)
            # Then it's a move to a new label
            return {
                'weight': self.weight({node},{node},symmetric=False),
                **{
                    attr: self.G.nodes[node][attr]
                    for attr in self.agg_attrs
                }
            }
        if node in self.HC.regions[newlabel]:
            print("AggregationScorer.update_aggregates_before_enter error: node",node,"is in community",newlabel)
        aggs_before = self.community_aggregates[newlabel]
        newc = self.HC.regions[newlabel]
        return {
            'weight': aggs_before['weight'] + self.weight({node},newc,symmetric=True) + self.weight({node},{node},symmetric=False),
            **{
                attr: aggs_before[attr] + self.G.nodes[node][attr]
                for attr in self.agg_attrs
            }
        }

    # Assuming node in HC.regions[oldlabel]
    def update_aggregates_before_leave(self,oldlabel,node):
        if not node in self.HC.regions[oldlabel]:
            print("AggregationScorer.update_aggregates_before_leave error: node",node,"is not in community",oldlabel)
        aggs_before = self.community_aggregates[oldlabel]
        oldc = self.HC.regions[oldlabel]-{node}

        return {
            'weight': aggs_before['weight'] - self.weight({node},oldc,symmetric=True) - self.weight({node},{node},symmetric=False),
            **{
                attr: aggs_before[attr] - self.G.nodes[node][attr]
                for attr in self.agg_attrs
            }
        }

    def weight(self,c_start,c_end,symmetric=False,weight_attr='weight'):
        # Note that even if weight_attr!='weight', then still the constructor
        # of LouvainScorer gives self.G edge-weights with label 'weight'
        weight = sum(
            self.G[i][j][weight_attr]
            for i in c_start
            for j in set(c_end).intersection(self.G[i])
        )
        if symmetric:
            weight += sum(
                self.G[i][j][weight_attr]
                for i in c_end
                for j in set(c_start).intersection(self.G[i])
            )
        return weight

class ModularityScorer(ProductFormScorer):
    def __init__(self, G, C0=None, res=1):
        in_weight,out_weight = ModularityScorer.in_out_weights(G,res=res)
        super().__init__(G, in_weight, out_weight, C0=C0)

    def in_out_weights(G,res=1):
        in_degree = G.degree
        out_degree = G.degree
        if G.is_directed():
            in_degree = G.in_degree
            out_degree = G.out_degree
        total = sum([d for _,d in in_degree(G.nodes, weight='weight')])
        scaling = math.sqrt(res / total)
        in_weight = {
            i: scaling * in_degree(i, weight='weight')
            for i in G.nodes
        }
        out_weight = {
            i: scaling * out_degree(i, weight='weight')
            for i in G.nodes
        }
        return (in_weight,out_weight)

    def value(self,normalized=False):
        if normalized:
            return self.current_value/self.total
        return self.current_value

class AdaptiveMobilityRegionsScorer(ProductFormScorer):
    def __init__(self, G,
                 C0=None, res=1):
        self.total_inhabitants = sum(
            G.nodes[i]['inhabitant']
            for i in G.nodes
        )
        self.res=res
        in_weight,out_weight = AdaptiveMobilityRegionsScorer.in_out_weights(G,self.total_inhabitants,res=res)
        super().__init__(G, in_weight, out_weight, C0=C0)

    @staticmethod
    def in_out_weights(G,total_inhabitants,res=1):
        scaling = math.sqrt(res / total_inhabitants)
        in_weight = {
            i: scaling * G.nodes[i]['infective']
            for i in G.nodes
        }
        out_weight = {
            i: scaling * G.nodes[i]['susceptible']
            for i in G.nodes
        }
        return (in_weight,out_weight)

    def value(self):
        return self.current_value

class ObjectiveScorer(LouvainScorer):
    def __init__(self,init_df,horizon=30,tradeoff=1500,init_hc=None,mobility=mobility):
        from .mobility_seir import MobilitySEIR
        if init_hc is None:
            init_hc = Division(dict(zip(init_df.index,init_df.index)))
        if not isinstance(init_hc,Division):
            init_hc = Division(init_hc)
        self.tradeoff=tradeoff
        self.horizon=horizon
        self.seir = MobilitySEIR(init_df, horizon, division=init_hc.flatClustering())
        self.block2areas = {
            l: init_hc.flatitems({l})
            for l in init_hc.keys()
        }
        G = self.aggregate_graph(mobility.G,Division.FromRegions(self.block2areas))
        super().__init__(G,C0=Division(init_hc))

    def flattenblocks(self,blocks):
        return set().union(*(self.block2areas[block] for block in blocks))

    def relabel_value(self,newlabel,node):
        oldlabel = self.HC[node]
        oldlabel_oldscore = self.region2score[oldlabel]
        flatnodes = self.flattenblocks(self.HC.flatitems({node}))
        # find the areas corresponding to the new destination region
        r_dest = flatnodes
        if newlabel in self.seir.division.regions:
            r_dest |= self.seir.division.regions[newlabel]
        newlabel_oldscore = self.region2score[newlabel] if newlabel in self.region2score else 0
        newlabel_newscore = self.seir.objective(self.tradeoff,region=r_dest)
        # Lazy-load oldlabel_newscore
        if not oldlabel in self.region2block2leave_score:
            self.region2block2leave_score[oldlabel] = {}
        if not node in self.region2block2leave_score[oldlabel]:
            r_orig = self.seir.division.regions[oldlabel] - flatnodes
            self.region2block2leave_score[oldlabel][node] = self.seir.objective(self.tradeoff,region=r_orig)
        oldlabel_newscore = self.region2block2leave_score[oldlabel][node]
        return self.current_value + (
            newlabel_newscore - newlabel_oldscore
        ) + (
            oldlabel_newscore - oldlabel_oldscore
        )

    def recalculate_value(self):
        # Clean lazy-loader
        self.region2block2leave_score = {}
        # Simulate and store results
        self.seir.simulate_all()
        self.region2score = {
            r: self.seir.objective(self.tradeoff,region=r)
            for r in self.seir.division.labels()
        }
        self.current_value = sum(self.region2score.values())
        return self.current_value

    def revisit_nodes_after_relabel(self,relabel_action):
        return self.HC.nodes()

    def will_flatten(self):
        # Clean lazyloader
        self.region2block2leave_score = {}

    def will_relabel(self,newlabel,node):
        # update region2score. Note that we do this super inefficiently by repeating the simulation.
        oldlabel = self.HC[node]

        # Clean lazyloader
        if newlabel in self.region2block2leave_score:
            del self.region2block2leave_score[newlabel]
        if oldlabel in self.region2block2leave_score:
            del self.region2block2leave_score[oldlabel]

        # Move all nodes
        nodes = self.flattenblocks(self.HC.flatitems({node}))
        for i in nodes:
            self.seir.division[i] = newlabel
        self.seir.region2states[newlabel] = self.seir.simulate_region(newlabel)
        self.region2score[newlabel] = self.seir.objective(self.tradeoff,region=newlabel)
        # Check whether the old region still exists
        if oldlabel in self.seir.division.regions:
            self.seir.region2states[oldlabel] = self.seir.simulate_region(oldlabel)
            self.region2score[oldlabel] = self.seir.objective(self.tradeoff,region=oldlabel)
        else:
            del self.seir.region2states[oldlabel], self.region2score[oldlabel]
'''
    Definition of all the different actions/moves that the optimization algorithm makes.
'''
class Action:
    # Performs the action on the (hierarchical) clustering HC and return the resulting clustering
    def perform(self,C):
        pass
    # For reproducability purposes
    def __hash__(self):
        return str(self).__hash__()
    # For reproducability purposes
    def __eq__(self,obj):
        return type(self)==type(obj) and hash(self)==hash(obj)

# A non-action (doesn't change the clustering)
class SkipAction(Action):
    def __init__(self,description=""):
        self.description = str(description)
    def __str__(self):
        if len(self.description)>0:
            return "skip {}".format(self.description)
        return "skip"
    # For reproducability purposes
    def __hash__(self):
        return (type(self),self.description).__hash__()

# A non-action (doesn't change the clustering)
class Done(Action):
    def __str__(self):
       return "optimization finished"

# A non-action (doesn't change the clustering)
class Restart(Action):
    def __str__(self):
       return "Start a new iteration"

class RelabelNode(Action):
    def __init__(self, node, newlabel):
        self.node = node
        self.newlabel = newlabel
    def perform(self,HC):
        HC[self.node] = self.newlabel
        return HC
    def __str__(self):
        return "relabeling node {} to label {}".format(self.node,self.newlabel)
    def __hash__(self):
        return (type(self),self.node,self.newlabel).__hash__()

# Aggregate the clusters into new nodes by going up a level in the hierarchy
class Aggregate(Action):
    print_str = "aggregating clusters"
    def perform(self,HC):
        self.print_str = "aggregating {} clusters".format(len(HC.regions))
        return HC.nextlevel()
    def __str__(self):
        return self.print_str
# Flatten the Hierarchical clustering into a regular clustering
class Flatten(Action):
    def perform(self,HC):
        return HC.flatClustering()
    def __str__(self):
        return "flattening clusters"

'''
    Definition of the optimization algorithm.
'''
class Louvain:
    '''
    I like to view Louvain as a state diagram with the states represented by
    the enum LouvainState. The actions cause the following transitions:
        RelabelNode: MOVENODES_UNIMPROVED/MOVENODES_IMPROVED -> MOVENODES_IMPROVED
        Aggregate:   PASS_HAS_IMPROVED -> MOVENODES_UNIMPROVED
        Flatten: PASS_WAS_STABLE -> ITERATION_ENDED
        SkipAction:
            if the queue is empty:
                MOVENODES_UNIMPROVED -> PASS_WAS_STABLE
                MOVENODES_IMPROVED -> PASS_HAS_IMPROVED
            else we remain in state MOVENODES_UNIMPROVED/MOVENODES_IMPROVED
    Louvain can be applied iteratively by transitioning from ITERATION_ENDED to
    MOVENODES_UNIMPROVED (whenever the iteration performed RelabelNode at least once).
    After each occurrence of RelabelNode or SkipAction, current_node is taken from
    the queue.
    '''
    class LouvainState(Enum):
        MOVENODES_UNIMPROVED=0
        MOVENODES_IMPROVED=1
        PASS_WAS_STABLE=2
        PASS_HAS_IMPROVED=3
        ITERATION_ENDED=4

    def __init__(self,scorer,iterative=False,random=False,seed=None):
        self.scorer = scorer
        # Keep track of the action_index at the start of the iteration so that
        # At the end of the iteration, we can see whether we can terminate.
        self.action_index_start = len(scorer.performed_actions)

        self.iterative = iterative
        if random:
            self.random = np.random.RandomState(seed=seed)
        else:
            self.random = False

        self.state = Louvain.LouvainState.MOVENODES_UNIMPROVED
        self.start_pass()

    # Iterator that returns the next action at each point.
    # NOTE THAT IT DOES NOT PERFORM THE ACTION so that iterating over this
    # without performing results in an infinite loop
    def action_sequence(self):
        while True:
            action = self.select_action(self.candidate_actions())
            yield (action)
            if isinstance(action, Done):
                break  # Optimization finished

    def optimize(self,verbose=False):
        for action in self.action_sequence():
            if verbose:
                print(action)
            if not type(action) in (Done,SkipAction,Restart):
                self.scorer.perform_action(action)
            self.update_state(action)

    def next(self):
        if not hasattr(self,'action_iterator'):
            self.action_iterator = self.action_sequence()

        action = next(self.action_iterator)
        if not type(action) in (Done,SkipAction):
            self.scorer.perform_action(action)
        self.update_state(action)
        print(action)

    def start_pass(self):
        nodes = self.scorer.HC.keys()
        # Randomize order
        if self.random:
            nodes = list(nodes)
            self.random.shuffle(nodes)
        self.queue = deque(nodes)
        self.stable_nodes = []
        self.current_node = self.queue.popleft()

    def iteration_ended(self):
        if not self.iterative:
            return # We are done, the algorithm has terminated
        # Check if RelabelNode was performed since the first aggregation of this iteration
        new_actions = self.scorer.performed_actions[self.action_index_start:]
        aggregated = False
        for action in new_actions:
            if type(action)==Aggregate:
                aggregated = True
            if aggregated and type(action)==RelabelNode:
                self.action_index_start = len(self.scorer.performed_actions)
                self.state = Louvain.LouvainState.MOVENODES_UNIMPROVED
                self.start_pass()
                return
        # Else, The algorithm has terminated

    def candidate_actions(self):
        if self.state == Louvain.LouvainState.ITERATION_ENDED:
            # Terminate
            return set()
        if self.state == Louvain.LouvainState.PASS_WAS_STABLE:
            return {Flatten()}
        if self.state == Louvain.LouvainState.PASS_HAS_IMPROVED:
            return {Aggregate()}

        # Else we return the RelabelNode candidates
        labels = {self.scorer.HC.newlabel()}
        # Don't allow to move to new community if it is in a singleton
        if len(self.scorer.HC.regions[self.scorer.HC[self.current_node]])==1:
            labels=set()
        # Check if the graph is aggregated
        if self.current_node in self.scorer.G.nodes:
            # If so, return the labels of the neighbors of i
            labels = labels.union({
                self.scorer.HC[j]
                for j in it.chain(
                        self.scorer.G.neighbors(self.current_node),
                        self.scorer.G.predecessors(self.current_node))
            } - {self.scorer.HC[self.current_node]})
        else:
            # If not, we simply return all other labels
            labels = labels.union(self.scorer.HC.labels() - {self.scorer.HC[self.current_node]})
        # We add the SkipAction, so that there is always at least one action returned.
        # If only one action is returned (a forced action), it will be the SkipAction.
        return {
            RelabelNode(self.current_node, label)
            for label in labels
        }.union({SkipAction(self.current_node)})

    def select_action(self,candidates):
        candidates = {
            c: self.scorer.action_value(c)
            for c in candidates
        }
        if len(candidates) == 0:
            # Iteration ends if there are no candidates left
            return Done()
        # If candidate is a dict, we assume its values are the scores.
        if isinstance(candidates,dict):
            best_action = max(candidates,key=candidates.get)
        # If there is only one candidate, then we consider it a 'forced move'
        # Otherwise, only perform an action if it increases the score
        best_action = max(candidates, key=self.scorer.action_value)
        if len(candidates) == 1 or self.scorer.action_value(best_action) > self.scorer.value():
            return best_action
        else:
            return SkipAction("tie")

    def update_state(self,action):
        if type(action)==RelabelNode:
            self.state = Louvain.LouvainState.MOVENODES_IMPROVED
            # We add all nodes that were previously stable, since this change
            # may make them unstable again.
            # The speed could be significantly improved by efficiently identifying
            # which nodes definitely remained stable so that they do not need to
            # requeued.
            self.queue.extend(self.stable_nodes)
            self.stable_nodes = [self.current_node]
            # len(queue)==len(HC.keys())-len(stable_nodes)==len(HC.keys())-1
            # while an occurrance of RelabelNode implies len(HC.keys())>1.
            # Therefore the queue must be nonempty at this point
            self.current_node = self.queue.popleft()
        if type(action)==SkipAction:
            self.stable_nodes.append(self.current_node)
            if len(self.queue)>0:
                self.current_node = self.queue.popleft()
            elif self.state==Louvain.LouvainState.MOVENODES_IMPROVED:
                self.state = Louvain.LouvainState.PASS_HAS_IMPROVED
            else: # self.state==Louvain.LouvainState.MOVENODES_UNIMPROVED
                self.state = Louvain.LouvainState.PASS_WAS_STABLE

        if type(action)==Aggregate:
            self.state = Louvain.LouvainState.MOVENODES_UNIMPROVED
            self.start_pass()
        if type(action)==Flatten:
            self.state = Louvain.LouvainState.ITERATION_ENDED
            self.iteration_ended()