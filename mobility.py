'''
    Usage
        from .mobility import Mobility
        mobility = Mobility.AverageGemeenteMobility()
        mobility[origin,destination]
    gives the number of visits from people living in origin to destination.
    To get the destinations corresponding to an origin
        mobility.destinations(origin)
    To get the origins corresponding to a destination
        mobility.origins(destination)
'''
class Mobility(dict):
    # For a networkx.DiGraph G that has edge attributes named 'weight'
    def __init__(self,G):
        self.G = G

        # For modularity
        self.total_weight = sum(self.G[i][j]['weight'] for i,j in self.G.edges)

    # Assume pair is a tuple (origin,destination).
    def __getitem__(self,pair):
        if pair in self.G.edges:
            origin,destination=pair
            return self.G[origin][destination]['weight']
        # Consider mobility zero else
        return 0

    def destinations(self,origin):
        return self.G.successors(origin)
    def origins(self,destination):
        return self.G.predecessors(destination)
    def neighbors(self,i):
        return set(self.destinations(i)).union(self.origins(i))

    def __iter__(self):
        return iter(self.G.edges)
    def keys(self):
        return self.G.edges

    # Returns modularity contribution of pair i,j
    def modularity(self,i,j,res=1):
        return self[i,j] - res * (
            self.G.out_degree(i,weight='weight')*self.G.in_degree(j,weight='weight')
        ) / self.total_weight

    @staticmethod
    def AverageGemeenteMobility():
        from .mezuro_preprocessing import gemeente2gemeente,gemeente_shapes
        from .networks import avg_graph
        return Mobility(avg_graph(
            gemeente_shapes,
            gemeente2gemeente,
            col='totaal_aantal_bezoekers'
        ))

    def subset_mobility(self,r):
        return sum([
            self[i,j]
            for i in r
            for j in r.intersection(self.destinations(i))
        ])

    def mobility_between(self,r1,r2):
        return sum([self[i,j]+self[j,i] for i in r1 for j in r2])

    def intra_region_mobility(self,division):
        from .division import Division
        division = Division(division)
        return sum(
            self.subset_mobility(c)
            for c in division.regions.values()
        )

    def percentage_mobility(self,division):
        return 100 * self.intra_region_mobility(division) / self.total_weight

mobility = Mobility.AverageGemeenteMobility()
