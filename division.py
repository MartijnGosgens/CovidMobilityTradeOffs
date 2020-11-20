import pandas as pd
import numpy as np

class Division(dict):
    def __init__(self, division_dict, previouslevel=None):
        super().__init__(division_dict)
        self.regions = {}
        for i, c in division_dict.items():
            if not c in self.regions:
                self.regions[c] = set()
            self.regions[c].add(i)
        if previouslevel:
            previouslevel = previouslevel.copy()
        self.previouslevel = previouslevel

    def FromRegions(regions):
        return Division({
            i: l
            for l,c in regions.items()
            for i in c
        })

    def copy(self):
        c = super().copy()
        previous = self.previouslevel
        if previous != None:
            previous = previous.copy()
        return Division(c,previouslevel=previous)

    def __setitem__(self, key, value):
        self.regions[self[key]].remove(key)
        if len(self.regions[self[key]]) == 0:
            del self.regions[self[key]]
        if not value in self.regions:
            self.regions[value] = set()
        self.regions[value].add(key)
        super().__setitem__(key, value)

    def __str__(self):
        return self.regions.__str__()

    def labels(self):
        return self.regions.keys()

    def newlabel(self):
        return 1+max([-1,*(l for l in self.labels() if isinstance(l,int))])


    def relabel_regions(self,random=False):
        newlabels = list(range(len(self.labels())))
        if random:
            np.random.shuffle(newlabels)
        mapping = dict(zip(self.labels(), newlabels))
        return Division({
            a: mapping[r]
            for a,r in self.items()
        })

    def to_csv(self,file_name,index_label='item',cluster_label='cluster'):
        import pandas as pd
        pd.Series(self).to_csv(file_name,index_label=index_label,header=[cluster_label])

    @staticmethod
    def FromCSV(path):
        return Division(pd.read_csv(path, index_col='item')['cluster'].to_dict())
    '''
        The following functions allow for the division to be viewed as a hierarchical object.
    '''
    def flatClustering(self):
        if self.previouslevel == None:
            return self
        return Division({
            i: self[c]
            for i, c in self.previouslevel.flatClustering().items()
        })

    def flatitems(self,topitems):
        if self.previouslevel == None:
            return topitems
        return set().union(*[
            self.previouslevel.flatitems(self.previouslevel.regions[i])
            for i in topitems
        ])

    # Add an additional level to the hierarchy
    def nextlevel(self):
        return Division(dict(zip(self.labels(),self.labels())),self)

