import os

import pandas as pd
from .constants import CoronaConstants
from .division import Division


class MobilitySEIR:
    compartments = [
        'susceptible',
        'exposed',
        'infected_tested',
        'infected_nottested',
        'removed_tested',
        'removed_nottested'
    ]

    '''
        init_df: DataFrame (indexed by areas) to initialize the state of the model.
            Can have columns for each of the compartments but should otherwise at
            least have a column 'inhabitant'.
    '''

    def __init__(self, init_df, horizon=0,
                 mobility=None,
                 division=None,
                 constants=CoronaConstants):
        self.init_df = init_df[[c for c in init_df.columns if c in MobilitySEIR.compartments]].copy()
        self.horizon = horizon
        if not set(self.compartments).issubset(init_df.columns):
            # assume 'inhabitant' in init_df.columns, place the whole population
            # into the first compartment.
            self.init_df[self.compartments[0]] = init_df['inhabitant'].copy()
            for i in range(1, len(self.compartments)):
                self.init_df[self.compartments[i]] = 0

        self.region2states = dict()

        # We represent a division as a dict indexed by the areas. Areas i,j are
        # in the same region if division[i]==division[j].
        # The division representing no restrictions:
        self.division_no_restrictions = {
            area: 0
            for area in init_df.index
        }
        self.division = division if not division is None else self.division_no_restrictions
        self.division = Division(self.division)

        if not mobility:
            from .mobility import mobility as mob
            mobility = mob
        self.mobility = mobility

        # Copy all the constants
        self.constants = constants
        self.transmission_prob = constants.transmission_prob
        self.contacts_per_visit = constants.contacts_per_visit
        self.contacts_local = constants.contacts_local
        self.infectious_period = constants.infectious_period
        self.latent_period = constants.latent_period
        self.fraction_tested = constants.fraction_tested

    # Region should either be a label in self.division or a subset of the areas.
    def simulate_region(self, region):
        import collections
        if isinstance(region, collections.Hashable) and region in self.division.regions:
            region = self.division.regions[region]

        state = self.init_df.loc[region, :].to_dict(orient='index')
        states = {
            0: state
        }
        contacts = {
            area: region.intersection(self.mobility.origins(area))
                        .union(region.intersection(self.mobility.destinations(area)))
            for area in region
        }
        # Returns a dictionary from area to a tuple of the compartments
        unpack = lambda state, areas: [
            (area, tuple(state[area][c] for c in self.compartments))
            for area in areas
        ]

        for t in range(1, self.horizon + 1):
            # deep copy
            new_state = {
                i: s.copy()
                for i, s in state.items()
            }
            for area, (S, E, It, Int, Rt, Rnt) in unpack(state, region):
                N = S + E + It + Int + Rt + Rnt
                # Mobility exposures
                exposures = self.transmission_prob * self.contacts_per_visit * S * sum([
                    Int_mob * (self.mobility[area, mob_area] + self.mobility[mob_area, area])
                    / (S_mob + E_mob + Int_mob + Rt_mob + Rnt_mob)
                    for mob_area, (S_mob, E_mob, It_mob, Int_mob, Rt_mob, Rnt_mob) in unpack(state, contacts[area])
                ]) / (S + E + Int + Rt + Rnt)
                # Local exposures
                exposures += self.transmission_prob * self.contacts_local * S * (
                        It + Int
                ) / N
                # Exposed to infected_tested
                infections_tested = E * self.fraction_tested / self.latent_period
                # Exposed to infected_nottested
                infections_nottested = E * (1 - self.fraction_tested) / self.latent_period
                # infected_tested to removed_tested
                removals_tested = It / self.infectious_period
                # infected_nottested to removed_nottested
                removals_nottested = Int / self.infectious_period
                # Update
                new_state[area]['susceptible'] -= exposures
                new_state[area]['exposed'] += exposures - infections_tested - infections_nottested
                new_state[area]['infected_tested'] += infections_tested - removals_tested
                new_state[area]['infected_nottested'] += infections_nottested - removals_nottested
                new_state[area]['removed_tested'] += removals_tested
                new_state[area]['removed_nottested'] += removals_nottested
            state = new_state
            states[t] = new_state
        return states

    def state_at_horizon(self):
        df = pd.DataFrame.from_dict({
            area: self.region2states[self.division[area]][self.horizon][area]
            for area in self.init_df.index
        },orient='index')
        df.index.name = 'name'
        return df

    def daily_reported_infections(self):
        cumulative = lambda x: x['infected_tested']+x['removed_tested']
        cumulatives = {
            area: [
                cumulative(self.region2states[self.division[area]][t][area])
                for t in range(self.horizon+1)
            ]
            for area in self.init_df.index
        }
        df = pd.DataFrame(index=self.init_df.index)
        df.index.name = 'name'
        for t in range(1,self.horizon+1):
            df['Infections on day {:02d}'.format(t)] = [
                cumulatives[area][t]-cumulatives[area][t-1]
                for area in df.index
            ]
        return df

    def simulate_all(self):
        self.region2states = {
            r: self.simulate_region(r)
            for r in self.division.labels()
        }

    def cumulative_infections(self, t=None, start_t=None, region=None):
        if region is None:
            return sum([
                self.cumulative_infections(t, start_t=start_t, region=r)
                for r in self.division.labels()
            ])
        if t is None:
            t = self.horizon
        if not isinstance(region, dict):
            # assume region is a key of region2states
            region_states = self.region2states[region]
        else:
            region_states = region
        susceptible_end = sum([s['susceptible'] for s in region_states[t].values()])
        if start_t == None:
            start_t = 0
        susceptible_start = sum([s['susceptible'] for s in region_states[start_t].values()])
        return susceptible_start - susceptible_end

    def objective(self, tradeoff, horizon=None, region=None, return_performance=False):
        import collections
        if horizon is None:
            horizon = self.horizon
        objective = 0
        movements = 0
        infections = 0
        if region is None:
            for r in self.division.labels():
                o, m, i = self.objective(tradeoff, horizon, region=r, return_performance=True)
                objective += o
                movements += m
                infections += i
        elif isinstance(region, collections.Hashable) and region in self.division.regions:
            movements = horizon * self.mobility.subset_mobility(self.division.regions[region])
            infections = self.cumulative_infections(region=region)
        else:
            # Assume region is a set of areas, perform the simulation
            movements = horizon * self.mobility.subset_mobility(region)
            states = self.simulate_region(region)
            infections = self.cumulative_infections(region=states)
        objective = movements - tradeoff * infections
        if return_performance:
            return (objective, movements, infections)
        return objective

    def merge_increase(self, r1, r2, tradeoff):
        movements = self.horizon * self.mobility.mobility_between(
            self.division.regions[r1],
            self.division.regions[r2]
        )
        infections_old = sum(self.cumulative_infections(region=r) for r in [r1, r2])
        r_new = self.division.regions[r1].union(self.division.regions[r2])
        states_new = self.simulate_region(r_new)
        infections_new = self.cumulative_infections(region=states_new)
        return movements - tradeoff * (infections_new - infections_old)

    '''
        Measures the heterogeneity of the distribution of infections throughout
        the country by computing the Kullback-Leibler divergence between two
        distributions:
            * The distribution of the area that a randomly chosen inhabitant lives in.
            * The distribution of the area that a randomly chosen infective inhabitant lives in.
        (i.e., the KL-divergence between the population distribution and the
        population distribution conditioned on infectiveness)

        This also relates to mutual information: Let the r.v. A denote the area
        of a randomly chosen inhabitant and let I denote the indicator that this
        inhabitant is infected. Then this corresponds to the mutual information
        between I and A.

        Returns zero for a homogeneous distribution of infections and is maximized
        when all infections reside in the smallest area.
        Note that this quantity does not depend on the absolute amount of infections.
    '''
    def initial_divergence(self):
        import numpy as np
        state = self.init_df.to_dict(orient='index')
        # Compute population distribution infective population distribution.
        population = [sum(a.values()) for a in state.values()]
        population_dist = np.array(population) / sum(population)

        infective = [a['exposed'] + a['infected_tested'] + a['infected_nottested'] for a in state.values()]
        infective_dist = np.array(infective) / sum(infective)

        # We use the convention that 0*log(0)=0
        divergence = sum([
            p_inf * np.log(p_inf / p_pop)
            for p_inf, p_pop in zip(infective_dist, population_dist)
            if p_inf > 0
        ])
        # Numerical errors may cause the resulting value to be negative
        # (of the order e-16), therefore we take the maximum with zero.
        return abs(divergence)

    def initial_concentration(self):
        import numpy as np
        return 1 - np.exp(-self.initial_divergence())

    def basic_reproduction_number(self):
        from numpy import linalg as LA
        K = next_generation_matrix(
            population=sum(self.init_df[c] for c in self.compartments),
            const=self.constants,
            mob=self.mobility,
            areas=self.init_df.index
        )
        w, v = LA.eig(K.transpose())
        basic_reproduction_number = max(abs(w))
        return basic_reproduction_number


def load_initialization(init_name):
    return pd.read_csv(init_path(init_name), index_col='name')

def init_path(init_name):
    return os.path.join(
        os.path.dirname(__file__),
        'initializations',
        init_name
    )

def generate_initializations():
    from .rivm_loader import rivm
    rivm0310 = rivm.SEI2R2_init('0310')
    rivm0421 = rivm.SEI2R2_init('0421')
    concentrated = rivm.concentrated_init(1000,'Uden')
    evenlydistributed = rivm.evenlydistributed_init(1000)
    seir_conc,seir_evenly = (MobilitySEIR(df,horizon=10) for df in [concentrated,evenlydistributed])
    seir_conc.simulate_all(),seir_evenly.simulate_all()
    concentrated,evenlydistributed = (seir.state_at_horizon() for seir in [seir_conc,seir_evenly])
    concentrated.to_csv(init_path('concentrated.csv'))
    evenlydistributed.to_csv(init_path('evenlydistributed.csv'))
    rivm0310.to_csv(init_path("historical0310.csv"))
    rivm0421.to_csv(init_path("historical0421.csv"))

def next_generation_matrix(population=None, susceptible=None, mob=None, const=CoronaConstants, areas=None):
    import numpy as np
    import itertools as it
    if mob is None:
        from .mobility import mobility
        mob = mobility
    if population is None:
        from .mezuro_preprocessing import gemeente_shapes
        population = gemeente_shapes['inhabitant']
    if susceptible is None:
        susceptible = population
    if areas is None:
        areas = list(mob.G.nodes)
    # Diagonal entries in Next Generation Matrix are:
    local_basic_reproduction_number = (
            const.contacts_local * const.transmission_prob * const.infectious_period
    )
    # Off-diagonal entries depend on mobility times a constant factor
    off_diagonal_factor = (
            (1 - const.fraction_tested) * const.contacts_per_visit * const.transmission_prob * const.infectious_period
    )
    # Next Generation Matrix:
    K = np.zeros((len(areas),len(areas)))
    for (i,area),(j,other_area) in it.product(enumerate(areas),enumerate(areas)):
        if i == j:
            K[i, j] = local_basic_reproduction_number * susceptible[area] / population[area]
        else:
            K[i, j] = off_diagonal_factor * (
                (mob[area,other_area] + mob[other_area,area]) / population[area]
            ) * susceptible[other_area] / population[other_area]
    return K
