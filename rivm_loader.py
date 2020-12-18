import os
import pandas as pd

from .constants import CoronaConstants
from datetime import datetime, timedelta
from .mezuro_preprocessing import koppeltabel, gemeente_shapes

class RivmLoader:
    def __init__(self):
        fn = os.path.join(os.path.dirname(__file__), 'data/COVID-19_aantallen_gemeente_per_dag.csv')
        rivm = pd.read_csv(fn, delimiter=';')
        mmdd2area2reported = {}
        for _, (date, area, reported) in rivm[~rivm.Municipality_name.isnull()][
            ['Date_of_publication', 'Municipality_name', 'Total_reported']].iterrows():
            mmdd = ''.join(date.split('-')[1:])
            if not mmdd in mmdd2area2reported:
                mmdd2area2reported[mmdd] = {}
            mmdd2area2reported[mmdd][area] = reported
        self.mmdd2area2reported = {
            mmdd: pd.Series(area2reported)
            for mmdd, area2reported in mmdd2area2reported.items()
        }
        cumulative = 0
        mmdd2cumulatives = {}
        for mmdd, infections in self.mmdd2area2reported.items():
            cumulative = infections + cumulative
            mmdd2cumulatives[mmdd] = cumulative
        self.mmdd2cumulatives = mmdd2cumulatives

    def rivm_corona(self, mmdd):
        if not mmdd in self.mmdd2cumulatives:
            print(mmdd, 'not in rivm data')
            return 0
        data = self.mmdd2cumulatives[mmdd]
        return pd.Series({
            i: data[koppeltabel.loc[i, 'rivm']] * koppeltabel.loc[i, 'inhabitant_frac_new']
            for i in gemeente_shapes.index
        })

    def susceptible(self,
                    mmdd,
                    totals=gemeente_shapes['inhabitant'],
                    undetected_multiplier=1 / CoronaConstants.fraction_tested,
                    latent_period=CoronaConstants.latent_period):
        mmdd_end_exposure = RivmLoader.date2mmdd(RivmLoader.mmdd2date(mmdd) + timedelta(days=latent_period))
        return totals - undetected_multiplier * self.rivm_corona(mmdd_end_exposure)

    def exposed(self,
                mmdd,
                undetected_multiplier=1 / CoronaConstants.fraction_tested,
                latent_period=CoronaConstants.latent_period):
        date = RivmLoader.mmdd2date(mmdd)
        mmdd_end_exposure = RivmLoader.date2mmdd(date + timedelta(latent_period))
        return undetected_multiplier * (
                self.rivm_corona(mmdd_end_exposure) - self.rivm_corona(mmdd)
        )

    def infected_tested(self,
                        mmdd,
                        infectious_period=CoronaConstants.infectious_period):
        date = RivmLoader.mmdd2date(mmdd)
        mmdd_start_infections = RivmLoader.date2mmdd(date - timedelta(infectious_period))
        return self.rivm_corona(mmdd) - self.rivm_corona(mmdd_start_infections)

    def infected_nottested(self,
                           mmdd,
                           infectious_period=CoronaConstants.infectious_period,
                           undetected_multiplier=1 / CoronaConstants.fraction_tested):
        date = RivmLoader.mmdd2date(mmdd)
        mmdd_start_infections = RivmLoader.date2mmdd(date - timedelta(infectious_period))
        return (undetected_multiplier - 1) * (
                self.rivm_corona(mmdd) - self.rivm_corona(mmdd_start_infections)
        )

    def removed_tested(self,
                       mmdd,
                       infectious_period=CoronaConstants.infectious_period):
        date = RivmLoader.mmdd2date(mmdd)
        mmdd_start_infections = RivmLoader.date2mmdd(date - timedelta(infectious_period))
        return self.rivm_corona(mmdd_start_infections)

    def removed_nottested(self,
                          mmdd,
                          infectious_period=CoronaConstants.infectious_period,
                          undetected_multiplier=1 / CoronaConstants.fraction_tested):
        date = RivmLoader.mmdd2date(mmdd)
        mmdd_start_infections = RivmLoader.date2mmdd(date - timedelta(infectious_period))
        return (undetected_multiplier-1) * self.rivm_corona(mmdd_start_infections)

    def SEI2R2_init(self,
                    mmdd='0421',
                    infectious_period=CoronaConstants.infectious_period,
                    undetected_multiplier=1 / CoronaConstants.fraction_tested,
                    latent_period=CoronaConstants.latent_period,
                    return_integer=False
                    ):
        df = pd.DataFrame()
        df["susceptible"] = self.susceptible(mmdd=mmdd, latent_period=latent_period,
                                             undetected_multiplier=undetected_multiplier)
        df["exposed"] = self.exposed(mmdd,undetected_multiplier = undetected_multiplier,latent_period=latent_period)
        df["infected_tested"] = self.infected_tested(mmdd=mmdd, infectious_period=infectious_period)
        df["infected_nottested"] = self.infected_nottested(mmdd=mmdd, infectious_period=infectious_period,
                                                           undetected_multiplier=undetected_multiplier)
        df["removed_tested"] = self.removed_tested(mmdd=mmdd, infectious_period=infectious_period)
        df["removed_nottested"] = self.removed_nottested(mmdd=mmdd, infectious_period=infectious_period,
                                            undetected_multiplier=undetected_multiplier)
        # replace negative values in dataframe by 0, negatives occur due to corrections in RIVM data
        num = df._get_numeric_data()
        num[num < 0] = 0
        if return_integer:
            df = df.round(0)
        df["inhabitant"] = df["susceptible"] + df["exposed"] + df["infected_tested"] + df["infected_nottested"] + df["removed_tested"] + df["removed_nottested"]
        return df

    def concentrated_init(self, people_exposed, municipality):
        df = pd.DataFrame()
        df["susceptible"] = gemeente_shapes['inhabitant']
        df["exposed"]=0
        df["infected_tested"]=0
        df["infected_nottested"]=0
        df["removed_tested"]=0
        df["removed_nottested"]=0
        # Test if people exposed < inhabitants:
        if df.loc[municipality]['susceptible']>=people_exposed:
            df.loc[df.index == municipality, 'exposed'] += people_exposed
            df.loc[df.index == municipality, 'susceptible'] -= people_exposed
            return df
        else:
            print("Mistake: Exposed can not be larger than Inhabitants")
            return None

    def evenlydistributed_init(self, people_exposed, return_integer=False):
        import numpy as np
        total_population = np.sum(gemeente_shapes['inhabitant'])
        df = pd.DataFrame()
        df["susceptible"] = gemeente_shapes['inhabitant']-people_exposed/total_population*gemeente_shapes['inhabitant']
        df["exposed"]=people_exposed/total_population*gemeente_shapes['inhabitant']
        df["infected_tested"]=0
        df["infected_nottested"]=0
        df["removed_tested"]=0
        df["removed_nottested"]=0
        if return_integer:
            return df.round(0)
        return df

    @classmethod
    def mmdd2date(cls,mmdd,year=2020,hour=14):
        return datetime(day=int(mmdd[2:]),month=int(mmdd[:2]),year=year,hour=hour)
    @classmethod
    def date2mmdd(cls,date):
        return f'{date.month:02d}{date.day:02d}'

rivm = RivmLoader()