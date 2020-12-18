class CoronaConstants:
    # Disease constants
    latent_period = 4               # nu 
    infectious_period = 5           # omega  
    fraction_tested = 1/15          # alpha
    contacts_average = 13.85        # c
    fraction_local_contacts = 0.5   # rho
    basic_reproduction_number = 1.25   
    transmission_prob = 1.25*0.01390958192240339      # epsilon
    
    # Mobility constants
    population_nl = 17242149          # population in the Netherlands
    average_total_mobility = 10069260 # sum of all mezuro data divided by number of days
    
    contacts_local = fraction_local_contacts*contacts_average   # rho*c
    contacts_per_visit = (1-fraction_local_contacts)*contacts_average*population_nl/(2*average_total_mobility) # c_m : see supplementary material

    # dictionary containing the changed attributes, set in __init__ so that CoronaConstants.changed
    # is always {}
    changed = {}

    # To change constants, simply pass keyword-arguments. E.g. CoronaConstants(latent_period=5)
    def __init__(self,**kwargs):
        self.changed = kwargs
        for key,value in kwargs.items():
            setattr(self, key, value)

        # Recompute dependent constants
        if 'transmission_prob' not in kwargs:
            self.transmission_prob = self.basic_reproduction_number * 0.01390958192240339  # epsilon
        if 'contacts_local' not in kwargs:
            self.contacts_local = self.fraction_local_contacts * self.contacts_average  # rho*c
        if 'contacts_per_visit' not in kwargs:
            self.contacts_per_visit = (1 - self.fraction_local_contacts) * self.contacts_average * self.population_nl / (
                2 * self.average_total_mobility
            )  # c_m : see supplementary material
