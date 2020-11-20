from .mobility_seir import MobilitySEIR
from .rivm_loader import rivm
from .constants import CoronaConstants
from .mezuro_preprocessing import gemeente_shapes
init_df = rivm.SEI2R2_init('1030')
for r,r_desc in zip([0.9,1.1],['09','11']):
    for contacts,c_desc in zip([{},{'contacts_average': 13.85*0.75, 'fraction_local_contacts': 2/3}],['normal','mobility_halved']):
        seir = MobilitySEIR(init_df,14,constants=CoronaConstants(basic_reproduction_number=r,**contacts))
        seir.simulate_all()
        infections = seir.daily_reported_infections()
        infections['gem_id'] = [
            gemeente_shapes.loc[a,'gem_id']
            for a in infections.index
        ]
        infections.to_csv('predicted_cases_r_{}_contacts_{}.csv'.format(r_desc,c_desc))

