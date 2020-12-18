import pandas as pd
import geopandas as gdp
from collections import Counter
import os
from .division import Division

# Inside data/mezuro_data, the zip-files 'nederland' and 'Mezuro shapefiles 2019'
# should be unpacked. These will be ignored by the git.
def mezuro_path(relative_path):
    return os.path.join(
        os.path.dirname(__file__),
        'data/mezuro_data',
        relative_path
    )
def data_path(relative_path):
    return os.path.join(
        os.path.dirname(__file__),
        'data',
        relative_path
    )

'''
This file gathers the relevant data and stores them in the following variables:
    - GeoDataFrames, containing name, population and province:
        - provincie_shapes
        - gemeente_shapes
            (with a column 'veiligheidsregio')
        - mezuro_shapes

    - Mobility dataframes
        - provincie2provincie
        - gemeente2gemeente
        - mezuro2mezuro

    - koppeltabel: a mapping from 2019 gemeentes to 2020 gemeentes

    - mezuro_hierarchy: represents
    the hierarchical clustering of mezuro shapes into gemeente shapes into provinces.
'''

provincie_shapes=(
    gdp.read_file(mezuro_path('Mezuro shapefiles 2019/shapefiles/provincies_2019.shp'))
    .rename(columns={'prov_name':'name'})
    .set_index('name')
)

gemeente_shapes=(
    gdp.read_file(mezuro_path('Mezuro shapefiles 2019/shapefiles/gemeenten_2019.shp'))
    .rename(columns={'gem_name':'name'})
    .set_index('name')
)

provincie2provincie=(
    pd.read_csv(mezuro_path(
        'nederland/provincie niveau/per dag/Herkomsten bezoekers per dag - van provincie naar provincie - 01-03-2019 tm 14-03-2019 - ilionx.csv'
    ), delimiter=';')
    .rename(columns={'woon_provincie_naam': 'woon','bezoek_provincie_naam': 'bezoek'})
    .set_index(['woon','bezoek','datum'])
)
gemeente2gemeente=(
    pd.read_csv(mezuro_path(
        'nederland/gemeente niveau/per dag/Herkomsten bezoekers per dag - van gemeente naar gemeente - 01-03-2019 tm 14-03-2019 - ilionx.csv'
    ), delimiter=';')
    .rename(columns={'woon_gemeente_naam': 'woon','bezoek_gemeente_naam': 'bezoek'})
    .set_index(['woon','bezoek','datum'])
)

### Mezuro shapes
mezuro_shapes=gdp.read_file(mezuro_path('Mezuro shapefiles 2019/shapefiles/mezurogebieden_2019.shp')).set_index('mzr_id')
# The names of the mezuro areas are not unique, we append the name with the provincie name to make them unique.
counts = Counter(mezuro_shapes['mzr_name'])
mzr_unique = {
    mzr: mezuro_shapes.loc[mzr,'mzr_name'] if counts[mezuro_shapes.loc[mzr,'mzr_name']]==1 else "{} ({})".format(mezuro_shapes.loc[mzr,'mzr_name'],mezuro_shapes.loc[mzr,'prov_name'])
    for mzr in mezuro_shapes.index
}
mezuro_shapes['mzr_unique'] = mzr_unique.values()
# 'Overige gebieden' is not in the shapefiles
mzr_unique[9999] = 'Overige gebieden'

mezuro_shapes=mezuro_shapes.set_index('mzr_unique')

### Mezuro mobility
mezuro2mezuro=pd.read_csv(mezuro_path(
    'nederland/mezurogebied niveau/per dag/Herkomsten bezoekers per dag - van mezurogebied naar mezurogebied - 01-03-2019 tm 14-03-2019 - ilionx.csv'
), delimiter=';')
# Give unique names to mezuro areas
mezuro2mezuro['woon'] = [
    mzr_unique[woon_id]
    for woon_id in mezuro2mezuro['woon_mezurogebied_id'].values
]
mezuro2mezuro['bezoek'] = [
    mzr_unique[bezoek_id]
    for bezoek_id in mezuro2mezuro['bezoek_mezurogebied_id'].values
]
mezuro2mezuro=mezuro2mezuro.set_index(['woon','bezoek','datum'])

def create_koppeltabel():
    from collections import defaultdict
    code2018to2019 = pd.read_csv(
        data_path('koppeltabel_gemeenten_2018_20192020.csv'), index_col='gem_id_2018'
    )['gem_id_20192020'].to_dict()
    code2020to_name = pd.read_csv('COVID_tools/data/Gemeenten2020.csv', index_col='Gemeentecode')[
        'Gemeentenaam'].to_dict()
    koppeltabel = pd.DataFrame()
    koppeltabel['code2018'] = gemeente_shapes['gem_id']
    koppeltabel['code2019'] = [
        code2018to2019[c]
        for a, c in koppeltabel['code2018'].items()
    ]
    koppeltabel['rivm'] = [
        code2020to_name[c]
        for c in koppeltabel['code2019']
    ]
    # Collect inhabitant_frac
    inhabitant_old = gemeente_shapes['inhabitant'].to_dict()
    inhabitant_new = defaultdict(float)
    for a, inh in inhabitant_old.items():
        inhabitant_new[koppeltabel['rivm'][a]] += inh
    koppeltabel['inhabitant_frac_new'] = [
        inh_old / inhabitant_new[koppeltabel['rivm'][a]]
        for a, inh_old in inhabitant_old.items()
    ]
    koppeltabel.to_csv(data_path('koppeltabel.csv'))

create_koppeltabel()
koppeltabel = pd.read_csv(data_path('koppeltabel.csv'),index_col='name')

# Add a column 'veiligheidsregio' to gemeente_shapes
def get_veiligheidsregios():
    koppel_dict=koppeltabel['rivm'].to_dict()
    df_veiligheidsregio=pd.read_csv(data_path('Veiligheidsregios.csv'),delimiter=';')
    df_veiligheidsregio[' gemeentenaam']=df_veiligheidsregio[' gemeentenaam'].str.strip()
    veiligheidsregio=df_veiligheidsregio.set_index(' gemeentenaam')['veiligheidsregio'].to_dict()
    gemeente_shapes['veiligheidsregio'] = [
        veiligheidsregio[koppel_dict[g]]
        for g in gemeente_shapes.index
    ]
#get_veiligheidsregios()

druktebeeld_provincie=pd.read_csv(mezuro_path(
    'nederland/provincie niveau/per dag/Druktebeeld bewoners en bezoekers per dag - provincie - 01-03-2019 tm 14-03-2019 - ilionx.csv'
), delimiter=';', index_col='bezoek_provincie_naam')

druktebeeld_gemeente=pd.read_csv(mezuro_path(
    'nederland/gemeente niveau/per dag/Druktebeeld bewoners en bezoekers per dag - gemeente - 01-03-2019 tm 14-03-2019 - ilionx.csv'
), delimiter=';', index_col='bezoek_gemeente_naam')


mezuro_hierarchy = Division(gemeente_shapes['prov_name'].to_dict(),previouslevel=Division(mezuro_shapes['gem_name'].to_dict()))

