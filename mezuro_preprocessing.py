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
    # Take koppeltabel
    koppeltabel=pd.read_csv(
        data_path('koppeltabel_gemeenten_2018_20192020.csv')
    )[['gem_id_2018','gem_id_20192020']].set_index('gem_id_2018')

    # translate to names instead of indices
    id2name=gemeente_shapes.reset_index().set_index('gem_id',drop=False)['name'].to_dict()
    koppeltabel['name_old'] = [
        id2name[i]
        for i in koppeltabel.index
    ]
    name2newid=koppeltabel.set_index('name_old')['gem_id_20192020'].to_dict()
    # Obtained from https://nl.wikipedia.org/wiki/Gemeentelijke_herindelingen_in_Nederland
    id2name.update({
        name2newid["'s-Gravenhage"]: 's-Gravenhage',
        name2newid['Bedum']: 'Het Hogeland',
        name2newid['Dongeradeel']: 'Noardeast-Fryslân',
        name2newid['Geldermalsen']: 'West Betuwe',
        name2newid['Leek']: 'Westerkwartier',
        name2newid['Winsum']: 'Westerkwartier',
        name2newid['Nuth']: 'Beekdaelen',
        name2newid['Aalburg']: 'Altena',
        name2newid['Leerdam']: 'Vijfheerenlanden',
        name2newid['Binnenmaas']: 'Hoeksche Waard',
        name2newid['Giessenlanden']: 'Molenlanden',
    })
    koppeltabel['name_new'] = [
        id2name[i]
        for i in koppeltabel['gem_id_20192020']
    ]
    old2new=koppeltabel.set_index('name_old')['name_new'].to_dict()
    gemeente_shapes['rivm'] = [
        old2new[i]
        for i in gemeente_shapes.index
    ]

    # Compute the fraction of the inhabitants that the old gemeente has in the new gemeente
    inhabitant_new = gemeente_shapes.groupby('rivm').agg({'inhabitant': 'sum'})['inhabitant'].to_dict()
    gemeente_shapes['inhabitant_frac_new'] = [
        row['inhabitant'] / inhabitant_new[row['rivm']]
        for g,row in gemeente_shapes.iterrows()
    ]
    gemeente_shapes[['rivm','inhabitant_frac_new']].to_csv(data_path('koppeltabel.csv'))

# create_koppeltabel()
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
get_veiligheidsregios()

druktebeeld_provincie=pd.read_csv(mezuro_path(
    'nederland/provincie niveau/per dag/Druktebeeld bewoners en bezoekers per dag - provincie - 01-03-2019 tm 14-03-2019 - ilionx.csv'
), delimiter=';', index_col='bezoek_provincie_naam')

druktebeeld_gemeente=pd.read_csv(mezuro_path(
    'nederland/gemeente niveau/per dag/Druktebeeld bewoners en bezoekers per dag - gemeente - 01-03-2019 tm 14-03-2019 - ilionx.csv'
), delimiter=';', index_col='bezoek_gemeente_naam')


mezuro_hierarchy = Division(gemeente_shapes['prov_name'].to_dict(),previouslevel=Division(mezuro_shapes['gem_name'].to_dict()))

