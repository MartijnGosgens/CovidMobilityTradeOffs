# COVID tools

The mobility data can be obtained by running `from COVID_tools.mobility import mobility`. 

The file `constants.py` contains the class `CoronaConstants`, which stores the parameters of the epidemiological model. By default, it will take the values used in the paper, but other parameter settings can be used by passing them to the constructor. For example,

```python
from COVID_tools.constants import CoronaConstants
constants = CoronaConstants(basic_reproduction_number=0.9)
```

sets the basic reproduction number ($R_0$) to $0.9$. All other related constants are automatically updated.

The epidemiologic model can be found in `mobility_SEIR.py`. The constructor of the `MobilitySEIR` class takes an optional input `constants`, which can be used to change the epidemiological parameters.

The initializations of this model that are used in the paper are stored in the folder `initializations/`. These were generated by 

```python
from COVID_tools.mobility_seir import generate_initializations
generate_initializations()
```

The folder `divisions/` contains all the divisions obtained by the various heuristics. The mobility regions and adaptive mobility regions were obtained by

```python
from COVID_tools.compute_heuristic_divisions import compute_mobility_regions,compute_adaptive_mobility_regions
compute_mobility_regions()
compute_adaptive_mobility_regions()
```

These heuristic divisions can be optimized by the script `optimize_divisions.py`. This optimization is computationally heavy.

```python
from COVID_tools.optimize_divisions import optimize_init
optimize_init('concentrated')
optimize_init('evenlydistributed')
optimize_init('historical0310')
optimize_init('historical0421')
```

All the figures that are used in the paper are be generated by `generate_figures.py`. This can be done by

```python
from COVID_tools.generate_figures import plot_all
plot_all(load_performances=True)
```

This will load all the performances of the divisions from the folder `results/`. When `load_performances` is set to `False`, the performances will be re-computed by simulating the model for each division.

