# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

subCh = pyna.rates.Library()

all_reactants = [(("he4", "he4", "he4"), ("c12")),
        (("c12", "he4"), ("o16"))]

for r, p in all_reactants:
    if not isinstance(p, tuple):
        p = p,
    rfilter = pyna.rates.RateFilter(reactants=r, products=p)
    _library = mylibrary.filter(rfilter)
    subCh += _library


net = StarKillerNetwork(libraries=[subCh])
net.write_network()
