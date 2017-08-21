#Module to perform power-law fits to wp(rp)

import numpy as np
import emcee
import matplotlib.pyplot as plt
from scipy.special import gamma as gammafn

from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import AssembiasZheng07Cens
from halotools.empirical_models import TrivialPhaseSpace
from halotools.empirical_models import AssembiasZheng07Sats
from halotools.empirical_models import NFWPhaseSpace
from halotools.mock_observables import return_xyz_formatted_array
from halotools.mock_observables import wp
from halotools.sim_manager import CachedHaloCatalog

from Corrfunc import _countpairs

#default settings of all parameters
from abhodfit_defaults import *

#A class to hold the attributes of a power-law model for wp(rp)
class ABHodFitModel():
#HOD with assembly bias fit model class for wp(rp)

##############################################3
    #initialize an instance of the ABHodFitModel
    def __init__(self,**kwargs):
        if 'priors' in kwargs.keys():
            self.set_prior(kwargs['priors'])
        else:
            self.set_prior(default_priors)

        self.param_names=
