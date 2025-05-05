print("Loading dependencies...")

from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy as sch

import warnings
warnings.filterwarnings("ignore")  # Suppress warnings

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 100  # Optional: Set nicer resolution for plots

from matplotlib import pyplot as plt

import scanpy as sc
sc.settings.verbosity = 0  # Suppress Scanpy logging

import squidpy as sq
import gcsfs

from clustergrammer2 import net, Network, CGM2

import zlib, json, base64
def json_zip(j):
    zip_json_string = base64.b64encode(
        zlib.compress(
            json.dumps(j).encode('utf-8')
        )
    ).decode('ascii')
    return zip_json_string

import os
import anndata as ad