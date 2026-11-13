"""Residue-matched backbone comparisons; scores are not optimised TM-scores."""

import re

import gemmi
import numpy as np
from scipy.spatial.distance import pdist, squareform


