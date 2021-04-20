import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time
import pickle
import networkx as nx
import argparse

import rdkit.Chem as Chem
from rdkit.Chem import AllChem, QED, DataStructs
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def read_smilesset(path):
    smiles_list = []
    with open(path) as f:
        for smiles in f:
            smiles_list.append(smiles.rstrip())

    return smiles_list
