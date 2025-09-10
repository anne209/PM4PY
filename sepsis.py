import pm4py
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn
import csv
from typing import List, Dict
import xml.etree.ElementTree as ET
from pm4py.algo.discovery.temporal_profile import algorithm as temporal_profile_discovery
from pm4py.algo.discovery.temporal_profile import algorithm as temporal_profile_discovery # doppelt
from pm4py.visualization.petri_net import visualizer as pn_visualizer # für Zusatzinfos Performance
from pm4py.algo.conformance.tokenreplay import algorithm as token_based_replay
from pm4py.algo.conformance.tokenreplay.diagnostics import duration_diagnostics # doppelt
from pm4py.algo.conformance.tokenreplay.diagnostics import root_cause_analysis
from pm4py.visualization.decisiontree import visualizer as dt_vis
from pm4py.algo.evaluation.generalization import algorithm as generalization_evaluator
from pm4py.algo.evaluation.simplicity import algorithm as simplicity_evaluator
from pm4py.util import constants
from pm4py.statistics.traces.generic.log import case_statistics
from pm4py.visualization.graphs import visualizer as graphs_visualizer
import Funktionen_Projekt as fkp

# Import des Datensatzes
log_sepsis = fkp.import_xes('sepsis_case.xes')
log_sepsis.to_csv('log_sepsis.csv', index=False)

fkp.sum_up_log(log_sepsis)

cases_no_sepsis = fkp.get_cases_events(log_sepsis)
