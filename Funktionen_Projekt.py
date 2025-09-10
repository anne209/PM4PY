# Libraries importieren
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

# Funktion zum Import von .xes-Dateien
def import_xes(path_to_xes: str):
    log = pm4py.read_xes(path_to_xes)
    return log

# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
def sum_up_log (log):
    print(f'Beschreibung des Datensatzes: \n{log.describe()}')
    print(f'Head des Datensatzes: \n{log.head(10)}\nTail des Datensatzes: \n{log.tail(10)}')

# Anzahl der Fälle und Ereignisse anzeigen
def get_cases_events(log):
    cases_no = len(log['case:concept:name'].unique())
    events_no = len(log)
    
    print(f'Anzahl Fälle: {cases_no}\nAnzahl Ereignisse: {events_no}')
    return cases_no

# Start- und Endaktivitäten anzeigen
def get_start_end_act(log):
    start_act = pm4py.get_start_activities(log)
    end_act = pm4py.get_end_activities(log)
    
    print(f'Startaktivitäten: {start_act}\nEndaktivitäten: {end_act}')
    return start_act, end_act