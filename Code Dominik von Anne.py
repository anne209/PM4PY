
# Libraries importieren
import pm4py
import pandas as pd
import os
import matplotlib.pyplot as plt

# Funktion zum Aufrufen von .xes-Dateien
def import_xes(path_to_xes: str):
    """Liest ein XES-Log ein"""
    log = pm4py.read_xes(path_to_xes)
    return log # noch in df umwandeln?


# Aufruf des Sepsis-Datensatzes und Ausgabe einer Zusammenfassung
log_sepsis = import_xes('sepsis_case.xes')
log_sepsis

# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
def sum_up_log (log):
    print(log.describe())
    print('Head des Datensatzes: {}\nTail des Datensatzes: {}'.format(log.head(10), log.tail(10)))
    
    
sum_up_log(log_sepsis)

# Anzahl der Fälle und Ereignisse anzeigen
cases_sepsis = len(log_sepsis['case:concept:name'].unique())
events_sepsis = len(log_sepsis)
print('Anzahl Fälle: {}\nAnzahl Ereignisse: {}'.format(cases_sepsis, events_sepsis))

# Start- und Endaktivitäten anzeigen
def get_start_end_act(log):
    start_act = pm4py.get_start_activities(log)
    end_act = pm4py.get_end_activities(log)
    print('Startaktivitäten: {}\nEndaktivitäten: {}'.format(start_act, end_act))


get_start_end_act(log_sepsis)

# DFG erstellen
def create_dfg_from_log(log_input):
    dfg, start_activities, end_activities = pm4py.discover_dfg(log_input)
    dfg_output = pm4py.view_dfg(dfg, start_activities, end_activities)
    return dfg_output
    
dfg_sepsis_vor_filtern = create_dfg_from_log(log_sepsis)

# Filtern nach Startaktivitäten
def filter_by_start_activities(log, starts_to_keep):
    log_start = pm4py.get_start_activities(log)  # {activity: count}
    filtered_log_start  = pm4py.filter_start_activities(log, starts_to_keep)
    return log_start, filtered_log_start


log_start, filtered_log_start = filter_by_start_activities(log_sepsis, ['ER Registration'])

# Alpha Miner anwenden
def run_alpha_miner(log):
    net, initial_marking, final_marking = pm4py.discover_petri_net_alpha(log)
    pm4py.view_petri_net(net, initial_marking, final_marking)


run_alpha_miner(log_sepsis)

# Heuristic Miner anwenden
def run_heuristic_miner(log):
    heu_net = pm4py.discover_heuristics_net(log, dependency_threshold=0.99)
    pm4py.view_heuristics_net(heu_net)


run_heuristic_miner(log_sepsis)

# Inductive Miner anwenden
def run_inductive_miner(log):
    net_inductive, initial_marking_inductive, final_marking_inductive = pm4py.discover_petri_net_inductive(log)
    pm4py.view_petri_net(net_inductive, initial_marking_inductive, final_marking_inductive)


run_inductive_miner(log_sepsis)

# Varianten auflisten und danach filtern
def show_filter_vars(log):
    variants = pm4py.get_variants(log_sepsis)
    print('Varianten: {}'.format(variants))


def filter_by_vars(log):
    k = 5 # evtl. überflüssig?
    filtered_log_var = pm4py.filter_variants_top_k(log, k)
    return filtered_log_var


filtered_log_var_sepsis = filter_by_vars(log_sepsis)