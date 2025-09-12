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


# BEGINN DES SEPSIS-DATENSATZES --------------------------------------------


# Import des Sepsis-Datensatzes
log_sepsis = fkp.import_xes('D:/Users/domin/Desktop/Projekt PM&Py/Datensätze/Sepsis/Sepsis Cases - Event Log.xes')
log_sepsis.to_csv('log_sepsis.csv', index=False)

# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
fkp.sum_up_log(log_sepsis)

# Anzahl der Cases und Events anzeigen
cases_no_sepsis = fkp.get_cases_events(log_sepsis)

# Start- und Endaktivitäten anzeigen
start_act_sepsis, end_act_sepsis = fkp.get_start_end_act(log_sepsis)

# DFG aus ungefiltertem Log erstellen und anzeigen
dfg_sepsis_unfiltered = fkp.create_dfg_from_log(log_sepsis)

# Filterkriterien definieren
criteria_end_sepsis = ['IV Antibiotics', 'ER Sepsis Triage', 'Leucocytes', 'IV Liquid', 'CRP', 'LacticAcid', 'Admission NC', 'ER Triage'] # Endaktivitäten, die herausgefiltert werden sollen
deleted_activities_sepsis = ['Return ER'] # Aktivitäten die generell gelöscht werden sollen im Log
activity_check_sepsis = ['Leucocytes', 'LacticAcid', 'CRP'] # Aktivitäten die auf ihre Vollständigkeit im Log überprüft werden sollen
column_filter_sepsis = ['lifecycle:transition'] # alle complete

filtered_log_sepsis = fkp.filter_log(start_act_sepsis, end_act_sepsis, log_sepsis, cases_no_sepsis, end_crit = criteria_end_sepsis, 
                                 delete_activities = deleted_activities_sepsis, check_value_activities=activity_check_sepsis, 
                                 col_filter=column_filter_sepsis) # Übergabe der Parameter an die Filter-Methode

# Beschreibung des gefilterten Logs inklusive Case- und Eventzahlen, Start- und Endaktivitäten
fkp.sum_up_log(filtered_log_sepsis)
fkp.get_cases_events(filtered_log_sepsis)
fkp.get_start_end_act(filtered_log_sepsis)

# DFG aus gefiltertem Log erstellen und anzeigen
dfg_sepsis_filtered = fkp.create_dfg_from_log(filtered_log_sepsis)

# Spaltenweises Zählen der Häufigkeit der Werte
fkp.get_column_count(filtered_log_sepsis, ['InfectionSuspected', 'org:group', 'DiagnosticBlood', 'DisfuncOrg', 'SIRSCritTachypnea', 'Hypotensie', 
                                       'SIRSCritHeartRate', 'Infusion', 'DiagnosticArtAstrup', 'DiagnosticIC', 'DiagnosticOther', 'SIRSCriteria2OrMore', 
                                       'DiagnosticXthorax', 'SIRSCritTemperature', 'DiagnosticUrinaryCulture', 'SIRSCritLeucos', 'Oligurie', 'DiagnosticLacticAcid',
                                       'Hypoxie', 'DiagnosticUrinarySediment', 'DiagnosticECG'])

# Erstellen eines Histogramms der Startaktivitäten
fkp.plot_start_activities(filtered_log_sepsis)

# Erstellen eines Histogramms der Endaktivitäten
fkp.plot_end_activities(filtered_log_sepsis)

# Erstellen eines Diagramms der Aktivitätshäufigkeiten
fkp.plot_activity_frequencies(filtered_log_sepsis)

# Erstellen eines Histogramms der Durchlaufzeiten der Cases
fkp.plot_case_duration_hist(filtered_log_sepsis)

# Erstellen eines Diagramms der Events pro Tag
fkp.plot_events_per_day(filtered_log_sepsis)

# Alpha Miner anwenden, Modell erstellen und anzeigen
alpha_net_sepsis, initial_marking_alpha_sepsis, final_marking_alpha_sepsis = fkp.run_alpha_miner(filtered_log_sepsis)

# Heuristic Miner anwenden, Modell anzeigen und Umwandlung in Petri-Netz
heur_net_sepsis, initial_marking_heur_sepsis, final_marking_heur_sepsis = fkp.run_heuristic_miner(filtered_log_sepsis, dependency_threshold=0.2)

# Inductive Miner für Sepsis-Datensatz anwenden, Modell erstellen und anzeigen
ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis = fkp.run_inductive_miner_sepsis(filtered_log_sepsis)

# Varianten auflisten
fkp.show_filter_variants(filtered_log_sepsis)

# Filtern und Ausgeben eines Logs mit den 5 häufigsten Varianten bzw. der häufigsten Variante, jeweils Visualisierung mit Inductive Miner und als DFG
for i in [1, 5]:
    filtered_log_var_sepsis = fkp.filter_by_variants(filtered_log_sepsis, i)
    fkp.run_inductive_miner_sepsis(filtered_log_var_sepsis)
    fkp.create_dfg_from_log(filtered_log_var_sepsis)

# Mittlere Durchlaufzeit pro Case berechnen
mean_throughput_sepsis = fkp.get_mean_throughput(filtered_log_sepsis)

# Performance-DFG erzeugen
fkp.get_performance_dfg(filtered_log_sepsis)

# Performance-Informationen zu Netz aus Inductive Miner hinzufügen
fkp.get_performance_net(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# Process Tree mit Inductive Miner erstellen und anzeigen
fkp.get_process_tree_ind(filtered_log_sepsis)

# TBR mit Modell aus Inductive Miner durchführen
df_tbr_sepsis = fkp.tbr_ind(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# Visualisierungen zum TBR erstellen (Histogramm der Fitness und Balkenplot fit/unfit)
fkp.plot_tbr_fitness_hist(df_tbr_sepsis)
fkp.plot_tbr_fit_flag(df_tbr_sepsis)

# TBR mit angepassten Einstellungen zur Vorbereitung der Throughput Analysis durchführen
log_diagnostics_sepsis = pm4py.convert_to_event_log(filtered_log_sepsis)
replayed_traces_sepsis, place_fitness_sepsis, trans_fitness_sepsis, unwanted_activities_sepsis = fkp.tbr_throughput(log_diagnostics_sepsis, heur_net_sepsis, initial_marking_heur_sepsis, final_marking_heur_sepsis)

# Throughput Analysis der falsch ausgeführten Transitionen durchführen und diese ausgeben
fkp.througput_trans(log_diagnostics_sepsis, trans_fitness_sepsis)

# Throughput Analysis der nicht im Modell enthaltenen Aktivitäten durchführen und diese ausgeben
fkp.throughput_act(log_diagnostics_sepsis, unwanted_activities_sepsis)

# Root Cause Analysis der falsch ausgeführten Transitionen durchführen und Decision Trees erstellen
fkp.rca_trans(log_diagnostics_sepsis, trans_fitness_sepsis)

# RCA der ausgeführten Aktivitäten, die nicht im Prozessmodell enthalten sind, durchführen und Decision Trees erstellen
fkp.rca_act(log_diagnostics_sepsis, unwanted_activities_sepsis)

# Alignments bestimmen mit Modell aus Inductive Miner und ausgeben
aligned_traces_ind_sepsis = fkp.alignments_inductive(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# Visualisierung der Alignments erstellen (Histogramm der Alignment-Kosten und Move-Typen)
df_stats_sepsis = fkp.extract_alignment_stats(aligned_traces_ind_sepsis)
fkp.plot_alignment_cost_hist(df_stats_sepsis)
fkp.plot_move_type_bars(df_stats_sepsis)

# Fitness zwischen Log und Modell berechnen (TBR und Alignments)
fitness_tbr_sepsis = fkp.get_replay_fitness_tbr(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)
fitness_align_sepsis = fkp.get_replay_fitness_align(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# Precision zwischen Log und Modell berechnen (TBR und Alignments(dauert sehr lange))
precision_tbr_sepsis = fkp.get_precision_tbr(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)
# precision_align_sepsis = fkp.get_precision_align(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis) # dauert zu lange

# F-Score zwischen Log und Modell berechnen, Fitness und Precision gegeben
f_score_tbr_sepsis = fkp.get_f_score(fitness_tbr_sepsis, precision_tbr_sepsis)
# f_score_align_sepsis = fkp.get_f_score(fitness_align_sepsis, precision_align_sepsis) # nicht möglich (s.o.)

# Generalization und Simplicity berechnen
fkp.get_generalization(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)
fkp.get_simplicity(ind_net_sepsis)

# Spalte bei Sepsis-Datensatz umbenennen für Social Network Analysis
log_sepsis_sna = filtered_log_sepsis.rename(columns = {'org:group' : 'org:resource'})

# Übergabe von Arbeit ermitteln und anzeigen
handover_sepsis = fkp.get_handover_of_work(log_sepsis_sna)

# Ermitteln und Anzeigen, wie oft Subcontracting vorkommt
fkp.get_subcontracting(log_sepsis_sna)

# Ermitteln und Anzeigen, wie oft zusammengearbeitet wird
fkp.get_working_together(log_sepsis_sna)

# Ähnlichkeiten der Arbeitsmuster zwischen Individuen ermitteln und anzeigen
fkp.get_similar_activities(log_sepsis_sna)

# Orginisationale Rollen entdecken und ausgeben
fkp.get_orga_roles(log_sepsis_sna)

# Cluster-Analyse nach Übergabe von Arbeit durchführen und anzeigen
fkp.cluster_handover(handover_sepsis)

# Herausfiltern aller Cases des Sepsis-Datensatzes, die mit ER-Registration beginnen (zur Überprüfung der eigenen Hypothese)
filtered_log_start_sepsis = fkp.filter_by_start_activities(log_sepsis, start_act_sepsis, ['ER Registration'])
fkp.sum_up_log(filtered_log_start_sepsis)
fkp.get_cases_events(filtered_log_start_sepsis)


'''# Funktionen aufrufen für nach Startaktivitäten gefiltertes Event Log Sepsis
heur_net_sepsis_start, initial_marking_heur_sepsis_start, final_marking_heur_sepsis_start = fkp.run_heuristic_miner(filtered_log_start_sepsis)

ind_net_sepsis_start, initial_marking_ind_sepsis_start, final_marking_ind_sepsis_start = fkp.run_inductive_miner_sepsis(filtered_log_start_sepsis)

mean_throuhput_start_sepsis = fkp.get_mean_throughput(filtered_log_start_sepsis)

fkp.get_performance_dfg(filtered_log_start_sepsis)

fkp.get_performance_net(filtered_log_start_sepsis, ind_net_sepsis_start, initial_marking_ind_sepsis_start, final_marking_ind_sepsis_start)

fkp.plot_start_activities(filtered_log_start_sepsis)

fkp.plot_end_activities(filtered_log_start_sepsis)

fkp.plot_activity_frequencies(filtered_log_start_sepsis)

fkp.plot_case_duration_hist(filtered_log_start_sepsis)

log_diagnostics_sepsis_start = pm4py.convert_to_event_log(filtered_log_start_sepsis)

replayed_traces_sepsis_start, place_fitness_sepsis_start, trans_fitness_sepsis_start, unwanted_activities_sepsis_start = fkp.tbr_throughput(log_diagnostics_sepsis_start, ind_net_sepsis_start, initial_marking_ind_sepsis_start, final_marking_ind_sepsis_start)

fkp.througput_trans(log_diagnostics_sepsis_start, trans_fitness_sepsis_start)

fkp.throughput_act(log_diagnostics_sepsis_start, unwanted_activities_sepsis_start)

fkp.rca_trans(log_diagnostics_sepsis_start, trans_fitness_sepsis_start)

fkp.rca_act(log_diagnostics_sepsis_start, unwanted_activities_sepsis_start)'''



# BEGINN DES IACS-DATENSATZES ---------------------------------------------


# Import des IACS-Datensatzes
log_iacs = fkp.import_xes('D:/Users/domin/Desktop/Projekt PM&Py/Datensätze/IACS/BPI Challenge 2018 (x0.05).xes') # Datei aktuell nicht hochgeladen
log_iacs.to_csv('log_iacs.csv', index=False)

# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
fkp.sum_up_log(log_iacs)

# Anzahl der Cases und Events anzeigen
cases_no_iacs = fkp.get_cases_events(log_iacs)

# Start- und Endaktivitäten anzeigen
start_act_iacs, end_act_iacs = fkp.get_start_end_act(log_iacs)

# DFG aus ungefiltertem Log erstellen
dfg_iacs_unfiltered = fkp.create_dfg_from_log(log_iacs)

# Filterkriterien definieren
column_filter_iacs = ['lifecycle:transition', 'case:program-id', 'case:penalty_BGKV', 'case:penalty_BGP', 'case:penalty_AVUVP', 
                      'case:greening', 'case:basic payment', 'case:penalty_B5F', 'case:penalty_JLP7', 
                      'case:penalty_JLP5', 'activity'] # complete, 215, False, False, False, True, True, False, False, False
criteria_end_iacs = ['remove document', 'begin editing', 'refuse', 'decide', 'calculate', 'insert document', 'withdraw', 
                     'change department', 'check admissibiliy', 'check', 'prepare offline', 'initialize', 'finish pre-check', 
                     'take original document', 'calculate protocol', 'restart editing', 'abort payment', 'begin admissibility check']

filtered_log_iacs = fkp.filter_log(start_act_iacs, end_act_iacs, log_iacs, cases_no_iacs, end_crit=criteria_end_iacs, col_filter=column_filter_iacs)

# Beschreibung des gefilterten Logs inklusive Case- und Eventzahlen, Start- und Endaktivitäten
fkp.sum_up_log(filtered_log_iacs)
fkp.get_cases_events(filtered_log_iacs)
fkp.get_start_end_act(filtered_log_iacs)

# DFG aus gefiltertem Log erstellen
dfg_iacs_filtered = fkp.create_dfg_from_log(filtered_log_iacs)

# Spaltenweises Zählen der Werte
fkp.get_column_count(filtered_log_iacs, ['success', 'doctype', 'subprocess', 'note', 'case:young farmer', 'case:selected_random',
                                     'case:penalty_AJLP', 'case:penalty_BGKV', 'case:penalty_AUVP', 'case:risk_factor', 'case:small farmer', # Spalten noch anpassen (dazu/weg)
                                     'case:penalty_BGP', 'case:department', 'case:penalty_C16', 'case:penalty_BGK', 'case:penalty_AVUVP',
                                     'case:penalty_CC', 'case:penalty_AVJLP', 'case:penalty_C9', 'case:cross_compliance', 'case:rejected',
                                     'case:greening', 'case:penalty_C4', 'case:penalty_AVGP', 'case:penalty_ABP', 'case:penalty_B6', 'case:penalty_B4',
                                     'case:penalty_B5', 'case:penalty_AVBP', 'case:penalty_B2', 'case:selected_risk', 'case:penalty_B3',
                                     'case:selected_manually', 'case:penalty_AGP', 'case:penalty_B16', 'case:penalty_GP1', 'case:basic payment',
                                     'case:penalty_B5F', 'case:penalty_V5', 'case:redistribution', 'case:penalty_JLP6', 'case:penalty_JLP7', 'case:year',
                                     'case:penalty_JLP5', 'case:penalty_JLP2', 'case:penalty_JLP3', 'case:number_parcels', 'case:penalty_JLP1'])

# Erstellen eines Histogramms der Startaktivitäten
fkp.plot_start_activities(filtered_log_iacs)

# Erstellen eines Histogramms der Endaktivitäten
fkp.plot_end_activities(filtered_log_iacs)

# Erstellen eines Diagramms der Aktivitätshäufigkeiten
fkp.plot_activity_frequencies(filtered_log_iacs)

# Erstellen eines Histogramms der Durchlaufzeiten der Cases
fkp.plot_case_duration_hist(filtered_log_iacs)

# Erstellen eines Diagramms der Events pro Tag
fkp.plot_events_per_day(filtered_log_iacs)

# Alpha Miner anwenden, Modell erstellen und anzeigen
alpha_net_iacs, initial_marking_alpha_iacs, final_marking_alpha_iacs = fkp.run_alpha_miner(filtered_log_iacs)

# Heuristic Miner anwenden, Modell anzeigen und Umwandlung in Petri-Netz
heur_net_iacs, initial_marking_heur_iacs, final_marking_heur_iacs = fkp.run_heuristic_miner(filtered_log_iacs, dependency_threshold=0.0)

# Inductive Miner für IACS-Datensatz anwenden, Modell erstellen und anzeigen
ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs = fkp.run_inductive_miner_iacs(filtered_log_iacs)

# Varianten auflisten
fkp.show_filter_variants(filtered_log_iacs)

# Filtern und Ausgeben eines Logs mit den 5 häufigsten Varianten bzw. der häufigsten Variante sowie jeweils Visualisierung mit Inductive Miner und als DFG
for i in [1, 5]:
    filtered_log_var_iacs = fkp.filter_by_variants(filtered_log_iacs, i)
    fkp.run_inductive_miner_iacs(filtered_log_var_iacs)
    fkp.create_dfg_from_log(filtered_log_var_iacs)

# Mittlere Durchlaufzeit pro Case berechnen
mean_throughput_iacs = fkp.get_mean_throughput(filtered_log_iacs)

# Performance-DFG erzeugen
fkp.get_performance_dfg(filtered_log_iacs)

# Performance-Informationen zu Netz aus Inductive Miner hinzufügen
fkp.get_performance_net(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Process Tree mit Inductive Miner erstellen und anzeigen
fkp.get_process_tree_ind(filtered_log_iacs)

# Temporal Profile erstellen

temp_prof_iacs = filtered_log_iacs.copy()
temp_prof_iacs.insert(loc=temp_prof_iacs.columns.get_loc('time:timestamp'), column='start_timestamp', value=None)

temp_prof_iacs = fkp.merge_begin_finish(temp_prof_iacs, 'begin editing', 'finish editing', 'editing')
temp_prof_iacs = fkp.merge_begin_finish(temp_prof_iacs, 'begin preparations', 'finish preparations', 'preparations')
temp_prof_iacs = fkp.merge_begin_finish(temp_prof_iacs, 'begin payment', 'finish payment', 'payment', allow_abort=True)

temp_prof_iacs = pm4py.filter_event_attribute_values(temp_prof_iacs, 'concept:name', ['begin editing', 'finish editing', 
                                                                                      'begin preperations', 'finish preparations', 
                                                                                      'begin payment', 'finish payment', 
                                                                                      'abort payment'], level='event', retain=False)

fkp.get_temporal_profile(temp_prof_iacs)

# TBR mit Modell aus Inductive Miner durchführen
df_tbr_iacs = fkp.tbr_ind(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Visualisierungen zum TBR erstellen (Histogramm der Fitness und Balkenplot fit/unfit)
fkp.plot_tbr_fitness_hist(df_tbr_iacs)
fkp.plot_tbr_fit_flag(df_tbr_iacs)

# TBR mit angepassten Einstellungen zur Vorbereitung der Throughput Analysis durchführen
log_diagnostics_iacs = pm4py.convert_to_event_log(filtered_log_iacs)
replayed_traces_iacs, place_fitness_iacs, trans_fitness_iacs, unwanted_activities_iacs = fkp.tbr_throughput(filtered_log_iacs, heur_net_iacs, initial_marking_heur_iacs, final_marking_heur_iacs)

# Throughput Analysis der falsch ausgeführten Transitionen durchführen und diese ausgeben
fkp.througput_trans(filtered_log_iacs, trans_fitness_iacs)

# Throughput Analysis der nicht im Modell enthaltenen Aktivitäten durchführen und diese ausgeben
fkp.throughput_act(filtered_log_iacs, unwanted_activities_iacs)

# Root Cause Analysis der falsch ausgeführten Transitionen durchführen und Decision Trees erstellen
fkp.rca_trans(filtered_log_iacs, trans_fitness_iacs)

# RCA der ausgeführten Aktivitäten, die nicht im Prozessmodell enthalten sind, durchführen und Decision Trees erstellen
fkp.rca_act(filtered_log_iacs, unwanted_activities_iacs)

# Alignments bestimmen und ausgeben mit Modell aus Inductive Miner (dauert zu lange mit vollständigem Datensatz)
# aligned_traces_ind_iacs = fkp.alignments_inductive(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Visualisierung der Alignments erstellen (Histogramm der Alignment-Kosten und Move-Typen), nicht möglich für IACS
# df_stats_iacs = fkp.extract_alignment_stats(aligned_traces_ind_iacs)
# fkp.plot_alignment_cost_hist(df_stats_iacs)
# fkp.plot_move_type_bars(df_stats_iacs)

# Fitness zwischen Log und Modell berechnen (mit TBR und Alignments)
fitness_tbr_iacs = fkp.get_replay_fitness_tbr(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)
# fitness_align_iacs = fkp.get_replay_fitness_align(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Precision zwischen Log und Modell berechnen (mit TBR und Alignments)
precision_tbr_iacs = fkp.get_precision_tbr(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)
# precision_align_iacs = fkp.get_precision_align(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# F-Score zwischen Log und Modell berechnen, Fitness und Precision gegeben
f_score_tbr_iacs = fkp.get_f_score(fitness_tbr_iacs, precision_tbr_iacs)
# f_score_align_iacs = fkp.get_f_score(fitness_align_iacs, precision_align_iacs)

# Generalization und Simplicity berechnen
fkp.get_generalization(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)
fkp.get_simplicity(ind_net_iacs)

# Übergabe von Arbeit ermitteln und anzeigen
handover_iacs = fkp.get_handover_of_work(filtered_log_iacs)

# Ermitteln und Anzeigen, wie oft Subcontracting vorkommt
fkp.get_subcontracting(filtered_log_iacs)

# Ermitteln und Anzeigen, wie oft zusammengearbeitet wird
fkp.get_working_together(filtered_log_iacs)

# Ähnlichkeiten der Arbeitsmuster zwischen Individuen ermitteln und anzeigen
similar_activities_iacs = fkp.get_similar_activities(filtered_log_iacs)

# Orginisationale Rollen entdecken und ausgeben
fkp.get_orga_roles(filtered_log_iacs)

# Cluster-Analyse nach Übergabe von Arbeit durchführen und anzeigen
# fkp.cluster_handover(handover_iacs) # Kein sinnvoller Output