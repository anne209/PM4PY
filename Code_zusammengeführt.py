
# Libraries importieren
import pm4py
import pandas as pd
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


#HEADER = b'''<?xml version="1.0" encoding="UTF-8" ?>
#<!-- This file has been generated with the OpenXES library. It conforms -->
#<!-- to the XML serialization of the XES standard for log storage and -->
#<!-- management. -->
#<!-- XES standard version: 1.0 -->
#<!-- OpenXES library version: 1.0RC7 -->
#<!-- OpenXES is available from http://www.openxes.org/ -->
#<log xes.version="1.0" xes.features="nested-attributes" openxes.version="1.0RC7">
#	<extension name="Organizational" prefix="org" uri="http://www.xes-standard.org/org.xesext"/>
#	<extension name="MetaData_Time" prefix="meta_time" uri="http://www.xes-standard.org/meta_time.xesext"/>
#	<extension name="Time" prefix="time" uri="http://www.xes-standard.org/time.xesext"/>
#	<extension name="MetaData_3TU" prefix="meta_3TU" uri="http://www.xes-standard.org/meta_3TU.xesext"/>
#	<extension name="Lifecycle" prefix="lifecycle" uri="http://www.xes-standard.org/lifecycle.xesext"/>
#	<extension name="Concept" prefix="concept" uri="http://www.xes-standard.org/concept.xesext"/>
#	<global scope="trace">
#		<boolean key="penalty_JLP1" value="false"/>
#		<int key="number_parcels" value="0"/>
#		<boolean key="penalty_JLP3" value="false"/>
#		<boolean key="penalty_JLP2" value="false"/>
#		<boolean key="penalty_JLP5" value="false"/>
#		<string key="year" value="none"/>
#		<boolean key="penalty_JLP7" value="false"/>
#		<boolean key="penalty_JLP6" value="false"/>
#		<boolean key="redistribution" value="false"/>
#		<float key="amount_applied0" value="0.0"/>
#		<id key="identity:id" value="AAE527A9-1929-49A2-AB4A-147CD5D2ADD5"/>
#		<float key="payment_actual0" value="0.0"/>
#		<boolean key="penalty_V5" value="false"/>
#		<string key="concept:name" value="none"/>
#		<boolean key="penalty_B5F" value="false"/>
#		<boolean key="basic payment" value="false"/>
#		<boolean key="penalty_GP1" value="false"/>
#		<boolean key="penalty_B16" value="false"/>
#		<boolean key="penalty_AGP" value="false"/>
#		<boolean key="selected_manually" value="false"/>
#		<float key="area" value="0.0"/>
#		<boolean key="penalty_B3" value="false"/>
#		<boolean key="selected_risk" value="false"/>
#		<boolean key="penalty_B2" value="false"/>
#		<boolean key="penalty_AVBP" value="false"/>
#		<boolean key="penalty_B5" value="false"/>
#		<boolean key="penalty_B4" value="false"/>
#		<boolean key="penalty_B6" value="false"/>
#		<boolean key="penalty_ABP" value="false"/>
#		<boolean key="penalty_AVGP" value="false"/>
#		<boolean key="penalty_C4" value="false"/>
#		<boolean key="greening" value="false"/>
#		<boolean key="rejected" value="false"/>
#		<float key="cross_compliance" value="0.0"/>
#		<boolean key="penalty_C9" value="false"/>
#		<boolean key="penalty_AVJLP" value="false"/>
#		<boolean key="penalty_CC" value="false"/>
#		<boolean key="penalty_AVUVP" value="false"/>
#		<boolean key="penalty_BGK" value="false"/>
#		<boolean key="penalty_C16" value="false"/>
#		<string key="department" value="none"/>
#		<boolean key="penalty_BGP" value="false"/>
#		<boolean key="small farmer" value="false"/>
#		<float key="risk_factor" value="0.0"/>
#		<string key="applicant" value="none"/>
#		<boolean key="penalty_AUVP" value="false"/>
#		<boolean key="penalty_BGKV" value="false"/>
#		<string key="program-id" value="none"/>
#		<float key="penalty_amount0" value="0.0"/>
#		<string key="application" value="none"/>
#		<boolean key="penalty_AJLP" value="false"/>
#		<boolean key="selected_random" value="false"/>
#		<boolean key="young farmer" value="false"/>
#	</global>
#	<global scope="event">
#		<string key="eventid" value="none"/>
#		<string key="note" value="none"/>
#		<string key="activity" value="none"/>
#		<string key="docid" value="none"/>
#		<string key="subprocess" value="none"/>
#		<date key="time:timestamp" value="1970-01-01T01:00:00.000+01:00"/>
#		<id key="identity:id" value="9C25102D-21F3-4361-84C7-DA8BBF6953BC"/>
#		<string key="doctype" value="none"/>
#		<id key="docid_uuid" value="5302F7BA-6EA6-4E69-880B-C9C203EB86C9"/>
#		<string key="org:resource" value="none"/>
#		<string key="concept:name" value="none"/>
#		<boolean key="success" value="false"/>
#		<string key="lifecycle:transition" value="none"/>
#	</global>
#	<classifier name="Event Name" keys="concept:name"/>
#	<classifier name="Original Activity" keys="doctype subprocess activity"/>
#	<string key="meta_3TU:language" value="eng"/>
#	<float key="meta_time:duration_average" value="2.8977602621E7"/>
#	<date key="meta_3TU:creation_time" value="2018-02-07T11:49:52.625+01:00"/>
#	<float key="meta_time:log_duration" value="1.17208982E8"/>
#	<string key="meta_3TU:source_institute_type" value="Document management company"/>
#	<string key="meta_3TU:rights" value="Public"/>
#	<string key="meta_3TU:creation_place" value="Germany"/>
#	<float key="meta_time:duration_max" value="8.7385045817E7"/>
#	<string key="meta_3TU:source_program" value="profil c/s"/>
#	<id key="identity:id" value="86E187BA-C776-4459-B446-9C61FA2AF98C"/>
#	<string key="meta_3TU:creation_person" value="B.F. van Dongen &amp; Florian Borchert"/>
#	<float key="meta_time:duration_standard_deviation" value="1.3734464217E7"/>
#	<string key="meta_3TU:description" value="The European Union spends a large fraction of its budget on the Common Agricultural Policy (CAP). Among these spendings are direct payments, which are mainly aimed to provide a basic income for farmers decoupled from production. The rest of the CAP budget is spent for market related expenditures and rural development. The processes that govern the distribution of these funds are subject to complex regulations captured in EU and national law. The member states are required to operate an Integrated Administration and Control System (IACS), which includes IT systems to support the complex processes of subsidy distribution. The process considered in this dataset covers the handling of applications for EU direct payments for German farmers from the European Agricultural Guarantee Fund. The process repeats every year with minor changes due to changes in EU regulations. The dataset is extracted from the systems of data experts, Germany. Their tool profil c/s supports these processes at the level of federal ministries of agriculture and local departments. The workflows in profil c/s can be understood in terms of documents, where each document has a state that allows for certain actions. These actions can be executed manually at any point in time through document specific tools or they can be scheduled automatically. The latter may be either explicitly stated in the log or implicitly apparent if a large number of actions is performed by the same user at around the same time (batch processing). In total, the event log contains 2,514,266 events for 43,809 applications over a period of three years. The shortest case contains 24 events, the longest 2973 and on average there are 57 events per case referring to 14 activities."/>
#	<string key="meta_3TU:creation_institute" value="Eindhoven University of Technology"/>
#	<string key="concept:name" value="BPI Challenge 2018"/>
#	<date key="meta_time:log_start_time" value="2014-05-04T00:00:00.000+02:00"/>
#	<float key="meta_time:duration_total" value="1.269479793225363E12"/>
#	<string key="meta_3TU:source_institute" value="Eindhoven University of Technology &amp; data experts, Neubrandenburg, Germany"/>
#	<string key="meta_3TU:process_type" value="Implicitly structured"/>
#	<float key="meta_time:duration_min" value="9828893.063"/>
#	<string key="meta_3TU:log_type" value="Real-life"/>
#	<date key="meta_time:log_end_time" value="2018-01-19T13:03:02.000+01:00"/>
#	<string key="meta_3TU:source_model" value="EU Common Agricultural Policy"/>
#	<string key="meta_3TU:doi" value="doi:10.4121/uuid:3301445f-95e8-4ff0-98a4-901f1f204972"/>'''


#def clean_xes(input_file, output_file, remove_keys):
#    with open(output_file, 'wb') as out:
#        out.write(HEADER)
        
#        context = ET.iterparse(input_file, events=("end",))
#        for event, elem in context:
#            if elem.tag in ("trace", "event"):
#                for child in list(elem):
#                    if child.tag in ("string", "date", "int", "float", "boolean", "id"):
#                        key = child.attrib.get("key")
#                        if key in remove_keys:
#                            elem.remove(child)

#                out.write(ET.tostring(elem, encoding="utf-8"))
#                elem.clear()

        # Am Ende die Log-Tag schließen
#        out.write(b"</log>")


#remove_keys = {"lifecycle:transition", "case:program-id", "case:penalty_BGKV", "case:penalty_BGP", "case:penalty_AVUVP", 
#                      "case:greening", "case:basic payment", "case:penalty_B5F", "case:penalty_JLP7", 
#                      "case:penalty_JLP5"}
#clean_xes("log_iacs.xes", "flog_iacs.xes", remove_keys)


# Funktion zum Import von .xes-Dateien
def import_xes(path_to_xes: str):
    log = pm4py.read_xes(path_to_xes)
    return log

# Import des Datensatzes
log_sepsis = import_xes('sepsis_case.xes')
log_sepsis.to_csv('log_sepsis.csv', index=False)

log_iacs = import_xes('BPI Challenge 2018 (x0.05).xes') # Datei aktuell nicht hochgeladen
log_iacs.to_csv('log_iacs.csv', index=False)


# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
def sum_up_log (log):
    print(f'Beschreibung des Datensatzes: \n{log.describe()}')
    print(f'Head des Datensatzes: \n{log.head(10)}\nTail des Datensatzes: \n{log.tail(10)}')
    
    
sum_up_log(log_sepsis)

sum_up_log(log_iacs)

# Anzahl der Fälle und Ereignisse anzeigen
def get_cases_events(log):
    cases_no = len(log['case:concept:name'].unique())
    events_no = len(log)
    
    print(f'Anzahl Fälle: {cases_no}\nAnzahl Ereignisse: {events_no}')
    return cases_no

cases_no_sepsis = get_cases_events(log_sepsis)

cases_no_iacs = get_cases_events(log_iacs)

# Start- und Endaktivitäten anzeigen
def get_start_end_act(log):
    start_act = pm4py.get_start_activities(log)
    end_act = pm4py.get_end_activities(log)
    
    print(f'Startaktivitäten: {start_act}\nEndaktivitäten: {end_act}')
    return start_act, end_act

start_act_sepsis, end_act_sepsis = get_start_end_act(log_sepsis)

start_act_iacs, end_act_iacs = get_start_end_act(log_iacs)

# DFG aus ungefiltertem Log erstellen
def create_dfg_from_log(log_input):
    dfg, start_activities, end_activities = pm4py.discover_dfg(log_input)
    dfg_output = pm4py.view_dfg(dfg, start_activities, end_activities)
    return dfg_output
    

dfg_sepsis_unfiltered = create_dfg_from_log(log_sepsis)

dfg_iacs_unfiltered = create_dfg_from_log(log_iacs)


# Filteralgorithmus
def filter_log(start_acts, end_acts, log, no_of_cases, min_ratio=0.1, end_crit = None, delete_activities = None, check_value_activities=None, col_filter=None):
    selected_activities = set()
    selected_end_acts = set()
    filtered_log = log

    for activity, count in start_acts.items(): # Iteration im Dict über den Namen der Startaktivität und ihre Häufigkeit.
        if count / no_of_cases >= min_ratio: # Nur Startaktivitäten mit einem Mindestanteil von 10% an Gesamtfällen werden berücksichtigt.
            selected_activities.add(activity) # Sammeln der Aktivitäten in einer Liste
    
    if end_crit is not None:
        selected_end_acts = set(end_acts.keys()) - set(end_crit) # Löschen Fällen, die nicht ordungsgemäß mit Release oder Return ER geendet haben.
    
    
    if check_value_activities is not None:

        drop_mask = pd.Series(False, index=filtered_log.index) # Erstellen einer Maske mit Boolean-Werte

        for item in check_value_activities: # Sepsisspezifisch!! Es wird überprüft, ob die Laborwerte tatsächlich vorhanden sind, oder ob nur die Aktivität ohne Laborwert ausgeführt wurde

            is_activity = filtered_log['concept:name'] == item # True, wenn die Aktivität ausgeführt wrude
            is_empty = filtered_log[item].isna() # True, wenn kein Laborwert eingetragen wurde
            condition = is_activity & is_empty # True, nur wenn beide True sind
    
            drop_mask = drop_mask | condition # Boolean-Werte werden der drop_mask hinzugefügt

        problematic_groups = filtered_log.loc[drop_mask, 'org:group'].value_counts() # Zeigt an, welche Organisationsgruppe keine Laborwerte eingetragen hat
        print(f'Häufigkeit der Gruppen in fehlerhaften Zeilen: {problematic_groups}')
    
        filtered_log.drop(index=filtered_log[drop_mask].index, inplace=True)# Alle Zeilen mit Boolean-Wert true werden herausgefiltert

    if col_filter is not None:
        filtered_log.drop(col_filter, axis=1, inplace=True)# Alle Spalten mit Namen in spalten_filter werden herausgefiltert
    
    filtered_log = pm4py.filter_start_activities(filtered_log, selected_activities) # Filtern nach Startaktivitäten
    filtered_log = pm4py.filter_end_activities(filtered_log, selected_end_acts)# Filtern nach Endaktivitäten

    if delete_activities is not None:
        filtered_log = pm4py.filter_event_attribute_values(filtered_log, 'concept:name', delete_activities, level='event', retain=False) # filtern nach Aktivitäten

    return filtered_log

# Filterkriterien definieren
criteria_end_sepsis = ['IV Antibiotics', 'ER Sepsis Triage', 'Leucocytes', 'IV Liquid', 'CRP', 'LacticAcid', 'Admission NC', 'ER Triage'] # Endaktivitäten, die herausgefiltert werden sollen
deleted_activities_sepsis = ['Return ER'] # Aktivitäten die generell gelöscht werden sollen im Log
activity_check_sepsis = ['Leucocytes', 'LacticAcid', 'CRP'] # Aktivitäten die auf ihre Vollständigkeit im Log überprüft werden sollen
column_filter_sepsis = ['lifecycle:transition'] # alle complete

filtered_log_sepsis = filter_log(start_act_sepsis, end_act_sepsis, log_sepsis, cases_no_sepsis, end_crit = criteria_end_sepsis, 
                                 delete_activities = deleted_activities_sepsis, check_value_activities=activity_check_sepsis, 
                                 col_filter=column_filter_sepsis) # Übergabe der Parameter an die Filter-Methode
sum_up_log(filtered_log_sepsis)
get_cases_events(filtered_log_sepsis)



column_filter_iacs = ['lifecycle:transition', 'case:program-id', 'case:penalty_BGKV', 'case:penalty_BGP', 'case:penalty_AVUVP', 
                      'case:greening', 'case:basic payment', 'case:penalty_B5F', 'case:penalty_JLP7', 
                      'case:penalty_JLP5'] # complete, 215, False, False, False, True, True, False, False, False
criteria_end_iacs = []

filtered_log_iacs = filter_log(start_act_iacs, end_act_iacs, log_iacs, cases_no_iacs, end_crit=criteria_end_iacs, col_filter=column_filter_iacs)
sum_up_log(filtered_log_iacs)
get_cases_events(filtered_log_iacs)


# DFG aus gefiltertem Log erstellen
dfg_sepsis_filtered = create_dfg_from_log(filtered_log_sepsis)

dfg_iacs_filtered = create_dfg_from_log(filtered_log_iacs)


# Alpha Miner anwenden
def run_alpha_miner(log):
    net, initial_marking, final_marking = pm4py.discover_petri_net_alpha(log)
    pm4py.view_petri_net(net, initial_marking, final_marking)
    return net, initial_marking, final_marking # Kann evtl. weg

net_sepsis, initial_marking_sepsis, final_marking_sepsis = run_alpha_miner(filtered_log_sepsis) # Variablen evtl. umbenennen (auch bei TBR)

alpha_net_iacs, initial_marking_alpha_iacs, final_marking_alpha_iacs = run_alpha_miner(filtered_log_iacs)

# Heuristic Miner anwenden
def run_heuristic_miner(log):
    heuristic_net = pm4py.discover_heuristics_net(log, dependency_threshold=0.99)
    pm4py.view_heuristics_net(heuristic_net)


run_heuristic_miner(filtered_log_sepsis)

# run_heuristic_miner(filtered_log_iacs, dependency_threshold=0.0) # Anpassen

# Inductive Miner anwenden
def run_inductive_miner(log, noise_threshold=0.0):
    inductive_net, initial_marking_inductive, final_marking_inductive = pm4py.discover_petri_net_inductive(log, noise_threshold)
    pm4py.view_petri_net(inductive_net, initial_marking_inductive, final_marking_inductive)
    return inductive_net, initial_marking_inductive, final_marking_inductive

ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis = run_inductive_miner(filtered_log_sepsis)

# ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs = run_inductive_miner(filtered_log_iacs, noise_threshold=0.0) # Anpassen

# Varianten auflisten
def show_filter_variants(log):
    variants = pm4py.get_variants(log)
    with open('variants_output.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in variants.items():
            writer.writerow([key, value])
    #print(f'Varianten: {variants}') # Ausdruck zu lang


show_filter_variants(filtered_log_sepsis)

# show_filter_vars(filtered_log_iacs)

# Filtern und Ausgeben eines Logs mit den 5 häufigsten Varianten sowie Visualisierung
def filter_by_variants(log, k):
    filtered_log_var = pm4py.filter_variants_top_k(log, k)
    return filtered_log_var

for i in [1, 5]:
    filtered_log_var_sepsis = filter_by_variants(filtered_log_sepsis, i)
    run_inductive_miner(filtered_log_var_sepsis)
    create_dfg_from_log(filtered_log_var_sepsis)


# filtered_log_var_iacs = filter_by_vars(filtered_log_iacs, 5)
# run_inductive_miner(filtered_log_var_iacs)
# create_dfg_from_log(filtered_log_var_iacs)

# Temporal Profile erstellen

def get_temporal_profile(log):
    temporal_profile = temporal_profile_discovery.apply(log)
    with open('temporal_profile.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in temporal_profile.items():
            writer.writerow([key, value])

get_temporal_profile(filtered_log_sepsis)

# get_temporal_profile(filtered_log_iacs)

# Umwandeln in Dataframe und Datetime, falls erforderlich
def ensure_df(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    return pm4py.convert_to_dataframe(obj)


def parse_ts(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors='coerce', utc=True).dt.tz_convert(None)


# Mittlere Durchlaufzeit pro Fall berechnen
def get_mean_throughput(log):
    df = ensure_df(log)
    df['time:timestamp'] = parse_ts(df['time:timestamp']) # Zeitstempel in Datetime umwandeln, sofern erforderlich
    case_durations = df.groupby('case:concept:name').agg({'time:timestamp': ['min', 'max']}) # Gruppieren nach Case ID und Berechnen der Dauer
    durations = (case_durations['time:timestamp']['max'] - case_durations['time:timestamp']['min']).dt.total_seconds() / 3600 # Fallweise Dauer in Stunden pro Case berechnen
    
    mean_throughput = durations.mean()
    print(f'Die mittlere Durchlaufzeit pro Fall beträgt: {mean_throughput} Stunden')
    return mean_throughput


mean_throughput_sepsis = get_mean_throughput(filtered_log_sepsis)

# mean_throughput_iacs = get_mean_throughput(filtered_log_iacs)

# Performance-DFG erzeugen
def get_performance_dfg(log):
    perf_dfg, start_activities, end_activities = pm4py.discover_performance_dfg(log)
    pm4py.view_performance_dfg(perf_dfg, start_activities, end_activities)


get_performance_dfg(filtered_log_sepsis)

# get_performance_dfg(filtered_log_iacs)

# Performance-Informationen zu Netz aus Inductive Miner hinzufügen
def get_performance_net(log, net, initial_marking, final_marking):
    parameters = {pn_visualizer.Variants.PERFORMANCE.value.Parameters.FORMAT: "png"}
    gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters, variant=pn_visualizer.Variants.PERFORMANCE, log=log)
    pn_visualizer.view(gviz)


get_performance_net(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# get_performance_net(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)








# Helfer TBR
def tbr_list_to_dataframe(replayed_traces: List[Dict]) -> pd.DataFrame: #(Variablen)namen für TBR ändern
    rows = []
    for i, d in enumerate(replayed_traces):
        row = {
            'trace_index': i,
            'trace_is_fit': d.get('trace_is_fit', None),
            'produced_tokens': d.get('produced_tokens', 0),
            'remaining_tokens': d.get('remaining_tokens', 0),
            'missing_tokens': d.get('missing_tokens', 0),
            'consumed_tokens': d.get('consumed_tokens', 0)
        }
        p = row['produced_tokens']
        r = row['remaining_tokens']
        m = row['missing_tokens']
        c = row['consumed_tokens']
        term1 = (1 - (r / p)) if p else 0.0
        term2 = (1 - (m / c)) if c else 0.0
        row['tbr_fitness'] = 0.5 * term1 + 0.5 * term2
        rows.append(row)
    return pd.DataFrame(rows)


# TBR mit Modell aus Inductive Miner
def tbr_ind(log, net, initial_marking, final_marking):
    replayed_traces = pm4py.conformance_diagnostics_token_based_replay(log, net, initial_marking, final_marking)
    df_tbr = tbr_list_to_dataframe(replayed_traces)
    print(f'Head des Token Based Replay mit Inductive Miner: \n{df_tbr.head(10)}')
    return df_tbr


df_tbr_sepsis = tbr_ind(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# df_tbr_iacs = tbr_ind(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Visualisierungen (Histogramm der Fitness und Balkenplot fit/unfit)
def plot_tbr_fitness_hist(df_tbr: pd.DataFrame, bins: int = 20):
    plt.figure()
    plt.hist(df_tbr['tbr_fitness'].dropna(), bins=bins)
    plt.xlabel('TBR-Fitness (pro Trace)')
    plt.ylabel('Anzahl')
    plt.title('Histogramm: TBR-Fitness')
    plt.show()


def plot_tbr_fit_flag(df_tbr: pd.DataFrame):
    counts = df_tbr['trace_is_fit'].value_counts(dropna=False)
    labels = [str(k) for k in counts.index]
    values = counts.values
    plt.figure()
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels)
    plt.xlabel('trace_is_fit')
    plt.ylabel('Anzahl Traces')
    plt.title('Fit vs. Unfit (TBR)')
    plt.show()


plot_tbr_fitness_hist(df_tbr_sepsis)
plot_tbr_fit_flag(df_tbr_sepsis)

# plot_tbr_fitness_hist(df_tbr_iacs)
# plot_tbr_fit_flag(df_tbr_iacs)







# Throughput Analysis: TBR mit angepassten Einstellungen
def tbr_throughput(log, net, initial_marking, final_marking):
    parameters_tbr = {
        token_based_replay.Variants.TOKEN_REPLAY.value.Parameters.DISABLE_VARIANTS: True,
        token_based_replay.Variants.TOKEN_REPLAY.value.Parameters.ENABLE_PLTR_FITNESS: True
    }
    replayed_traces, place_fitness, trans_fitness, unwanted_activities = token_based_replay.apply(
        log, net, initial_marking, final_marking, parameters=parameters_tbr
    )

    df_tbr = tbr_list_to_dataframe(replayed_traces)
    print(f'Head des Token Based Replay mit Throughput Analysis: \n{df_tbr.head(10)}')
    print(f'Place Fitness {place_fitness}')
    print(f'Trans Fitness {trans_fitness}')
    print(f'Unwanted Activities {unwanted_activities}')

    return replayed_traces, place_fitness, trans_fitness, unwanted_activities


log_diagnostics_sepsis = pm4py.convert_to_event_log(filtered_log_sepsis)
replayed_traces_sepsis, place_fitness_sepsis, trans_fitness_sepsis, unwanted_activities_sepsis = tbr_throughput(log_diagnostics_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# log_diagnostics_iacs = pm4py.convert_to_event_log(filtered_log_iacs)
# replayed_traces_iacs, place_fitness_iacs, trans_fitness_iacs, unwanted_activities_iacs = tbr_throughput(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Throughput Analysis der falsch ausgeführten Transitionen und Ausgabe
def througput_trans(log, trans_fitness):
    trans_diagnostics = duration_diagnostics.diagnose_from_trans_fitness(log, trans_fitness)
    for trans in trans_diagnostics:
        print(f'Throughput Analysis der falsch ausgeführten Transitionen{trans, trans_diagnostics[trans]}')


througput_trans(log_diagnostics_sepsis, trans_fitness_sepsis)

# througput_trans(filtered_log_iacs, trans_fitness_iacs)

# Throughput Analysis der nicht im Modell enthaltenen Aktivitäten und Ausgabe
def throughput_act(log, unwanted_activities):
    act_diagnostics = duration_diagnostics.diagnose_from_notexisting_activities(log, unwanted_activities)
    for act in act_diagnostics:
        print(act, act_diagnostics[act])


throughput_act(log_diagnostics_sepsis, unwanted_activities_sepsis)

# throughput_act(filtered_log_iacs, unwanted_activities_iacs)

# Vorbereitung Root Cause Analysis
string_attributes = ['org:group']
numeric_attributes = []
parameters = {'string_attributes': string_attributes, 'numeric_attributes': numeric_attributes}

# Root Cause Analysis der falsch ausgeführten Transitionen
def rca_trans(log, trans_fitness):
    trans_root_cause = root_cause_analysis.diagnose_from_trans_fitness(log, trans_fitness, parameters=parameters)
    for trans in trans_root_cause:
        clf = trans_root_cause[trans]['clf']
        feature_names = trans_root_cause[trans]['feature_names']
        classes = trans_root_cause[trans]['classes']
        # Visualization can be called
        gviz = dt_vis.apply(clf, feature_names, classes)
        dt_vis.view(gviz)


rca_trans(log_diagnostics_sepsis, trans_fitness_sepsis)

# rca_trans(filtered_log_iacs, trans_fitness_iacs)

# RCA der ausgeführten Aktivitäten, die nicht im Prozessmodell enthalten sind, und Ausgabe
def rca_act(log, unwanted_activities):
    act_root_cause = root_cause_analysis.diagnose_from_notexisting_activities(log, unwanted_activities, parameters=parameters)
    for act in act_root_cause:
        clf = act_root_cause[act]["clf"]
        feature_names = act_root_cause[act]["feature_names"]
        classes = act_root_cause[act]["classes"]
        # Visualization can be called
        gviz = dt_vis.apply(clf, feature_names, classes)
        dt_vis.view(gviz)


rca_act(log_diagnostics_sepsis, unwanted_activities_sepsis)

# rca_act(filtered_log_iacs, unwanted_activities_iacs)






# Alignments bestimmen mit Modell aus Inductive Miner
def alignments_inductive(log, net, initial_marking, final_marking):
    aligned_traces_ind = pm4py.conformance_diagnostics_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Alignments mit Inductive Miner: {aligned_traces_ind}') # Ausdruck zu lang
    return aligned_traces_ind

aligned_traces_ind_sepsis = alignments_inductive(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# aligned_traces_ind_iacs = alignments_inductive(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Visualisierung der Alignments (Histogramm der Alignment-Kosten und Move-Typen)
def extract_alignment_stats(aligned_traces) -> pd.DataFrame: # Variablennamen etc. evtl. ändern
    rows = []
    for i, d in enumerate(aligned_traces):
        row = {
            'trace_index': i,
            'cost': d.get('cost', None),
            'fitness': d.get('fitness', None)
        }
        # Moves zählen
        align_pairs = d.get('alignment', [])
        sync = 0
        move_log = 0
        move_model = 0
        for a, b in align_pairs:
            if a != '>>' and b != '>>':
                sync += 1
            elif a != '>>' and b == '>>':
                move_log += 1
            elif a == '>>' and b != '>>':
                move_model += 1
        row['sync_moves'] = sync
        row['move_on_log'] = move_log
        row['move_on_model'] = move_model
        rows.append(row)
    return pd.DataFrame(rows)

def plot_alignment_cost_hist(df_stats: pd.DataFrame, bins: int = 20):
    plt.figure()
    plt.hist(df_stats['cost'].dropna(), bins=bins)
    plt.xlabel('Alignment-Kosten')
    plt.ylabel('Anzahl')
    plt.title('Histogramm: Alignment-Kosten')
    plt.show()

def plot_move_type_bars(df_stats: pd.DataFrame):
    agg = df_stats[['sync_moves', 'move_on_log', 'move_on_model']].sum()
    labels = list(agg.index)
    values = list(agg.values)
    plt.figure()
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels, rotation=0)
    plt.ylabel('Gesamtanzahl Moves (über alle Traces)')
    plt.title('Move-Typen (aggregiert)')
    plt.show()

df_stats_sepsis = extract_alignment_stats(aligned_traces_ind_sepsis)
plot_alignment_cost_hist(df_stats_sepsis)
plot_move_type_bars(df_stats_sepsis)

# df_stats_iacs = extract_alignment_stats(aligned_traces_ind_iacs)
# plot_alignment_cost_hist(df_stats_iacs)
# plot_move_type_bars(df_stats_iacs)











# Fitness zwischen Log und Modell berechnen (TBR und Alignments)
def get_replay_fitness_tbr(log, net, initial_marking, final_marking):
    rp_fitness_tbr = pm4py.fitness_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Token-Based Replay: {rp_fitness_tbr}')
    return rp_fitness_tbr['log_fitness']


def get_replay_fitness_align(log, net, initial_marking, final_marking):
    rp_fitness_align = pm4py.fitness_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Alignments: {rp_fitness_align}')
    return rp_fitness_align['log_fitness']


fitness_tbr_sepsis = get_replay_fitness_tbr(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)
fitness_align_sepsis = get_replay_fitness_align(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# fitness_tbr_iacs = get_replay_fitness_tbr(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)
# fitness_align_iacs get_replay_fitness_align(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# Precision zwischen Log und Modell berechnen (ETConformance (TBR) und Align-ETConformance(Alignments))
def get_precision_tbr(log, net, initial_marking, final_marking):
    precision_tbr = pm4py.precision_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Token-Based Replay beträgt: {precision_tbr}') # evtl. nur bestimmten Wert ausgeben
    return precision_tbr


def get_precision_align(log, net, initial_marking, final_marking):
    precision_align = pm4py.precision_alignments(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Alignments beträgt: {precision_align}')
    return precision_align


precision_tbr_sepsis = get_precision_tbr(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)
# precision_align_sepsis = get_precision_align(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# precision_tbr_iacs = get_precision_tbr(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)
# precision_align_iacs = get_precision_align(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)

# F1-Score zwischen Log und Modell berechnen, Fitness und Precision gegeben
def get_f1_score(fitness, precision):
    if fitness + precision == 0: # Division durch 0 vermeiden
        return 0.0
    f1_score = 2 * (fitness * precision) / (fitness + precision)
    print(f'Der F1-Score beträgt: {f1_score}')
    return f1_score


f1_score_tbr_sepsis = get_f1_score(fitness_tbr_sepsis, precision_tbr_sepsis)
# f1_score_align_sepsis = get_f1_score(fitness_align_sepsis, precision_align_sepsis)

# f1_score_tbr_iacs = get_f1_score(fitness_tbr_iacs, precision_tbr_iacs)
# f1_score_align_iacs = get_f1_score(fitness_align_iacs, precision_align_iacs)

# Generalization und Simplicity zwischen Log und Modell berechnen
def get_generalization(log, net, initial_marking, final_marking):
    generalization = generalization_evaluator.apply(log, net, initial_marking, final_marking)
    print(f'Die Generalization beträgt: {generalization}')


def get_simplicity(log, net, initial_marking, final_marking):
    simplicity = simplicity_evaluator.apply(net)
    print(f'Die Simplicity beträgt: {simplicity}')


get_generalization(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)
get_simplicity(filtered_log_sepsis, ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis)

# get_generalization(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)
# get_simplicity(filtered_log_iacs, ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs)








# Erstellen eines Histogramms der Startaktivitäten
def plot_start_activities(log):
    df = ensure_df(log)
    first_acts = df.groupby('case:concept:name')['concept:name'].first() # Erste Aktivität pro Fall ermitteln
    counts = first_acts.value_counts() # Zählen der Startaktivitäten
    
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(counts)), counts.values)
    plt.xticks(range(len(counts)), counts.index, rotation=45, ha='right')
    plt.xlabel('Startaktivität')
    plt.ylabel('Anzahl der Fälle')
    plt.title('Startaktivitäten')
    plt.tight_layout()
    plt.show()


plot_start_activities(filtered_log_sepsis)

# plot_start_activities(filtered_log_iacs)
    
# Erstellen eines Histogramms der Endaktivitäten
def plot_end_activities(log):
    df = ensure_df(log)
    last_acts = df.groupby('case:concept:name')['concept:name'].last() # Letzte Aktivität pro Fall ermitteln
    counts = last_acts.value_counts() # Zählen der Endaktivitäten
    
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(counts)), counts.values)
    plt.xticks(range(len(counts)), counts.index, rotation=45, ha='right')
    plt.xlabel('Endaktivität')
    plt.ylabel('Anzahl der Fälle')
    plt.title('Endaktivitäten')
    plt.tight_layout()
    plt.show()


plot_end_activities(filtered_log_sepsis)

# plot_end_activities(filtered_log_iacs)







def plot_activity_frequencies(obj):
    df = ensure_df(obj)
    counts = df['concept:name'].value_counts()
    labels = counts.index.tolist()
    values = counts.values.tolist()
    plt.figure()
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.ylabel('Häufigkeit (Events)')
    plt.title('Aktivitätshäufigkeiten')
    plt.tight_layout()
    plt.show()

plot_activity_frequencies(filtered_log_sepsis)

# plot_activity_frequencies(filtered_log_iacs)


def plot_case_duration_hist(obj, bins: int = 30): # Durchlaufzeit der Fälle
    df = ensure_df(obj)
    cid = 'case:concept:name'
    ts  = 'time:timestamp'
    df[ts] = parse_ts(df[ts])
    ag = df.groupby(cid)[ts].agg(['min', 'max'])
    dur_days = (ag['max'] - ag['min']).dt.total_seconds() / 86400.0
    plt.figure()
    plt.hist(dur_days.dropna(), bins=bins)
    plt.xlabel('Case-Dauer [Tage]')
    plt.ylabel('Anzahl Cases')
    plt.title('Histogramm der Case-Dauern')
    plt.show()

plot_case_duration_hist(filtered_log_sepsis)

# plot_case_duration_hist(filtered_log_iacs)


def plot_events_per_day(obj):
    df = ensure_df(obj)
    ts = 'time:timestamp'
    df[ts] = parse_ts(df[ts])
    daily = df.set_index(ts).resample('D').size()
    plt.figure()
    daily.plot()
    plt.xlabel('Datum')
    plt.ylabel('Ereignisse pro Tag')
    plt.title('Event-Zeitreihe (täglich)')
    plt.tight_layout()
    plt.show()

plot_events_per_day(filtered_log_sepsis)

# plot_events_per_day(filtered_log_iacs)





# Spalte bei Sepsis-Datensatz umbenennen für Social Network Analysis
log_sepsis_sna = filtered_log_sepsis.rename(columns = {'org:group' : 'org:resource'})

# Übergabe von Arbeit ermitteln und anzeigen
def get_handover_of_work(log):
    handover_values = pm4py.discover_handover_of_work_network(log)
    pm4py.view_sna(handover_values)


get_handover_of_work(log_sepsis_sna)

# get_handover_of_work(log_iacs_sna)

# Ermitteln und Anzeigen, wie oft Subcontracting vorkommt
def get_subcontracting(log):
    subcont_values = pm4py.discover_subcontracting_network(log)
    pm4py.view_sna(subcont_values)


get_subcontracting(log_sepsis_sna)

# get_subcontracting(log_iacs_sna)

# Ermitteln und Anzeigen, wie oft zusammengearbeitet wird
def get_working_together(log):
    work_together_values = pm4py.discover_working_together_network(log)
    pm4py.view_sna(work_together_values)


get_working_together(log_sepsis_sna)

# get_working_together(log_iacs_sna)

# Ähnlichkeiten der Arbeitsmuster zwischen Individuen ermitteln und anzeigen
def get_similar_activities(log):
    similar_act = pm4py.discover_activity_based_resource_similarity(log) # Benennung Var. wg Cluster-A.
    pm4py.view_sna(similar_act)
    return similar_act


similar_activities_sepsis = get_similar_activities(log_sepsis_sna)

# similar_activities_iacs = get_similar_activities(log_iacs_sna)

# Orginisationale Rollen entdecken und ausgeben
def get_orga_roles(log):
    roles = pm4py.discover_organizational_roles(log)
    print([x[0] for x in roles])


get_orga_roles(log_sepsis_sna)

# get_orga_roles(log_iacs_sna)

# Cluster-Analyse nach ähnlichen Aktivitäten der Individuen
'''from pm4py.algo.organizational_mining.sna import util # wahrscheinlich nicht gut

def cluster_similar_act(log, sim_act): # Evtl. noch nach anderen Sachen clustern
    clustering = util.cluster_affinity_propagation(sim_act)


cluster_similar_act(log_sepsis_sna, similar_activities_sepsis_sna)

# cluster_similar_act(log_iacs_sna, similar_activities_iacs_sna)''' 


