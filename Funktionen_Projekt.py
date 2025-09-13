# Libraries importieren
import pm4py
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn # Cluster-Analyse SNA funktioniert sonst nicht
import csv
import string
import re
from typing import List, Dict
from pm4py.algo.discovery.temporal_profile import algorithm as temporal_profile_discovery
from pm4py.visualization.petri_net import visualizer as pn_visualizer # Für Zusatzinfos Performance
from pm4py.algo.conformance.tokenreplay import algorithm as token_based_replay
from pm4py.algo.conformance.tokenreplay.diagnostics import duration_diagnostics
from pm4py.algo.conformance.tokenreplay.diagnostics import root_cause_analysis
from pm4py.visualization.decisiontree import visualizer as dt_vis
from pm4py.algo.evaluation.generalization import algorithm as generalization_evaluator
from pm4py.algo.evaluation.simplicity import algorithm as simplicity_evaluator
from pm4py.util import constants
from pm4py.statistics.traces.generic.log import case_statistics
from pm4py.visualization.graphs import visualizer as graphs_visualizer
from pm4py.algo.organizational_mining.sna import util
from string import ascii_uppercase
from itertools import product



# Funktion zum Import von .xes-Dateien
def import_xes(path_to_xes: str):
    log = pm4py.read_xes(path_to_xes)
    return log

# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
def sum_up_log (log):
    print(f'Beschreibung des Datensatzes: \n{log.describe()}')
    print(f'Head des Datensatzes: \n{log.head(10)}\nTail des Datensatzes: \n{log.tail(10)}') # Jeweils Ausgabe von 10 Zeilen

# Anzahl der Cases und Events anzeigen
def get_cases_events(log):
    cases_no = len(log['case:concept:name'].unique()) # Zählen der Cases
    events_no = len(log) # Zählen der Events
    
    print(f'Anzahl Cases: {cases_no}\nAnzahl Events: {events_no}')
    return cases_no

# Start- und Endaktivitäten anzeigen
def get_start_end_act(log):
    start_act = pm4py.get_start_activities(log)
    end_act = pm4py.get_end_activities(log)
    
    print(f'Startaktivitäten: {start_act}\nEndaktivitäten: {end_act}')
    return start_act, end_act

# DFG aus Log erstellen und anzeigen
def create_dfg_from_log(log_input):
    dfg, start_activities, end_activities = pm4py.discover_dfg(log_input)
    dfg_output = pm4py.view_dfg(dfg, start_activities, end_activities)
    return dfg_output

# Filteralgorithmus
def filter_log(start_acts, end_acts, log, no_of_cases, min_ratio=0.1, end_crit = None, delete_activities = None, check_value_activities=None, col_filter=None):
    selected_activities = set()
    selected_end_acts = set()
    filtered_log = log

    for activity, count in start_acts.items(): # Iteration im Dict über den Namen der Startaktivität und ihre Häufigkeit.
        if count / no_of_cases >= min_ratio: # Nur Startaktivitäten mit einem Mindestanteil von 10% an Gesamt-Cases werden berücksichtigt.
            selected_activities.add(activity) # Sammeln der Aktivitäten in einer Liste
    
    if end_crit is not None:
        selected_end_acts = set(end_acts.keys()) - set(end_crit) # Löschen von Cases, die nicht ordungsgemäß mit Release oder Return ER geendet haben.
    
    
    if check_value_activities is not None:

        drop_mask = pd.Series(False, index=filtered_log.index) # Erstellen einer Maske mit Boolean-Werte

        for item in check_value_activities: # Sepsisspezifisch!! Es wird überprüft, ob die Laborwerte tatsächlich vorhanden sind, oder ob nur die Aktivität ohne Laborwert ausgeführt wurde

            is_activity = filtered_log['concept:name'] == item # True, wenn die Aktivität ausgeführt wrude
            is_empty = filtered_log[item].isna() # True, wenn kein Laborwert eingetragen wurde
            condition = is_activity & is_empty # True, nur wenn beide True sind
    
            drop_mask = drop_mask | condition # Boolean-Werte werden der drop_mask hinzugefügt

        problematic_groups = filtered_log.loc[drop_mask, 'org:group'].value_counts() # Zeigt an, welche Organisationsgruppe keine Laborwerte eingetragen hat
        print(f'Häufigkeit der Gruppen in fehlerhaften Zeilen: {problematic_groups}')
    
        filtered_log.drop(index=filtered_log[drop_mask].index, inplace=True) # Alle Zeilen mit Boolean-Wert true werden herausgefiltert

    if col_filter is not None:
        filtered_log.drop(col_filter, axis=1, inplace=True) # Alle Spalten mit Namen in spalten_filter werden herausgefiltert
    
    filtered_log = pm4py.filter_start_activities(filtered_log, selected_activities) # Filtern nach Startaktivitäten
    filtered_log = pm4py.filter_end_activities(filtered_log, selected_end_acts) # Filtern nach Endaktivitäten

    if delete_activities is not None:
        filtered_log = pm4py.filter_event_attribute_values(filtered_log, 'concept:name', delete_activities, level='event', retain=False) # Filtern nach Aktivitäten

    return filtered_log

# Umwandeln in Dataframe, falls erforderlich
def ensure_df(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    return pm4py.convert_to_dataframe(obj) # Vor allem für Event Log zu Dataframe

# Umwandeln von Zeitstempeln in Datetime, falls erforderlich
def parse_ts(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors='coerce', utc=True).dt.tz_convert(None) #

# Spaltenweises Zählen der Häufigkeit der Werte
def get_column_count(log, columns: list): # Übergabe muss eine Liste sein
    df = ensure_df(log)
    
    for col in columns: # Durchzählen und Ausgabe für alle übergebenen Spalten
        if col in df.columns:
            value_counts = df[col].value_counts()
            print(f'\nAnzahlen für Spalte {value_counts}')
        else:
            print(f'\nSpalte {col} existiert nicht') # Fehlermeldung bei nicht existierender Spalte

# Erstellen eines Histogramms der Startaktivitäten
def plot_start_activities(log):
    df = ensure_df(log)
    first_acts = df.groupby('case:concept:name')['concept:name'].first() # Erste Aktivität pro Case ermitteln
    counts = first_acts.value_counts() # Zählen der Startaktivitäten
    
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(counts)), counts.values)
    plt.xticks(range(len(counts)), counts.index, rotation=45, ha='right')
    plt.xlabel('Startaktivität')
    plt.ylabel('Anzahl der Cases')
    plt.title('Startaktivitäten')
    plt.tight_layout()
    plt.show()

# Erstellen eines Histogramms der Endaktivitäten
def plot_end_activities(log):
    df = ensure_df(log)
    last_acts = df.groupby('case:concept:name')['concept:name'].last() # Letzte Aktivität pro Case ermitteln
    counts = last_acts.value_counts() # Zählen der Endaktivitäten
    
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(counts)), counts.values)
    plt.xticks(range(len(counts)), counts.index, rotation=45, ha='right')
    plt.xlabel('Endaktivität')
    plt.ylabel('Anzahl der Cases')
    plt.title('Endaktivitäten')
    plt.tight_layout()
    plt.show()

# Erstellen eines Diagramms der Aktivitätshäufigkeiten
def plot_activity_frequencies(obj):
    df = ensure_df(obj)
    counts = df['concept:name'].value_counts() # Durchzählen der Aktivitäten, Series
    labels = counts.index.tolist() # Indizes in Liste umwandeln
    values = counts.values.tolist() # Werte in Liste umwandeln
    plt.figure()
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.ylabel('Häufigkeit (Events)')
    plt.title('Aktivitätshäufigkeiten')
    plt.tight_layout()
    plt.show()

# Erstellen eines Histogramms der Durchlaufzeiten der Cases
def plot_case_duration_hist(obj, bins: int = 30):
    df = ensure_df(obj)
    cid = 'case:concept:name'
    ts  = 'time:timestamp'
    df[ts] = parse_ts(df[ts])
    ag = df.groupby(cid)[ts].agg(['min', 'max']) # Gruppieren nach Case-ID/minimalem und maximalem Zeitstempel
    dur_days = (ag['max'] - ag['min']).dt.total_seconds() / 86400.0 # Berechnung der Dauer in Tagen
    plt.figure()
    plt.hist(dur_days.dropna(), bins=bins)
    plt.xlabel('Case-Dauer [Tage]')
    plt.ylabel('Anzahl Cases')
    plt.title('Histogramm der Case-Dauern')
    plt.show()

# Erstellen eines Diagramms der Events pro Tag
def plot_events_per_day(obj):
    df = ensure_df(obj)
    ts = 'time:timestamp'
    df[ts] = parse_ts(df[ts])
    daily = df.set_index(ts).resample('D').size() # Series, enthält Anzahl der Events pro Tag
    plt.figure()
    daily.plot()
    plt.xlabel('Datum')
    plt.ylabel('Events pro Tag')
    plt.title('Event-Zeitreihe (täglich)')
    plt.tight_layout()
    plt.show()

# # Alpha Miner anwenden, Modell erstellen und anzeigen
def run_alpha_miner(log):
    net, initial_marking, final_marking = pm4py.discover_petri_net_alpha(log)
    pm4py.view_petri_net(net, initial_marking, final_marking)
    return net, initial_marking, final_marking

# Heuristic Miner anwenden, Modell anzeigen und Umwandlung in Petri-Netz
def run_heuristic_miner(log, dependency_threshold=0.99):
    heuristic_net = pm4py.discover_heuristics_net(log, dependency_threshold)
    pm4py.view_heuristics_net(heuristic_net)
    heuristic_to_petri = pm4py.objects.conversion.heuristics_net.variants.to_petri_net.apply(heuristic_net, parameters=None) # Umwandeln in Petri-Netz
    return heuristic_to_petri

# Inductive Miner für Sepsis-Datensatz anwenden, Modell erstellen und anzeigen
def run_inductive_miner_sepsis(log): # Programm funktionert nicht bei Übergabe von noise_threshold, daher getrennte Funktion
    inductive_net, initial_marking_inductive, final_marking_inductive = pm4py.discover_petri_net_inductive(log, noise_threshold=0.1)
    pm4py.view_petri_net(inductive_net, initial_marking_inductive, final_marking_inductive)
    return inductive_net, initial_marking_inductive, final_marking_inductive

# Inductive Miner für IACS-Datensatz anwenden, Modell erstellen und anzeigen
def run_inductive_miner_iacs(log):
    inductive_net, initial_marking_inductive, final_marking_inductive = pm4py.discover_petri_net_inductive(log, noise_threshold=0.5)
    pm4py.view_petri_net(inductive_net, initial_marking_inductive, final_marking_inductive)
    return inductive_net, initial_marking_inductive, final_marking_inductive

# Varianten auflisten
def show_filter_variants(log):
    variants = pm4py.get_variants(log)
    with open('variants_output.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in variants.items():
            writer.writerow([key, value])
    length_dict = {}
    for variant, count in variants.items():
        l =len(variant)
        if l in length_dict:
            length_dict[l] += count
        else:
            length_dict[l] = count
    
    total_variants = sum(length_dict.values())
    avg_length = sum(length *count for length, count in length_dict.items()) / total_variants
    min_length = min(length_dict)
    max_length = max(length_dict)
    mode_length = max(length_dict, key=length_dict.get)

    print('Längen-Häufigkeiten', dict(sorted(length_dict.items())))
    print(f'Durchschnittliche Länge: {avg_length}')
    print(f'Kürzeste Variante: {min_length}')
    print(f'Längste Variante: {max_length}')
    print(f'Häufigste Variante: {mode_length} (Häufigkeit: {length_dict[mode_length]})')

# Filtern und Ausgeben eines Logs mit den k häufigsten Varianten
def filter_by_variants(log, k):
    filtered_log_var = pm4py.filter_variants_top_k(log, k)
    return filtered_log_var

# Mittlere Durchlaufzeit pro Case berechnen
def get_mean_throughput(log):
    durations = pm4py.get_all_case_durations(log) # Liste aller Case-Dauern
    mean_throughput = np.mean(durations) / 3600 # Umrechnung in Stunden
    print(f'Die mittlere Durchlaufzeit pro Case beträgt: {mean_throughput} Stunden')
    return mean_throughput

# Performance-DFG erzeugen
def get_performance_dfg(log):
    perf_dfg, start_activities, end_activities = pm4py.discover_performance_dfg(log)
    pm4py.view_performance_dfg(perf_dfg, start_activities, end_activities)

# Temporal Profile erstellen
def get_temporal_profile(log):
    temporal_profile = temporal_profile_discovery.apply(log)
    with open('temporal_profile.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in temporal_profile.items():
            writer.writerow([key, value])

# Neues Log für Temporal Profile erstellen, in dem "begin" und "finish" zu einer Aktivität zusammengefasst werden
def merge_begin_finish(log, begin, finish, new_name, allow_abort=False):
    
    mask_begin = log['concept:name'] == begin
    mask_finish = log['concept:name'] == finish
    if allow_abort: # Speziell für die Endaktivität 'abort payment' im IACS-Datensatz
        mask_finish = mask_finish | (log['concept:name'] == 'abort payment')
    log_begin = log[mask_begin].copy() # Kopie der Zeilen mit der Startaktivität
    log_finish = log[mask_finish].copy() # Kopie der Zeilen mit der Endaktivität

    log_begin = log_begin.rename(columns={'time:timestamp': 'start_time'})
    log_finish = log_finish.rename(columns={'time:timestamp': 'end_time'})

    # Sortieren und Zusammenführen der DataFrames basierend auf Case-ID und Zeitstempeln
    merged = pd.merge_asof(
        log_begin.sort_values('start_time'),
        log_finish.sort_values('end_time'),
        by='case:identity:id',
        left_on='start_time',
        right_on='end_time',
        direction='forward'
    )

    result = log.copy()
    indices_to_drop = []

    # Iteration über die zusammengeführten Zeilen und Aktualisierung des Original-DataFrames
    for _, row in merged.iterrows():
        
        begin_index = result[
            (result['concept:name'] == begin) &
            (result['case:identity:id'] == row['case:identity:id']) &
            (result['time:timestamp'] == row['start_time'])
        ].index

        finish_index = result[
            (result['concept:name'].isin([finish] + (['abort payment'] if allow_abort else []))) &
            (result['case:identity:id'] == row['case:identity:id']) &
            (result['time:timestamp'] == row['end_time'])
        ].index

        if len(begin_index) > 0:
            i = begin_index[0]
        result.at[i, 'concept:name'] = new_name # Neue Aktivität
        result.at[i, 'start_timestamp'] = row['start_time'] # Neue Spalte für Startzeit
        result.at[i, 'time:timestamp'] = row['end_time']
        indices_to_drop.extend(finish_index)

    result = result.drop(index=indices_to_drop).reset_index(drop=True)
    return result

# Performance-Informationen zu Netz aus Inductive Miner hinzufügen
def get_performance_net(log, net, initial_marking, final_marking):
    parameters = {pn_visualizer.Variants.PERFORMANCE.value.Parameters.FORMAT: 'png'}
    gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters, variant=pn_visualizer.Variants.PERFORMANCE, log=log)
    pn_visualizer.view(gviz)

# Process Tree mit Inductive Miner erstellen und anzeigen
def get_process_tree_ind(log):
    tree = pm4py.discover_process_tree_inductive(log)
    pm4py.view_process_tree(tree, format='png')

# Helfer TBR erstellen
def tbr_list_to_dataframe(replayed_traces: List[Dict]) -> pd.DataFrame: #(Variablen)namen für TBR ändern
    rows = []
    for i, d in enumerate(replayed_traces): # Erstellen eines Dict mit Trace-Index, Fitness und Anzahlen der Tokens 
        row = {
            'trace_index': i,
            'trace_is_fit': d.get('trace_is_fit', None),
            'produced_tokens': d.get('produced_tokens', 0),
            'remaining_tokens': d.get('remaining_tokens', 0),
            'missing_tokens': d.get('missing_tokens', 0),
            'consumed_tokens': d.get('consumed_tokens', 0)
        }
        p = row['produced_tokens'] # Auslesen der Anzahlen aus dem Dict
        r = row['remaining_tokens']
        m = row['missing_tokens']
        c = row['consumed_tokens']
        term1 = (1 - (r / p)) if p else 0.0
        term2 = (1 - (m / c)) if c else 0.0
        row['tbr_fitness'] = 0.5 * term1 + 0.5 * term2 # Trace-Fitness anhand der Formel berechnen
        rows.append(row)
    return pd.DataFrame(rows)

# TBR durchführen
def tbr_ind(log, net, initial_marking, final_marking):
    replayed_traces = pm4py.conformance_diagnostics_token_based_replay(log, net, initial_marking, final_marking)
    df_tbr = tbr_list_to_dataframe(replayed_traces)
    print(f'Head des Token Based Replay mit Inductive Miner: \n{df_tbr.head(10)}')
    return df_tbr

# Visualisierungen zum TBR erstellen (Histogramm der Fitness und Balkenplot fit/unfit)
def plot_tbr_fitness_hist(df_tbr: pd.DataFrame, bins: int = 20):
    plt.figure()
    plt.hist(df_tbr['tbr_fitness'].dropna(), bins=bins)
    plt.xlabel('TBR-Fitness (pro Trace)')
    plt.ylabel('Anzahl')
    plt.title('Histogramm: TBR-Fitness')
    plt.show()

def plot_tbr_fit_flag(df_tbr: pd.DataFrame):
    counts = df_tbr['trace_is_fit'].value_counts(dropna=False) # Traces, die fit/unfit sind, zählen
    labels = [str(k) for k in counts.index] # List der Indizes als Strings
    values = counts.values # Liste der Werte erstellen
    plt.figure()
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels)
    plt.xlabel('Trace ist fit')
    plt.ylabel('Anzahl Traces')
    plt.title('Fit vs. Unfit (TBR)')
    plt.show()

# TBR mit angepassten Einstellungen zur Vorbereitung der Throughput Analysis durchführen
def tbr_throughput(log, net, initial_marking, final_marking): # Modell mit maximaler Fitness (Inductive Miner) funktionert nicht
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

# Throughput Analysis der falsch ausgeführten Transitionen durchführen und diese ausgeben
def througput_trans(log, trans_fitness):
    trans_diagnostics = duration_diagnostics.diagnose_from_trans_fitness(log, trans_fitness)
    for trans in trans_diagnostics:
        print(f'Throughput Analysis der falsch ausgeführten Transitionen{trans, trans_diagnostics[trans]}')

# Throughput Analysis der nicht im Modell enthaltenen Aktivitäten durchführen und diese ausgeben
def throughput_act(log, unwanted_activities):
    act_diagnostics = duration_diagnostics.diagnose_from_notexisting_activities(log, unwanted_activities)
    for act in act_diagnostics:
        print(act, act_diagnostics[act])

# Vorbereitung Root Cause Analysis durchführen
string_attributes = ['org:group']
numeric_attributes = []
parameters = {'string_attributes': string_attributes, 'numeric_attributes': numeric_attributes}

# Root Cause Analysis der falsch ausgeführten Transitionen durchführen und Decision Trees erstellen
def rca_trans(log, trans_fitness):
    trans_root_cause = root_cause_analysis.diagnose_from_trans_fitness(log, trans_fitness, parameters=parameters)
    for trans in trans_root_cause:
        clf = trans_root_cause[trans]['clf']
        feature_names = trans_root_cause[trans]['feature_names']
        classes = trans_root_cause[trans]['classes']
        # Visualization can be called
        gviz = dt_vis.apply(clf, feature_names, classes)
        dt_vis.view(gviz)

# RCA der ausgeführten Aktivitäten, die nicht im Prozessmodell enthalten sind, durchführen und Decision Trees erstellen
def rca_act(log, unwanted_activities):
    act_root_cause = root_cause_analysis.diagnose_from_notexisting_activities(log, unwanted_activities, parameters=parameters)
    for act in act_root_cause:
        clf = act_root_cause[act]['clf']
        feature_names = act_root_cause[act]['feature_names']
        classes = act_root_cause[act]['classes']
        # Visualization can be called
        gviz = dt_vis.apply(clf, feature_names, classes)
        dt_vis.view(gviz)

# Alignments bestimmen und ausgeben
def alignments_inductive(log, net, initial_marking, final_marking):
    aligned_traces_ind = pm4py.conformance_diagnostics_alignments(log, net, initial_marking, final_marking)
    return aligned_traces_ind

# Visualisierung der Alignments erstellen (Histogramm der Alignment-Kosten und Move-Typen)
def extract_alignment_stats(aligned_traces) -> pd.DataFrame:
    rows = []
    for i, d in enumerate(aligned_traces): # Aufzählen mit Trace-Index, Kosten, Fitness und Moves
        row = {
            'trace_index': i,
            'cost': d.get('cost', None),
            'fitness': d.get('fitness', None)
        }
        align_pairs = d.get('alignment', [])
        sync = 0 # Zähler der Move-Typen
        move_log = 0
        move_model = 0
        for a, b in align_pairs: # Erhöhen der Zähler gemäß beobachtetem Move-Typ
            if a != '>>' and b != '>>':
                sync += 1
            elif a != '>>' and b == '>>':
                move_log += 1
            elif a == '>>' and b != '>>':
                move_model += 1
        row['sync_moves'] = sync # gibt Anzahl synchroner Moves an
        row['move_on_log'] = move_log # gibt Anzahl der Log Moves an
        row['move_on_model'] = move_model # gibt Anzahl der Model Moves an
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
    agg = df_stats[['sync_moves', 'move_on_log', 'move_on_model']].sum() # Move-Typen aggregieren
    labels = list(agg.index) # Liste der Move-Typen erstellen
    values = list(agg.values) # Liste der Werte erstellen
    plt.figure()
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels, rotation=0)
    plt.ylabel('Gesamtanzahl Moves (über alle Traces)')
    plt.title('Move-Typen (aggregiert)')
    plt.show()

# Fitness zwischen Log und Modell berechnen (mit TBR und Alignments)
def get_replay_fitness_tbr(log, net, initial_marking, final_marking):
    rp_fitness_tbr = pm4py.fitness_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Token-Based Replay: {rp_fitness_tbr}')
    return rp_fitness_tbr['log_fitness']

def get_replay_fitness_align(log, net, initial_marking, final_marking):
    rp_fitness_align = pm4py.fitness_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Alignments: {rp_fitness_align}')
    return rp_fitness_align['log_fitness']

# Precision zwischen Log und Modell berechnen (mit TBR und Alignments)
def get_precision_tbr(log, net, initial_marking, final_marking):
    precision_tbr = pm4py.precision_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Token-Based Replay beträgt: {precision_tbr}')
    return precision_tbr

def get_precision_align(log, net, initial_marking, final_marking):
    precision_align = pm4py.precision_alignments(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Alignments beträgt: {precision_align}')
    return precision_align

# F-Score zwischen Log und Modell aus ermittelten Fitness- und Precision-Werten berechnen
def get_f_score(fitness, precision):
    if fitness + precision == 0: # Division durch 0 vermeiden
        return 0.0
    f_score = 2 * (fitness * precision) / (fitness + precision) # F-Score anhand von Formel berechnen
    print(f'Der F-Score beträgt: {f_score}')
    return f_score

# Generalization und Simplicity berechnen
def get_generalization(log, net, initial_marking, final_marking):
    generalization = generalization_evaluator.apply(log, net, initial_marking, final_marking)
    print(f'Die Generalization beträgt: {generalization}')

def get_simplicity(net):
    simplicity = simplicity_evaluator.apply(net)
    print(f'Die Simplicity beträgt: {simplicity}')

# Übergabe von Arbeit ermitteln und anzeigen
def get_handover_of_work(log):
    handover_values = pm4py.discover_handover_of_work_network(log)
    pm4py.view_sna(handover_values)
    return handover_values # Ausgabe für die Cluster-Analyse nötig

# Ermitteln und Anzeigen, wie oft Subcontracting vorkommt
def get_subcontracting(log):
    subcont_values = pm4py.discover_subcontracting_network(log)
    pm4py.view_sna(subcont_values)

# Ermitteln und Anzeigen, wie oft zusammengearbeitet wird
def get_working_together(log):
    work_together_values = pm4py.discover_working_together_network(log)
    pm4py.view_sna(work_together_values)

# Ähnlichkeiten der Arbeitsmuster zwischen Individuen ermitteln und anzeigen
def get_similar_activities(log):
    similar_act = pm4py.discover_activity_based_resource_similarity(log)
    pm4py.view_sna(similar_act)

# Orginisationale Rollen entdecken und ausgeben
def get_orga_roles(log):
    roles = pm4py.discover_organizational_roles(log)

    print('Organisationale Rollen:')
    for role in roles:
        print(f'Aktivitäten: {role.activities}')
        print(f'Verantwortliche Rolle: {role.originator_importance}')
        important_role = max(role.originator_importance, key=role.originator_importance.get)
        value = role.originator_importance[important_role]
        print(f'Verantwortliche Rolle: {important_role} (Häufigkeit: {value})')
        print('-' * 40)


# Name der anonymen Ressourcen in Buchstaben umwandeln und Mapping ausgeben
def rename_org_resource(log):

    def is_anonymous(val):
        return bool(re.fullmatch(r'[a-z0-9]+', str(val)))
    
    anonymous = [val for val in log['org:resource'].unique() if is_anonymous(val)] # Liste der anonymen Ressourcen
    
    def label_generator():
        
        n = 1
        while True:
            for s in product(ascii_uppercase, repeat=n): # Buchstaben-Kombinationen generieren
                yield ''.join(s)
            n += 1

    labels = label_generator()
    mapping = {val: next(labels) for val in anonymous} # Mapping der anonymen Ressourcen zu Buchstaben
    
    log['org:resource'] = log['org:resource'].map(lambda x: mapping[x] if x in mapping else x) # Umbenennen der anonymen Ressourcen im Log
    
    print("Mapping anonymer Namen zu Buchstaben:")
    print(mapping)
    return mapping


# Cluster-Analyse nach Übergabe von Arbeit durchführen und anzeigen

def cluster_handover(handover): # Funktioniert nur für den Sepsis-Datensatz
    clustering = util.cluster_affinity_propagation(handover)
    print(f'Cluster nach Übergabe von Arbeit: {clustering}')
    
    colours = ['green', 'yellow', 'orange', 'red', 'blue', 'purple', 'brown'] # Farben für die Cluster
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for colour, (x, ys) in zip(colours, clustering.items()): # Erstellen eines Scatterplots der Cluster
        ax.scatter([x] * len(ys), ys, c=colour, linewidth=0, s=50)
    plt.show()

# Herausfiltern aller Cases des Sepsis-Datensatzes, die mit ER-Registration beginnen (zur Überprüfung der 1. Hypothese)
def filter_by_start_activities(log, start_activities, starts_to_remove):
    starts_to_keep = [act for act in start_activities.keys() if act not in starts_to_remove] # Alle Startaktivitäten behalten, außer den zu entfernenden
    filtered_log_start  = pm4py.filter_start_activities(log, starts_to_keep)
    return  filtered_log_start