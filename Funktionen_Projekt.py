# Libraries importieren
import pm4py
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn # Cluster-Analyse SNA funktioniert sonst nicht
import csv
from typing import List, Dict
import xml.etree.ElementTree as ET
from pm4py.algo.discovery.temporal_profile import algorithm as temporal_profile_discovery
from pm4py.visualization.petri_net import visualizer as pn_visualizer # für Zusatzinfos Performance
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

# DFG aus Log erstellen
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

# Umwandeln in Dataframe und Datetime, falls erforderlich
def ensure_df(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    return pm4py.convert_to_dataframe(obj)

def parse_ts(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors='coerce', utc=True).dt.tz_convert(None)

# Spaltenweises Zählen der Werte
def get_column_count(log, columns: list): # Übergabe muss eine Liste sein
    df = ensure_df(log) # Umwandeln in Data Frame, sofern erforderlich
    
    for col in columns: # Durchzählen und Ausgabe für alle übergebenen Spalten
        if col in df.columns:
            value_counts = df[col].value_counts()
            print(f'\nAnzahlen für Spalte {value_counts}')
        else:
            print(f'\nSpalte {col} existiert nicht')

# Alpha Miner anwenden
def run_alpha_miner(log):
    net, initial_marking, final_marking = pm4py.discover_petri_net_alpha(log)
    pm4py.view_petri_net(net, initial_marking, final_marking)
    return net, initial_marking, final_marking # Kann evtl. weg

# Heuristic Miner anwenden und in Petri-Netz umwandeln
def run_heuristic_miner(log, dependency_threshold=0.99):
    heuristic_net = pm4py.discover_heuristics_net(log, dependency_threshold)
    pm4py.view_heuristics_net(heuristic_net)
    heuristic_to_petri = pm4py.objects.conversion.heuristics_net.variants.to_petri_net.apply(heuristic_net, parameters=None)
    return heuristic_to_petri

# Inductive Miner anwenden
def run_inductive_miner(log, noise_threshold=0.0):
    inductive_net, initial_marking_inductive, final_marking_inductive = pm4py.discover_petri_net_inductive(log, noise_threshold)
    pm4py.view_petri_net(inductive_net, initial_marking_inductive, final_marking_inductive)
    return inductive_net, initial_marking_inductive, final_marking_inductive

# Varianten auflisten
def show_filter_variants(log):
    variants = pm4py.get_variants(log)
    with open('variants_output.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in variants.items():
            writer.writerow([key, value])
    print(f'Varianten: {variants}') # Ausdruck zu lang

# Filtern und Ausgeben eines Logs mit den 5 häufigsten Varianten sowie Visualisierung
def filter_by_variants(log, k):
    filtered_log_var = pm4py.filter_variants_top_k(log, k)
    return filtered_log_var

# Mittlere Durchlaufzeit pro Fall berechnen
def get_mean_throughput(log):
    durations = pm4py.get_all_case_durations(log)
    mean_throughput = np.mean(durations) / 3600 # Umrechnung in Stunden
    print(f'Die mittlere Durchlaufzeit pro Fall beträgt: {mean_throughput} Stunden')
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

# Performance-Informationen zu Netz aus Inductive Miner hinzufügen
def get_performance_net(log, net, initial_marking, final_marking):
    parameters = {pn_visualizer.Variants.PERFORMANCE.value.Parameters.FORMAT: "png"}
    gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters, variant=pn_visualizer.Variants.PERFORMANCE, log=log)
    pn_visualizer.view(gviz)

# Process Tree Inductive
def get_process_tree_ind(log):
    tree = pm4py.discover_process_tree_inductive(log)
    pm4py.view_process_tree(tree, format='png')

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
    plt.xlabel('Trace ist fit')
    plt.ylabel('Anzahl Traces')
    plt.title('Fit vs. Unfit (TBR)')
    plt.show()

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

# Throughput Analysis der falsch ausgeführten Transitionen und Ausgabe
def througput_trans(log, trans_fitness):
    trans_diagnostics = duration_diagnostics.diagnose_from_trans_fitness(log, trans_fitness)
    for trans in trans_diagnostics:
        print(f'Throughput Analysis der falsch ausgeführten Transitionen{trans, trans_diagnostics[trans]}')

# Throughput Analysis der nicht im Modell enthaltenen Aktivitäten und Ausgabe
def throughput_act(log, unwanted_activities):
    act_diagnostics = duration_diagnostics.diagnose_from_notexisting_activities(log, unwanted_activities)
    for act in act_diagnostics:
        print(act, act_diagnostics[act])

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

# Alignments bestimmen mit Modell aus Inductive Miner
def alignments_inductive(log, net, initial_marking, final_marking):
    aligned_traces_ind = pm4py.conformance_diagnostics_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Alignments mit Inductive Miner: {aligned_traces_ind}') # Ausdruck zu lang
    return aligned_traces_ind

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

# Fitness zwischen Log und Modell berechnen (TBR und Alignments)
def get_replay_fitness_tbr(log, net, initial_marking, final_marking):
    rp_fitness_tbr = pm4py.fitness_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Token-Based Replay: {rp_fitness_tbr}')
    return rp_fitness_tbr['log_fitness']

def get_replay_fitness_align(log, net, initial_marking, final_marking):
    rp_fitness_align = pm4py.fitness_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Alignments: {rp_fitness_align}')
    return rp_fitness_align['log_fitness']

# Precision zwischen Log und Modell berechnen (TBR und Alignments)
def get_precision_tbr(log, net, initial_marking, final_marking):
    precision_tbr = pm4py.precision_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Token-Based Replay beträgt: {precision_tbr}')
    return precision_tbr

def get_precision_align(log, net, initial_marking, final_marking):
    precision_align = pm4py.precision_alignments(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Alignments beträgt: {precision_align}')
    return precision_align

# F1-Score zwischen Log und Modell berechnen, Fitness und Precision gegeben
def get_f1_score(fitness, precision):
    if fitness + precision == 0: # Division durch 0 vermeiden
        return 0.0
    f1_score = 2 * (fitness * precision) / (fitness + precision)
    print(f'Der F1-Score beträgt: {f1_score}')
    return f1_score

# Generalization und Simplicity zwischen Log und Modell berechnen
def get_generalization(log, net, initial_marking, final_marking):
    generalization = generalization_evaluator.apply(log, net, initial_marking, final_marking)
    print(f'Die Generalization beträgt: {generalization}')

def get_simplicity(log, net, initial_marking, final_marking):
    simplicity = simplicity_evaluator.apply(net)
    print(f'Die Simplicity beträgt: {simplicity}')

# Erstellen eines Histogramms der Startaktivitäten
def plot_start_activities(log):
    df = ensure_df(log)
    first_acts = df.groupby('case:concept:name')['concept:name'].first() # Erste Aktivität pro Fall ermitteln
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
    last_acts = df.groupby('case:concept:name')['concept:name'].last() # Letzte Aktivität pro Fall ermitteln
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

# Erstellen eines Histogramms der Durchlaufzeiten der Cases
def plot_case_duration_hist(obj, bins: int = 30):
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

# Erstellen eines Diagramms der Events pro Tag
def plot_events_per_day(obj):
    df = ensure_df(obj)
    ts = 'time:timestamp'
    df[ts] = parse_ts(df[ts])
    daily = df.set_index(ts).resample('D').size()
    plt.figure()
    daily.plot()
    plt.xlabel('Datum')
    plt.ylabel('Events pro Tag')
    plt.title('Event-Zeitreihe (täglich)')
    plt.tight_layout()
    plt.show()

# Übergabe von Arbeit ermitteln und anzeigen
def get_handover_of_work(log):
    handover_values = pm4py.discover_handover_of_work_network(log)
    pm4py.view_sna(handover_values)
    return handover_values

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
    similar_act = pm4py.discover_activity_based_resource_similarity(log) # Benennung Var. wg Cluster-A.
    pm4py.view_sna(similar_act)

# Orginisationale Rollen entdecken und ausgeben
def get_orga_roles(log):
    roles = pm4py.discover_organizational_roles(log)
    print(f'Organisationale Rollen: {roles}')

# Cluster-Analyse nach ähnlichen Aktivitäten der Individuen

def cluster_similar_act(handover): # Evtl. noch nach anderen Sachen clustern
    clustering = util.cluster_affinity_propagation(handover)
    print(f'Cluster nach Übergabe von Arbeit: {clustering}')
    
    colours = ['green', 'yellow', 'orange', 'red', 'blue', 'purple', 'brown']
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for colour, (x, ys) in zip(colours, clustering.items()):
        ax.scatter([x] * len(ys), ys, c=colour, linewidth=0, s=50)
    plt.show() # aktuell nicht geordnet

# Herausfiltern aller Cases des Sepsis-Datensatzes, die mit ER-Registration beginnen (zur Überprüfung der 1. Hypothese)
def filter_by_start_activities(log, start_activities, starts_to_remove):
    starts_to_keep = [act for act in start_activities.keys() if act not in starts_to_remove]
    filtered_log_start  = pm4py.filter_start_activities(log, starts_to_keep)
    return  filtered_log_start