
# Libraries importieren
import pm4py
import pandas as pd
import os
import matplotlib.pyplot as plt
from typing import List, Dict

# Funktion zum Import von .xes-Dateien
def import_xes(path_to_xes: str):
    log = pm4py.read_xes(path_to_xes)
    return log

# Import des Datensatzes
log_sepsis = import_xes('sepsis_case.xes')
log_sepsis.to_csv('log_sepsis.csv', index=False)
#log_iacs = import_xes('log_iacs.xes')
#log_iacs.to_csv('log_iacs.csv', index=False)


# Ausgabe statistischer Kennzahlen sowie von Head und Tail des Logs
def sum_up_log (log):
    print(f'Beschreibung des Datensatzes: \n{log.describe()}')
    print(f'Head des Datensatzes: \n{log.head(10)}\nTail des Datensatzes: \n{log.tail(10)}')
    
    
sum_up_log(log_sepsis)

# sum_up_log(log_iacs)

# Anzahl der Fälle und Ereignisse anzeigen
def get_cases_events(log):
    cases_no = len(log_sepsis['case:concept:name'].unique()) # case:... für iacs?
    events_no = len(log_sepsis)
    
    print(f'Anzahl Fälle: {cases_no}\nAnzahl Ereignisse: {events_no}')
    return cases_no

cases_no_sepsis = get_cases_events(log_sepsis)

# cases_no_iacs = get_cases_events(log_iacs)

# Start- und Endaktivitäten anzeigen
def get_start_end_act(log):
    start_act = pm4py.get_start_activities(log)
    end_act = pm4py.get_end_activities(log)
    
    print(f'Startaktivitäten: {start_act}\nEndaktivitäten: {end_act}')
    return start_act, end_act # Evtl. überflüssig

start_act_sepsis, end_act_sepsis = get_start_end_act(log_sepsis)

# start_act_iacs, end_act_iacs = get_start_end_act(log_iacs)

# DFG aus ungefiltertem Log erstellen
def create_dfg_from_log(log_input):
    dfg, start_activities, end_activities = pm4py.discover_dfg(log_input)
    dfg_output = pm4py.view_dfg(dfg, start_activities, end_activities)
    return dfg_output
    

dfg_sepsis_unfiltered = create_dfg_from_log(log_sepsis)

# dfg_iacs_unfiltered = create_dfg_from_log(log_iacs)


# Filteralgorithmus
def filter_log(start_acts, end_acts, log, no_of_cases, min_ratio=0.1, end_crit = None, delete_activities = None, check_value_activities=None):
    selected_activities = set()
    selected_end_acts = set()

    for activity, count in start_acts.items(): #Iteration im Dict über den Namen der Startaktivität und ihre Häufigkeit.
        if count / no_of_cases >= min_ratio: #Nur Startaktivitäten mit einem Mindestanteil von 10% an Gesamtfällen werden berücksichtigt.
            selected_activities.add(activity) #Sammeln der Aktivitäten in einer Liste
    
    if end_crit is not None:
        selected_end_acts = set(end_acts.keys()) - set(end_crit) #Löschen Fällen, die nicht ordungsgemäß mit Release oder Return ER geendet haben.
    
        drop_mask = pd.Series(False, index=log.index)
    
    for item in check_value_activities:

        is_activity = log["concept:name"] == item
        is_empty = log[item].isna()
        condition = is_activity & is_empty
        drop_mask = drop_mask | condition

    problematic_groups = log.loc[drop_mask, "org:group"].value_counts()
    print("Häufigkeit der Gruppen in fehlerhaften Zeilen:")
    print(problematic_groups)
    
    log.drop(index=log[drop_mask].index, inplace=True)
    

    filtered_log = pm4py.filter_event_attribute_values(log, 'concept:name', delete_activities, level='event', retain=False)
    filtered_log = pm4py.filter_start_activities(filtered_log, selected_activities)
    filtered_log = pm4py.filter_end_activities(filtered_log, selected_end_acts)

    return filtered_log

#Filterkriterien definieren
criteria_end_sepsis = ['IV Antibiotics', 'ER Sepsis Triage', 'Leucocytes', 'IV Liquid', 'CRP', 'LacticAcid', 'Admission NC', 'ER Triage'] #Endaktivitäten, die herausgefiltert werden sollen
deleted_activities_sepsis = ['Return ER'] #Aktivitäten die generell gelöscht werden sollen im Log
activity_check_sepsis = ["Leucocytes", "LacticAcid", "CRP"] #Aktivitäten die auf ihre Vollständigkeit im Log überprüft werden sollen

filtered_log_sepsis = filter_log(start_act_sepsis, end_act_sepsis, log_sepsis, cases_no_sepsis, end_crit = criteria_end_sepsis, 
                                 delete_activities = deleted_activities_sepsis, check_value_activities=activity_check_sepsis) #Übergabe der Parameter an die Filter-Methode
print(filtered_log_sepsis) # Vielleicht lieber nur Head/Tail etc.? (vor allem bei IACS-Datensatz)

filtered_log_sepsis.to_csv('filtered_log.csv', index=False)

# sum_up_log(filtered_log_sepsis)
# get_cases_events(filtered_log_sepsis)

'''criteria_end_iacs = ['PLATZHALTER']
deleted_activities_iacs = ["PLATZHALTER"]

filtered_log_iacs = filter_log(start_act_iacs, end_act_iacs, log_iacs, cases_no_iacs, end_crit = criteria_end_iacs, delete_activities = deleted_activities_iacs)
print(filtered_log_iacs)'''



# DFG aus gefiltertem Log erstellen
dfg_sepsis_filtered = create_dfg_from_log(filtered_log_sepsis)

# dfg_iacs_filtered = create_dfg_from_log(filtered_log_iacs)





# Filtern nach Startaktivitäten
def filter_by_start_activities(log, starts_to_keep):
    filtered_log_start  = pm4py.filter_start_activities(log, starts_to_keep)
    return filtered_log_start


filtered_log_start_sepsis = filter_by_start_activities(filtered_log_sepsis, ['ER Registration']) # aktuell nur für sepsis

# filtered_log_start_iacs = filter_by_start_activities(filtered_log_iacs, ['STARTAKTIVITÄT'])

# Alpha Miner anwenden
def run_alpha_miner(log):
    net, initial_marking, final_marking = pm4py.discover_petri_net_alpha(log)
    pm4py.view_petri_net(net, initial_marking, final_marking)
    return net, initial_marking, final_marking # Ausgabe für Conformance Checking notwendig

net_sepsis, initial_marking_sepsis, final_marking_sepsis = run_alpha_miner(filtered_log_sepsis) # Variablen evtl. umbenennen (auch bei TBR)

# net_iacs, initial_marking_iacs, final_marking_iacs = run_alpha_miner(filtered_log_iacs)

# Heuristic Miner anwenden
def run_heuristic_miner(log):
    heuristic_net = pm4py.discover_heuristics_net(log, dependency_threshold=0.99)
    pm4py.view_heuristics_net(heuristic_net)


run_heuristic_miner(filtered_log_sepsis)

# run_heuristic_miner(filtered_log_iacs)

# Inductive Miner anwenden
def run_inductive_miner(log):
    inductive_net, initial_marking_inductive, final_marking_inductive = pm4py.discover_petri_net_inductive(log)
    pm4py.view_petri_net(inductive_net, initial_marking_inductive, final_marking_inductive)
    return inductive_net, initial_marking_inductive, final_marking_inductive

ind_net_sepsis, initial_marking_ind_sepsis, final_marking_ind_sepsis = run_inductive_miner(filtered_log_sepsis)

# ind_net_iacs, initial_marking_ind_iacs, final_marking_ind_iacs = run_inductive_miner(filtered_log_iacs)

# Varianten auflisten
def show_filter_variants(log):
    variants = pm4py.get_variants(log)
    print(f'Varianten: {variants}')


show_filter_variants(filtered_log_sepsis)

# show_filter_vars(filtered_log_iacs)

# Filtern nach Varianten
def filter_by_variants(log, k):
    filtered_log_var = pm4py.filter_variants_top_k(log, k)
    return filtered_log_var


filtered_log_var_sepsis = filter_by_variants(filtered_log_sepsis, 5)

# filtered_log_var_iacs = filter_by_vars(filtered_log_iacs, 5)








# Helfer TBR
'''from typing import List, Dict
import pandas as pd'''

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


# TBR mit Modell aus Alpha Miner
def tbr_alpha(log, net, initial_marking, final_marking):
    replayed_traces = pm4py.conformance_diagnostics_token_based_replay(log, net, initial_marking, final_marking)
    df_tbr = tbr_list_to_dataframe(replayed_traces)
    print(f'Head des Token Based Replay mit Alpha Miner: \n{df_tbr.head(10)}')
    return df_tbr


df_tbr_sepsis = tbr_alpha(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)

# df_tbr_iacs = tbr_alpha(filtered_log_iacs, net_iacs, initial_marking_iacs, final_marking_iacs)

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

# Alignments bestimmen mit Modell aus Inductive Miner
def alignments_inductive(log, net, initial_marking, final_marking):
    aligned_traces_ind = pm4py.conformance_diagnostics_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Alignments mit Inductive Miner: {aligned_traces_ind}')
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


def get_replay_fitness_align(log, net, initial_marking, final_marking):
    rp_fitness_align = pm4py.fitness_alignments(log, net, initial_marking, final_marking)
    print(f'Ergebnisse für die Fitness mit Alignments: {rp_fitness_align}')


get_replay_fitness_tbr(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)
get_replay_fitness_align(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)

# get_replay_fitness_tbr(filtered_log_iacs, net_iacs, initial_marking_iacs, final_marking_iacs)
# get_replay_fitness_align(filtered_log_iacs, net_iacs, initial_marking_iacs, final_marking_iacs)

# Precision zwischen Log und Modell berechnen (ETConformance (TBR) und Align-ETConformance(Alignments))
def get_precision_tbr(log, net, initial_marking, final_marking):
    precision_tbr = pm4py.precision_token_based_replay(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Token-Based Replay beträgt: {precision_tbr}')


def get_precision_align(log, net, initial_marking, final_marking):
    precision_align = pm4py.precision_alignments(log, net, initial_marking, final_marking)
    print(f'Die Precision mit Alignments beträgt: {precision_align}')


get_precision_tbr(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)
get_precision_align(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)

# get_precision_tbr(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)
# get_precision_align(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)




'''# Generalization und Simplicity zwischen Log und Modell berechnen
from pm4py.algo.evaluation.generalization import algorithm as generalization_evaluator

def get_generalization(log, net, initial_marking, final_marking):
    generalization = generalization_evaluator.apply(log, net, initial_marking, final_marking)
    print(f'Die Generalization beträgt: {generalization}')


from pm4py.algo.evaluation.simplicity import algorithm as simplicity_evaluator

def get_simplicity(log, net, initial_marking, final_marking):
    simplicity = simplicity_evaluator.apply(net)
    print(f'Die Simplicity beträgt: {simplicity}')


get_generalization(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)
get_simplicity(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)

# get_generalization(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)
# get_simplicity(filtered_log_sepsis, net_sepsis, initial_marking_sepsis, final_marking_sepsis)'''












def ensure_df(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    return pm4py.convert_to_dataframe(obj)


def parse_ts(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors='coerce', utc=True).dt.tz_convert(None)







'''def plot_start_activities(obj):
    df = ensure_df(obj)'''
    










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





# Spalte umbenennen für Social Network Analysis
'''log_sepsis_sna = filtered_log_sepsis.rename(columns = {'org:group' : 'org:resource'})

# log_iacs_sna = filtered_log_iacs.rename(columns = {'ALTER SPALTENNAME' : 'org:resource'})'''

# Übergabe von Arbeit ermitteln und anzeigen
'''def get_handover_of_work(log):
    handover_values = pm4py.discover_handover_of_work_network(log)
    pm4py.view_sna(handover_values)


get_handover_of_work(log_sepsis_sna)

# get_handover_of_work(log_iacs_sna)''' # Fehler: Column not found: org:resource

# Ermitteln und Anzeigen, wie oft Subcontracting vorkommt
'''def get_subcontracting(log):
    subcont_values = pm4py.discover_subcontracting_network(log)
    pm4py.view_sna(subcont_values)


get_subcontracting(log_sepsis_sna)

# get_subcontracting(log_iacs_sna)''' # Fehler: Column not found: org:resource

# Ermitteln und Anzeigen, wie oft zusammengearbeitet wird
'''def get_working_together(log):
    work_together_values = pm4py.discover_working_together_network(log)
    pm4py.view_sna(work_together_values)


get_working_together(log_sepsis_sna)

# get_working_together(log_iacs_sna)''' # Fehler: Column not found: org:resource

# Ähnlichkeiten der Arbeitsmuster zwischen Individuen ermitteln und anzeigen
'''def get_similar_activities(log):
    similar_act = pm4py.discover_activity_based_resource_similarity(log) # Benennung Var. wg Cluster-A.
    pm4py.view_sna(similar_act)
    return similar_act


similar_activities_sepsis = get_similar_activities(log_sepsis_sna)

# similar_activities_iacs = get_similar_activities(log_iacs_sna)''' # Fehler: Column not found: org:resource

# Orginisationale Rollen entdecken und ausgeben
'''def get_orga_roles(log):
    roles = pm4py.discover_organizational_roles(log)
    print([x[0] for x in roles])


get_orga_roles(log_sepsis_sna)

# get_orga_roles(log_iacs_sna)''' # KeyError: org:resource

# Cluster-Analyse nach ähnlichen Aktivitäten der Individuen
'''from pm4py.algo.organizational_mining.sna import util # wahrscheinlich nicht gut

def cluster_similar_act(log, sim_act): # Evtl. noch nach anderen Sachen clustern
    clustering = util.cluster_affinity_propagation(sim_act)


cluster_similar_act(log_sepsis_sna, similar_activities_sepsis_sna)

# cluster_similar_act(log_iacs_sna, similar_activities_iacs_sna)''' # Var. von vorher nicht definiert


