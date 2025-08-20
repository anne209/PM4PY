import pm4py
import pandas as pd
from pm4py.objects.conversion.log import converter as log_converter



def import_xes(path_to_xes: str):
    """Liest ein XES-Log ein"""
    log = pm4py.read_xes(path_to_xes)
    return log

log_sepsis = import_xes('sepsis_case.xes')
print(log_sepsis)


df_sepsis = log_converter.apply(log_sepsis, variant=log_converter.Variants.TO_DATA_FRAME)
print(df_sepsis.head()) 

print("Anzahl der Spalten:", len(df_sepsis.columns))


start_activities = pm4py.get_start_activities(log_sepsis)
print("Startaktivitäten mit Häufigkeiten: ")
print(start_activities)

end_activities = pm4py.get_end_activities(log_sepsis)
print("Endaktivitäten und Häufigkeiten: ")
print(end_activities)

cases_sepsis = len(log_sepsis['case:concept:name'].unique())
print(cases_sepsis)

def filter_log(start_acts, log, min_ratio=0.1):
    selected_activities = set()

    for activity, count in start_acts.items():
        if count / cases_sepsis >= min_ratio:
            selected_activities.add(activity)
    
    filtered_start = pm4py.filter_start_activities(log, selected_activities)
    return filtered_start

filter_start = filter_log(start_activities, log_sepsis)
print(filter_start)

dfg, start_activities, end_activities = pm4py.discover_dfg(log_sepsis)
dfg_output = pm4py.view_dfg(dfg, start_activities, end_activities)