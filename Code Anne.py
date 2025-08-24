import pm4py
import pandas as pd
from pm4py.objects.conversion.log import converter as log_converter



def import_xes(path_to_xes: str):
    """Liest ein XES-Log ein"""
    log = pm4py.read_xes(path_to_xes)
    return log

log_sepsis = import_xes('sepsis_case.xes')
print(log_sepsis)

log_sepsis.to_csv('log_sepsis.csv', index=False)

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

def filter_log(start_acts, end_acts, log, min_ratio=0.1, end_crit=None, delete_activities=None, check_value_activities=None):
    selected_activities = set()
    selected_end_acts = set()



    for activity, count in start_acts.items():
        if count / cases_sepsis >= min_ratio:
            selected_activities.add(activity)
    
    if end_crit is not None:
        selected_end_acts = set(end_acts.keys()) - set(end_crit)

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
    log.to_csv('log_test.csv', index=False)



    filtered_log = pm4py.filter_event_attribute_values(log, "concept:name", delete_activities, level='event', retain=False)
    filtered_log = pm4py.filter_start_activities(filtered_log, selected_activities)
    filtered_log = pm4py.filter_end_activities(filtered_log, selected_end_acts)

    return filtered_log

criteria_end=["IV Antibiotics", "ER Sepsis Triage", "Leucocytes", "IV Liquid", "CRP", "LacticAcid", "Admission NC", "ER Triage"]
deleted_activities = ["Return ER"]
activity_check = ["Leucocytes", "LacticAcid", "CRP"]

filtered_log= filter_log(start_activities, end_activities, log_sepsis, end_crit=criteria_end, 
                         delete_activities=deleted_activities, check_value_activities=activity_check)
print(filtered_log)


start_activities = pm4py.get_start_activities(filtered_log)
print(start_activities)

end_activities = pm4py.get_end_activities(filtered_log)
print(end_activities)

dfg, start_activities, end_activities = pm4py.discover_dfg(log_sepsis)
dfg_output = pm4py.view_dfg(dfg, start_activities, end_activities)

filtered_log.to_csv('filtered_log.csv', index=False)
