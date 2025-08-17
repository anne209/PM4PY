
# Libraries importieren
import pm4py
# import pandas as pd
# import os

# Funktion zum Aufrufen von .xes-Dateien
def import_xes(path_to_xes: str):
    """Liest ein XES-Log ein"""
    log = pm4py.read_xes(path_to_xes)
    return log # noch in df umwandeln?


# Aufruf des Sepsis-Datensatzes und Ausgabe einer Zusammenfassung
log_sepsis = import_xes('D:/Users/domin/Desktop/Projekt PM&Py/Datensätze/Sepsis/Sepsis Cases - Event Log.xes')
log_sepsis
log_sepsis.describe()
print('Head des Datensatzes: {}\nTail des Datensatzes: {}'.format(log_sepsis.head(10), log_sepsis.tail(10)))


# Anzahl der Fälle und Ereignisse anzeigen
cases_sepsis = len(log_sepsis['case:concept:name'].unique())
events_sepsis = len(log_sepsis)
print('Anzahl Fälle: {}\nAnzahl Ereignisse: {}'.format(cases_sepsis, events_sepsis))

# Start- und Endaktivitäten anzeigen
start_act_sepsis = pm4py.get_start_activities(log_sepsis)
end_act_sepsis = pm4py.get_end_activities(log_sepsis)
print('Startaktivitäten: {}\nEndaktivitäten: {}'.format(start_act_sepsis, end_act_sepsis))

# Filtern nach Startaktivitäten
def filter_by_start_activities(log, starts_to_keep):
    log_start = pm4py.get_start_activities(log)  # {activity: count}
    filtered_log_start  = pm4py.filter_start_activities(log, starts_to_keep)
    return log_start, filtered_log_start


log_start, filtered_log_start = filter_by_start_activities(log_sepsis, ['ER Registration'])

# Varianten auflisten und filtern
variants = pm4py.get_variants(log_sepsis)
print('Varianten: {}'.format(variants))

k = 5
filtered_log_var = pm4py.filter_variants_top_k(log_sepsis, k)
filtered_log_var