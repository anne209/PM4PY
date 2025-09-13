# Introduction

This repository is part of a university project that uses Process Mining for Python (PM4PY) to obtain information about real data event logs.

Two event logs were used for analysis. One event log concerns a sepsis process in a hospital, the other concerns the processing of applications from farmers in the European Union.

The event log of the sepsis process is uploaded to this repository under “sepsis_case.xes.” As long as it is in the same path as the code, no changes are necessary because the path is already used in the code.

The event log of the farmers' application process is much more complex and has a size of about 2 GB. If you want to use the entire event log, please visit the following website: https://data.4tu.nl/articles/dataset/BPI_Challenge_2018/12688355/1. You will find the file for download at the bottom of the page. Please note that not every device can handle such large amounts of data.
You must adjust the path in line 232 to use the original event log.
If the data is too large for your device, you can use the event log “iacs_case.xes”, which is only 5% of the actual event log but is much easier to process. The path to this event log is already preset. You do not need to make any adjustments as long as it is located in the same path.


Please note that some imports must be made before executing the code. These can all be found in the first lines of the import section. Please check for yourself which imports are required on your device, as some may already be installed on your device.

In addition, many images and histograms are displayed as they are part of the analysis. You will also receive some CSV files, as some of the analysis results are too large to be output in the program in which you are running the code. If you would like further analysis, please check the CSV files, as they offer opportunities for further analysis.
