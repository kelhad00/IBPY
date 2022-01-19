# IBPY
This repository provides tools for Interaction Behavior Processing.

This work represents the expressions occurences in a list of tuples each tuple representing an event's (such as an expression) occurence by (starting time, ending time, label of the expression). The starting and ending times are expressed in ms.

The following modules are provided in this repo:
* db.py : interfaces to access data of specific datasets.
* extract_data.py : get eaf files and extract data from them.
* interaction_model.py : encapsulte expressions, subjects and an entier interaction in OOC.
* interaction_analysis.py : functions for the analysis of interactions (mirroring coefficients, sequences occurences, etc.)
* processing.py : transforming expressions list.
* utils.py : other useful functions.
* visualization.py : functions used for the visualization of annotation of sequencial data.