import json

def load_protocol(file:str | None):
    default = 'ensemble_analyser/parameters_file/default_protocol.json'
    return json.load(open(default if not file else file ))

def load_threshold(file:str | None):
    default = 'ensemble_analyser/parameters_file/default_threshold.json'
    return json.load(open(default if not file else file ))
