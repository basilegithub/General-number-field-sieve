# This is the file to parse the configuration file

import configparser

def parse_config(path):
    
    config = configparser.ConfigParser()
    
    config.read(path)
    
    return [config['DEFAULT'][key] for key in config['DEFAULT']]