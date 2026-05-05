#! /bin/python

from os import listdir
from os.path import splitext
from cobra.io import read_sbml_model, save_matlab_model, save_json_model
import sys

model_path = sys.argv[1]
prefix = splitext(model_path)[0]

# Converting model
model = read_sbml_model(prefix + '.xml')
save_matlab_model(model, prefix + '.mat')
save_json_model(model, prefix + '.json')
