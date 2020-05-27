"""
Functions for importing and exporting the yeast model using COBRA from anywhere in the repo.
"""

from cobra.io import read_sbml_model, write_sbml_model
from dotenv import find_dotenv
from os.path import dirname

# find .env + define paths:
dotenv_path = find_dotenv()
repo_path = dirname(dotenv_path)
MODEL_PATH = f"{repo_path}/ModelFiles/xml/yeastGEM.xml"

def read_yeast_model():
    """Reads the SBML file of the yeast model using COBRA.

    Returns
    -------
    cobra.core.Model
    """
    model = read_sbml_model(MODEL_PATH)
    return model

def write_yeast_model(model):
    """Writes the SBML file of the yeast model using COBRA.

    Parameters
    ----------
    model : cobra.core.Model
        Yeast model to be written
    """
    write_sbml_model(model, MODEL_PATH)
