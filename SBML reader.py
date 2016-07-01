import cobra
import json
import libsbml
import os



file_path = os.path.join('C:\\','Users','user','PycharmProjects','SimParam','models', 'BIOMD0000000001.xml')


def libSBML_read():

    with open(file_path, 'r') as file:
        SBML = file.read()
        SBML = libsbml.readSBMLFromString(SBML)

    model = SBML.getModel()

    species = model.getListOfSpecies()

    print(species.get(1))

def cobra_read():

    model = cobra.io.sbml.create_cobra_model_from_sbml_file(file_path, old_sbml=False, legacy_metabolite=False, print_time=False, use_hyphens=False)

    species = model.metabolites
    species_id = model.metabolites.list_attr("id")

    reactions = model.reactions
    reactions_id = reactions.list_attr("id")

    query = species.get_by_id(species_id[1])

    JSON = json.loads(cobra.io.to_json(model))

    print(json.dumps(JSON, sort_keys=True, indent=4))

cobra_read()