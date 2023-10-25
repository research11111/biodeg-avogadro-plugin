import argparse
import json
import sys
from rdkit import Chem
from biodeg import BioDegClassifier

def getOptions():

    return {
        'userOptions': {
            'Press OK to compute the biodegradability index': {
                'type': 'string',
                'default': '(unsued)'
            }
        },
        'inputMoleculeFormat': 'cjson'
    }

def jsonStrToMol(cjson_data):
    atom_numbers = cjson_data['atoms']['elements']['number']
    coordinates = cjson_data['atoms']['coords']['3d']
    bond_indices = cjson_data['bonds']['connections']['index']
    bond_orders = cjson_data['bonds']['order']

    mol = Chem.EditableMol(Chem.Mol())

    for atomic_num in atom_numbers:
        mol.AddAtom(Chem.Atom(atomic_num))

    for i in range(0, len(bond_indices), 2):
        mol.AddBond(bond_indices[i], bond_indices[i+1], Chem.BondType.values[bond_orders[i//2]])

    return mol.GetMol()

def run_command():
    stdinStr = sys.stdin.read()

    opts = json.loads(stdinStr)
    mol = jsonStrToMol(opts['cjson'])
    
    # This is a quick fix because the model has not been trained with enough data, any added information
    # make the compute false
    smiles = Chem.MolToSmiles(mol,allHsExplicit=False)
    mol = Chem.MolFromSmiles(smiles)

    c = BioDegClassifier.Prod()
    c.load()
    c.loadMols(mol) 
    result = c.guess()
    if 1 == len(result):
        key = next(iter(result.keys()))
        biodeg = c.biodeg_string_from_state(result[key])
        smiles = Chem.MolToSmiles(key,allHsExplicit=False)
        
        from PyQt5.QtWidgets import QDialog, QLabel, QVBoxLayout, QApplication

        app = QApplication([])

        dialog = QDialog()
        dialog.setWindowTitle("Biodegradability")

        layout = QVBoxLayout()
        label = QLabel("This molecule is %s\n" % (biodeg))
        layout.addWidget(label)

        dialog.setLayout(layout)
        dialog.show()

        app.exec_()

    else:
        # issue during processing
        pass
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser('BioDegradability')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    parser.add_argument('--run-command', action='store_true')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("BioDegradability")
    if args['menu_path']:
        print("Analytique")
    if args['print_options']:
        print(json.dumps(getOptions()))
    if args['run_command']:
        print(json.dumps(run_command()))

