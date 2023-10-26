import argparse
import json
import sys, os
from rdkit import Chem
from biodeg import BioDegClassifier
from PyQt5.QtWidgets import QDialog, QLabel, QVBoxLayout, QApplication, QPushButton, QFileDialog
from PyQt5.QtCore import QStandardPaths

data = {
    'classifier': None,
    'result': None
}

def getOptions():

    return {
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


class RunResultGui():
    def __init__(self,biodeg):
        app = QApplication([])

        dialog = QDialog()
        dialog.setWindowTitle("Biodegradability")

        layout = QVBoxLayout()
        label = QLabel("This molecule is %s\n" % (biodeg))
        layout.addWidget(label)

        error_output = QLabel("")
        layout.addWidget(error_output)

        export_button = QPushButton("Export")
        export_button.clicked.connect(self.export_molecules_from_current_data)
        layout.addWidget(export_button)

        self.app = app
        self.export_button = export_button
        self.error_output = error_output

        dialog.setLayout(layout)
        dialog.show()
        self.app.exec_()

    def export_molecules_from_current_data(self):
        fileLocation = os.path.join(QStandardPaths.writableLocation(QStandardPaths.StandardLocation.DownloadLocation),"export.csv")
        file_name, _ = QFileDialog.getSaveFileName(None, "Save File", fileLocation, "File (*.csv);;All Files (*)", options=QFileDialog.Option.DontUseNativeDialog)
        if file_name:
            try:
                data['classifier'].guess_result_to_csv(file_name,data['result'])
                self.error_output.setText("Result writted to " + file_name)
            except OSError as e:
                self.error_output.setText("Error %s" % (str(e)))
        else:
            self.error_output.setText("No filename provided")

def run_command():
    stdinStr = sys.stdin.read()

    opts = json.loads(stdinStr)
    mol = jsonStrToMol(opts['cjson'])
    
    # This is a quick fix because the model has not been trained with enough data, any added information
    # make the compute false
    smiles = Chem.MolToSmiles(mol,allHsExplicit=False)
    mol = Chem.MolFromSmiles(smiles)

    c = BioDegClassifier.Prod()
    data['classifier'] = c
    c.load()
    c.loadMols(mol) 
    result = c.guess()
    result_dict = {}
    if 1 == len(result):
        data['result'] = result
        key = next(iter(result.keys()))
        biodeg = c.biodeg_string_from_state(result[key])
        smiles = Chem.MolToSmiles(key,allHsExplicit=False)
        result_dict['cjson'] = {
            "properties": { "biodegradability": "{biodeg}" }
        }
        #RunResultGui(biodeg)

    else:
        # issue during processing
        result_dict['cjson'] = {
            "properties": { "biodegradability": "Error during computation" }
        }
        pass
    return result_dict
    

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

