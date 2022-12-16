import sys

import onnxruntime
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from chembl_webresource_client.new_client import new_client

FP_SIZE = 1024
RADIUS = 2


class PredictInteraction:
    def __init__(self, name):
        self.smiles = self.get_smile(name)
        self.ort_session = onnxruntime.InferenceSession("data/chembl_multitask_model/trained_models/chembl_31_model/chembl_31_multitask.onnx")
        self.predictions = []

    def get_smile(self, name):
        molecule = new_client.molecule
        mols = molecule.filter(pref_name__iexact=name.upper())
        if len(mols) > 1:
            print(f"More than one compound found for {name}")
        if len(mols) == 0:
            raise IndexError(f"No compounds for {name}")
        else:
            mol = mols[0]
            try:
                if mol["structure_type"] == "MOL":
                    smile = mol["molecule_structures"]["canonical_smiles"]
                    return smile
                else:
                    return False
            except KeyError as e:
                raise KeyError(f"No smile for {name}")

    ## Code from https://github.com/chembl/chembl_multitask_model
    def calc_morgan_fp(self):
        mol = Chem.MolFromSmiles(self.smiles)
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=FP_SIZE)
        a = np.zeros((0,), dtype=np.float32)
        Chem.DataStructs.ConvertToNumpyArray(fp, a)
        return a

    ## Code from https://github.com/chembl/chembl_multitask_model
    def format_preds(self, preds, targets):
        preds = np.concatenate(preds).ravel()
        np_preds = [(tar, pre) for tar, pre in zip(targets, preds)]
        dt = [('chembl_id', '|U20'), ('pred', '<f4')]
        np_preds = np.array(np_preds, dtype=dt)
        np_preds[::-1].sort(order='pred')
        return np_preds

    def run_prediction(self):
        if self.smiles:
            descs = self.calc_morgan_fp()
            ort_inputs  = {self.ort_session.get_inputs()[0].name: descs}
            preds = self.ort_session.run(None, ort_inputs)
            self.predictions = self.format_preds(preds, [o.name for o in self.ort_session.get_outputs()])

    def output_csv(self, output_file):
        with open(output_file, "w") as w:
            w.write("target\tprediction\n")
            for target, prediction in self.predictions:
                w.write(f"{target}\t{prediction}\n")


if __name__ == '__main__':
    args = sys.argv[1:]
    drug = args[0]
    output = args[1]
    pi = PredictInteraction(drug)
    pi.run_prediction()
    pi.output_csv(output)