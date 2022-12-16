import sys
import pandas
from chembl_webresource_client.new_client import new_client


class ExtraInfo:
    def __init__(self, targets):
        self.ids = targets
        self.targets = new_client.target
        self.target_gene_pairs = []

    def filter(self, target_id):
        filter_kwargs = {
            "target_chembl_id__exact": target_id,
            "organism__exact": "Homo sapiens",
        }
        print(f"Querying {target_id}")
        hits = self.targets.filter(**filter_kwargs)

        if len(hits) == 0:
            return None
        else:
            return hits[0]

    def unpack(self, protein):
        genes = []
        components = protein["target_components"]
        for component in components:
            synonyms = component["target_component_synonyms"]
            for synonym in synonyms:
                if synonym["syn_type"] == "GENE_SYMBOL":
                    genes.append(synonym["component_synonym"])

        if len(genes) == 0:
            raise NameError(f"No genes found.")
        return genes

    def append(self):
        for target_id in self.ids:
            hit = self.filter(target_id)
            if hit is None:
                self.target_gene_pairs.append((target_id, ["Non-human",]))
            else:
                genes = self.unpack(hit)
                self.target_gene_pairs.append((target_id, genes))

    def write_csv(self, output):
        with open(output, "w") as w:
            w.write(f"target\tgene\n")
            for target_id, genes in self.target_gene_pairs:
                for gene in genes:
                    w.write(f"{target_id}\t{gene}\n")


if __name__ == '__main__':
    args = sys.argv[1:]
    df = pandas.read_csv(args[0], sep="\t")
    ids = list(df["target"])
    ei = ExtraInfo(targets=ids)
    ei.append()
    ei.write_csv(args[1])
