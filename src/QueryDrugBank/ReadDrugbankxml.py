def get_string_ids(stringdb_name):
    translate_dict = dict()
    with open(stringdb_name, "r") as f:
        head = True
        for line in f:
            if head:
                head = False
            else:
                line = line.strip().split("\t")
                translate_dict.update({line[1]: line[0]})
    return translate_dict

def write_interactions(drug, output_location, translate_dict):
    drug_name = drug.find(f"{ns}name").text

    targets_hits = get_interacting_genes(
        drug.findall(f"{ns}targets/{ns}target")
    )
    enzyme_hits = get_interacting_genes(
        drug.findall(f"{ns}enzymes/{ns}enzyme")
    )
    carriers_hits = get_interacting_genes(
        drug.findall(f"{ns}carriers/{ns}carrier")
    )
    transporter_hits = get_interacting_genes(
        drug.findall(f"{ns}transporters/{ns}transporter")
    )

    _write_interaction("targets", drug_name, output_location, targets_hits, translate_dict)
    _write_interaction("enzymes", drug_name, output_location, enzyme_hits, translate_dict)
    _write_interaction("carriers", drug_name, output_location, carriers_hits, translate_dict)
    _write_interaction("transporter", drug_name, output_location, transporter_hits, translate_dict)


def _write_interaction(query_type, drug_name, output_location, hits, translate_dict):
    rename = {
        "CYP2B": "CYP2B6",

    }
    with open(f"{output_location}/{query_type}/{drug_name}.csv", "w") as w:
        for hit in hits:
            if hit in rename:
                hit = rename[hit]
            w.write(f"{translate_dict[hit]}\n")


def get_interacting_genes(queries):
    human_genes = []
    for match in queries:
        organism = match.find(f"{ns}organism")
        if organism.text == "Humans":
            polypeptide = match.find(f"{ns}polypeptide")
            if polypeptide:
                gene_name = _get_gene_name(polypeptide)
                human_genes.append(gene_name)
            else:
                print(match.text)

    return human_genes


def _get_gene_name(polypeptide):
    name = polypeptide.find(f"{ns}gene-name")
    return name.text


from lxml import etree as ET

translate_dict = get_string_ids("data/stringdb/9606.protein.info.v11.5.txt")

with open("data/selected_clusters/All_selected.txt", "r") as f:
    terms = [l.strip().capitalize() for l in f.readlines()]

ns = "{http://www.drugbank.ca}"


tree = ET.parse("data/drugbank/full_database.xml")
root = tree.getroot()
selected_drugs = []

for drug in root:
    name = drug.find(f"{ns}name").text
    if name in terms:
        selected_drugs.append(drug)


output = "input/drugbank"
for drug in selected_drugs:
    write_interactions(drug, output, translate_dict)
