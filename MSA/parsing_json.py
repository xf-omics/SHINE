# parsing json file from ensembl genetree
# extract common_name and scientific_name for species, gene and protein accession IDs, and aligned protein sequences
# output 5 lines per gene, the order is a mess

import json
import argparse
import pathlib


def create_parser():
    parser = argparse.ArgumentParser(
        description="parse aligned protein sequences from ensembl genetree (homology) for each gene"  # noqa
    )

    parser.add_argument(
        "gene_list",
        type=str,
        help="a list of ensembl gene stable id, starting with ENSG for human",
    )

    parser.add_argument(
        "input_dir",
        type=str,
        help="directory for each gene tree, downloaded from ensembl",
    )

#    parser.add_argument(
#        "gene_stable_id",
#        type=str,
#        help="ensembl gene stable id, starting with ENSG for human",
#    )

    parser.add_argument(
        "output_dir",
        type=str,
        help="directory for each alignment, parsed from json",
    )

    return parser


def item_generator(json_input, key, lookup_key, barcode):
    if isinstance(json_input, dict):
        if key in json_input.keys():
            v = json_input[key]
            yield from item_generator(v, key, lookup_key, barcode)
        else:
            for k, v in json_input.items():
                if k in lookup_key:
                    yield barcode, k, v
                else:
                    yield from item_generator(v, key, lookup_key, barcode)
    elif isinstance(json_input, list):
        i = 33
        for item in json_input:
            c = barcode + chr(i)
            if i<126:
                i = i + 1
            yield from item_generator(item, key, lookup_key, c)


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    dict_lookup = {'scientific_name': 'Homo sapiens', 'common_name': 'Human', 'accession': 'ENSP0', 'seq': "amino acids"}
    key = "children"

    file_input = open(args.gene_list,'r')

    for gene_id in file_input:
        gene_id = gene_id.rstrip('\n')
        file_json = args.input_dir + gene_id + ".json"
        f = open(file_json)

        try:
            data = json.load(f)
        except ValueError as err:
            print(gene_id)
            continue

        data = data["tree"]
        f.close()

        file_a3m = args.output_dir + gene_id + ".txt"
        fw = open(file_a3m,'w')

        for b, k, v in item_generator(data,key,dict_lookup.keys(),'0'):
            fw.write(k)
            fw.write(' ')
            fw.write(v)
            fw.write(' ')
            fw.write(b)
            fw.write('\n')

        fw.close()

    file_input.close()

