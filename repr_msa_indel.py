#!/usr/bin/env python3
# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import argparse
import pathlib
from os.path import exists
import os
import torch
import esm
from esm import pretrained
from Bio import SeqIO
import itertools
from typing import List, Tuple
import string
import numpy as np

os.environ['CUDA_LAUNCH_BLOCKING'] = '-1'
# os.environ["CUDA_VISIBLE_DEVICES"]="0,1,2"

def create_parser():
    parser = argparse.ArgumentParser(
        description="Extract per-token representations and model outputs for sequences in a FASTA file"  # noqa
    )

    parser.add_argument(
        "model_location",
        type=str,
        help="PyTorch model file OR name of pretrained model to download (see README for models)",
    )
    parser.add_argument(
        "transcript_file",
        type=pathlib.Path,
        help="transcript file on which to extract representations",
    )
    parser.add_argument(
        "input_dir",
        type=pathlib.Path,
        help="input a3m directory to generate representations",
    )
    parser.add_argument(
        "output_dir",
        type=pathlib.Path,
        help="output directory for extracted representations",
    )

    parser.add_argument(
        "--species_per_msa", type=int, default=1024, help="maximum number of species"
    )
    parser.add_argument(
        "--repr_layers",
        type=int,
        default=[-1],
        nargs="+",
        help="layers indices from which to extract representations (0 to num_layers, inclusive)",
    )
    parser.add_argument(
        "--truncate", 
        action="store_true", 
        help="Truncate sequences longer than 1024 to match the training setup"
    )

    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
    return parser


# This is an efficient way to delete lowercase characters and insertion characters from a string
deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)

def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    return sequence.translate(translation)

def read_msa(filename: str, nseq: int) -> List[Tuple[str, str]]:
    """ Reads the first nseq sequences from an MSA file, automatically removes insertions."""
    return [(record.description, remove_insertions(str(record.seq)))
            for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)]


def main(args):
    model, alphabet = pretrained.load_model_and_alphabet(args.model_location)
    # model.eval()
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")

    batch_converter = alphabet.get_batch_converter()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    assert all(
        -(model.num_layers + 1) <= i <= model.num_layers for i in args.repr_layers
    )
    repr_layers = [
        (i + model.num_layers + 1) % (model.num_layers + 1) for i in args.repr_layers
    ]

    with open(args.transcript_file) as file:
        transcripts = file.read().split()
        print(
           f"Processing {len(transcripts)} transcripts"
        )
    file.close()
    
    with torch.no_grad():
        for transcript in transcripts:
            i = 0
            while i>=0:
                filename = (args.input_dir / f"{transcript}_{i}.a3m")
                if not exists(filename):
                    break

                msa_data = [
                    read_msa(filename, args.species_per_msa),
                ]
            
                try:
                    labels, strs, toks = batch_converter(msa_data)
                except:
                    print(filename)
                    break

                if torch.cuda.is_available() and not args.nogpu:
                    toks = toks.cuda()
                # print(toks.size(), toks.dtype)
                # [1, #of species, #of amino acids]

                aa = toks.shape[2]
                if aa>1024:
                    print(f"{transcript}_{i} has {toks.shape[2]} amino acids")
                    continue

                out = model(toks, repr_layers=repr_layers)

#           logits = out["logits"].to(device="cpu")
                representations = {
                    layer: t.to(device="cpu") for layer, t in out["representations"].items()
                }

                args.output_file = (
                    args.output_dir / f"{transcript}_{i}.pt"
                )
                args.output_file.parent.mkdir(parents=True, exist_ok=True)
                result = {"label": transcript}
            # Call clone on tensors to ensure tensors are not views into a larger representation
            # See https://github.com/pytorch/pytorch/issues/1995
                    
                result["representations"] = {
                    layer: t[0, 0, 0 : toks.size()[2]].clone()
                    for layer, t in representations.items()
                }

                torch.save(
                    result,
                    args.output_file,
                )
                i = i + 1



if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)

