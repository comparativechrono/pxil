# pxil: sequences as pixels

pxil is a small, permissive file format for **one sequence as a picture**: a
gene, a transcript, a protein or any other string of symbols, one pixel per
residue, row after row like lines of text.

A pxil file is a PNG image. Its pixel values are the residues, its metadata
gives the alphabet, length and row width, and its palette colours them, so
the same file opens in any image software and reads back as the exact
sequence. That makes pxil handy for quick analysis, comparison and sharing
of single sequences. It is not meant for keeping genomes: that is what
[gpxil](https://github.com/comparativechrono/Gpxil), the genome pixel image
format, is for.

The [Gpxil browser](https://comparativechrono.github.io/Gpxil-browser/)
downloads genes, transcripts and proteins as pxil files:

| File | What it holds |
|---|---|
| `TP53_ENSG00000141510.gene.pxil.png` | the gene's genomic sequence, 5′ to 3′, introns and soft-masked repeats included |
| `TP53-201_ENST00000269305.transcript.pxil.png` | the transcript's spliced sequence, 5′ to 3′ |
| `TP53-201_ENSP00000269305.protein.pxil.png` | the protein, one pixel per amino acid |

These three files are in [`examples/`](examples/).

## The format in brief

- A PNG, 8-bit indexed colour. Each pixel is one residue; its palette index is
  the residue's position in the file's **alphabet**.
- Residues fill rows of `width` pixels; the rest of the last row is index 0.
- An `iTXt` chunk called `pxil` holds JSON metadata: the format version,
  `kind` (`gene`, `transcript`, `protein` or `sequence`), `alphabet`, `length`
  and `width`, and optionally identifiers, coordinates, a SHA-256 of the
  sequence and anything else useful.
- Two standard alphabets: DNA uses gpxil's reference codes
  (`-ACMGRSVTWYHKDBN-acmgrsvtwyhkdbn`, lower case for soft-masked repeats);
  protein uses `-ACDEFGHIKLMNPQRSTVWYUOBZJX*`. Any other alphabet of up to
  256 symbols is allowed.
- Files are named `<label>.<kind>.pxil.png`, with the database identifier in
  the label. The `.png` ending means they open as pictures anywhere.

The full description is [`SPEC.md`](SPEC.md), also as a PDF,
[`pxil-spec.pdf`](pxil-spec.pdf).

## Reading and writing pxil files in Python

[`pxil.py`](pxil.py) is a single module that needs only Pillow:

```sh
pip install git+https://github.com/comparativechrono/pxil
```

```python
import pxil

sequence, meta = pxil.read("examples/TP53-201_ENSP00000269305.protein.pxil.png")
print(meta["id"], len(sequence), sequence[:10])   # ENSP00000269305 394 MEEPQSDPSV

pxil.write("my.transcript.pxil.png", "ATGGAGGAGCCGCAGTCAGATCC", kind="transcript", id="my-1")
```

From the command line it prints pxil files as FASTA, or their metadata:

```sh
python -m pxil examples/*.pxil.png
python -m pxil --info examples/TP53-201_ENST00000269305.transcript.pxil.png
```

Without the module, any PNG reader will do: read the palette indices, keep
the first `length`, and look each up in `alphabet` (see
[`SPEC.md`](SPEC.md), section 8).

## Tests

```sh
pip install pillow pytest
python -m pytest tests
```

## Licence

Copyright © 2026 Chronomic Medicine Group, Department of Genomic Medicine,
University of Cambridge. Released under the [MIT licence](LICENSE).
