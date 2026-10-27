"""pxil: one sequence as a picture, one pixel per residue.

A pxil file is a PNG image whose pixel values are the residues of a sequence,
with the alphabet and other details in its metadata (see SPEC.md). This
module reads and writes them:

    import pxil

    sequence, meta = pxil.read("TP53-201_ENST00000269305.transcript.pxil.png")
    pxil.write("p53.protein.pxil.png", "MEEPQSDPSV", kind="protein", id="ENSP00000269305")

and prints them as FASTA from the command line:

    python -m pxil TP53-201_ENST00000269305.transcript.pxil.png
    python -m pxil --info TP53-201_ENSP00000269305.protein.pxil.png

Copyright (c) 2026 Chronomic Medicine Group, Department of Genomic Medicine,
University of Cambridge. MIT licence.
"""

from __future__ import annotations

import argparse
import colorsys
import contextlib
import datetime as _dt
import hashlib
import json
import sys

from PIL import Image, PngImagePlugin

__version__ = "0.1.0"
SPEC_VERSION = "0.1"
KINDS = ("gene", "transcript", "protein", "sequence")

# the standard alphabets (SPEC.md section 6)
DNA = "-ACMGRSVTWYHKDBN-acmgrsvtwyhkdbn"  # gpxil reference codes; +16 = soft-masked
PROTEIN = "-ACDEFGHIKLMNPQRSTVWYUOBZJX*"

REQUIRED = ("pxil", "kind", "alphabet", "length", "width")

# Taylor (1997) amino acid colours, as in the Gpxil browser
TAYLOR = {
    "A": (204, 255, 0), "R": (0, 0, 255), "N": (204, 0, 255), "D": (255, 0, 0), "C": (255, 255, 0),
    "Q": (255, 0, 204), "E": (255, 0, 102), "G": (255, 153, 0), "H": (0, 102, 255), "I": (102, 255, 0),
    "L": (51, 255, 0), "K": (102, 0, 255), "M": (0, 255, 0), "F": (0, 255, 102), "P": (255, 204, 0),
    "S": (255, 51, 0), "T": (255, 102, 0), "W": (0, 204, 255), "Y": (0, 255, 204), "V": (153, 255, 0),
    "U": (255, 120, 200), "*": (40, 40, 48),
}


class PxilError(ValueError):
    """Not a valid pxil file, or a sequence that cannot be written as one."""


def _dna_colours() -> list[tuple[int, int, int]]:
    """A red, C blue, G green, T yellow, N grey; ambiguity codes mixed; repeats paler."""
    base = {1: (255, 0, 0), 2: (0, 0, 255), 4: (0, 255, 0), 8: (255, 255, 0)}
    out = []
    for code in range(32):
        low = code & 15
        if low == 0:
            c = (255, 255, 255)
        elif low == 15:
            c = (128, 128, 128)
        else:
            parts = [base[b] for b in (1, 2, 4, 8) if low & b]
            c = tuple(round(sum(p[i] for p in parts) / len(parts)) for i in range(3))
        if code & 16 and low:
            c = tuple(round(v + (255 - v) * 0.45) for v in c)
        out.append(c)
    return out


def default_colours(alphabet: str) -> list[tuple[int, int, int]]:
    """Display colours for an alphabet: the standard ones get their usual colours."""
    if alphabet == DNA:
        return _dna_colours()
    if alphabet == PROTEIN:
        return [(255, 255, 255)] + [TAYLOR.get(ch, (110, 110, 120)) for ch in alphabet[1:]]
    n = max(1, len(alphabet) - 1)  # other alphabets: hues round the colour wheel
    return [(255, 255, 255)] + [tuple(round(255 * v) for v in colorsys.hsv_to_rgb(i / n, 0.8, 0.95)) for i in range(n)]


def write(path, sequence: str, *, kind: str = "sequence", alphabet: str | None = None, width: int | None = None, colours=None, **meta) -> dict:
    """Write `sequence` as a pxil file; returns the metadata written.

    `alphabet` defaults to the standard protein alphabet for proteins and the
    DNA one otherwise; `width` (residues per row) to 100. Extra keyword
    arguments (id, name, source, ...) go into the metadata.
    """
    if kind not in KINDS:
        raise PxilError(f"kind must be one of {', '.join(KINDS)}, not {kind!r}")
    alphabet = alphabet or (PROTEIN if kind == "protein" else DNA)
    if not 2 <= len(alphabet) <= 256:
        raise PxilError("an alphabet has 2 to 256 symbols")
    index: dict[str, int] = {}
    for i, ch in enumerate(alphabet):
        index.setdefault(ch, i)
    try:
        codes = bytes(index[ch] for ch in sequence)
    except KeyError as e:
        raise PxilError(f"{e.args[0]!r} is not in the alphabet {alphabet!r}") from None
    n = len(codes)
    width = int(width or min(max(n, 1), 100))
    if width < 1:
        raise PxilError("width must be positive")
    height = max(1, -(-n // width))
    im = Image.frombytes("P", (width, height), codes + bytes(width * height - n))
    pal = list(colours or default_colours(alphabet))
    if len(pal) < len(alphabet):
        raise PxilError("give a colour for every symbol of the alphabet")
    pal = (pal + [(0, 0, 0)] * 256)[:256]  # 256 entries: the PNG is always 8-bit
    im.putpalette([v for c in pal for v in c])
    out = {"pxil": SPEC_VERSION, "kind": kind, "alphabet": alphabet, "length": n, "width": width}
    out.update({k: v for k, v in meta.items() if v is not None})
    out.setdefault("created", _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat())
    out.setdefault("generator", f"pxil.py {__version__}")
    out.setdefault("sha256", hashlib.sha256(sequence.encode()).hexdigest())
    info = PngImagePlugin.PngInfo()
    info.add_itxt("pxil", json.dumps(out, ensure_ascii=False), zip=False)
    title = " ".join(str(out[k]) for k in ("name", "id") if out.get(k))
    if title:
        info.add_text("Title", title)
    im.save(path, format="PNG", pnginfo=info)
    return out


def read(path) -> tuple[str, dict]:
    """The sequence in a pxil file, and its metadata."""
    with Image.open(path) as im:
        if im.format != "PNG":
            raise PxilError(f"{path}: not a PNG file")
        raw = im.info.get("pxil")
        if raw is None:
            raise PxilError(f"{path}: no pxil metadata (an iTXt chunk with the keyword 'pxil')")
        try:
            meta = json.loads(raw)
        except json.JSONDecodeError as e:
            raise PxilError(f"{path}: the pxil metadata is not JSON: {e}") from None
        missing = [k for k in REQUIRED if k not in meta]
        if missing:
            raise PxilError(f"{path}: metadata lacks {', '.join(missing)}")
        if str(meta["pxil"]).split(".")[0] != SPEC_VERSION.split(".")[0]:
            raise PxilError(f"{path}: pxil version {meta['pxil']} is not one this module reads")
        alphabet, n, w = meta["alphabet"], int(meta["length"]), int(meta["width"])
        if w < 1 or im.size != (w, max(1, -(-n // w))):
            raise PxilError(f"{path}: a {im.size[0]} x {im.size[1]} image cannot hold {n} residues {w} to a row")
        if im.mode not in ("P", "L"):
            raise PxilError(f"{path}: pixels must be palette indices or grey levels, not {im.mode}")
        data = im.tobytes()[:n]
    if data and max(data) >= len(alphabet):
        raise PxilError(f"{path}: index {max(data)} is outside the alphabet ({len(alphabet)} symbols)")
    table = {i: ch for i, ch in enumerate(alphabet)}
    return "".join(table[i] for i in data), meta


def main(argv=None) -> int:
    """The command line: ``pxil FILE…`` (or ``python -m pxil FILE…``)."""
    prog = "python -m pxil" if sys.argv and str(sys.argv[0]).endswith(".py") else "pxil"
    ap = argparse.ArgumentParser(prog=prog, description="Print pxil files as FASTA, or their metadata.")
    ap.add_argument("files", nargs="+", help="pxil files")
    ap.add_argument("--info", action="store_true", help="print the metadata instead")
    ap.add_argument("--line", type=int, default=60, help="FASTA line length (default 60)")
    ap.add_argument("-o", "--output", help="write to this file instead of printing")
    a = ap.parse_args(argv)
    with open(a.output, "w", encoding="utf-8") if a.output else contextlib.nullcontext(sys.stdout) as out:
        for f in a.files:
            seq, meta = read(f)
            if a.info:
                print(json.dumps(meta, indent=2, ensure_ascii=False), file=out)
                continue
            label = " ".join(str(meta[k]) for k in ("id", "name", "kind") if meta.get(k))
            print(f">{label or f}", file=out)
            for i in range(0, len(seq), a.line):
                print(seq[i : i + a.line], file=out)
    if a.output:
        print(f"wrote {a.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
