"""Round trips and checks for the pxil reference module (see SPEC.md)."""

from __future__ import annotations

import json
import os
import struct
import subprocess
import sys
import zlib

import pytest
from PIL import Image, PngImagePlugin

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pxil  # noqa: E402

EXAMPLES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def chunks(path):
    """(type, data) of each PNG chunk, checking the CRCs."""
    with open(path, "rb") as fh:
        b = fh.read()
    assert b[:8] == b"\x89PNG\r\n\x1a\n"
    i, out = 8, []
    while i < len(b):
        (n,) = struct.unpack(">I", b[i : i + 4])
        kind, data = b[i + 4 : i + 8], b[i + 8 : i + 8 + n]
        (crc,) = struct.unpack(">I", b[i + 8 + n : i + 12 + n])
        assert crc == zlib.crc32(kind + data)
        out.append((kind.decode(), data))
        i += 12 + n
    return out


def test_dna_round_trip(tmp_path):
    seq = "ACGTNacgtnRYKMSWBDHV" * 7 + "ACG"
    p = tmp_path / "x.gene.pxil.png"
    meta = pxil.write(p, seq, kind="gene", id="ENSG00000000001", name="PIXL1", width=30)
    got, m = pxil.read(p)
    assert got == seq
    assert m == meta and m["length"] == len(seq) and m["width"] == 30 and m["alphabet"] == pxil.DNA
    with Image.open(p) as im:
        assert im.size == (30, 5) and im.mode == "P"
        assert im.tobytes()[len(seq) :] == bytes(30 * 5 - len(seq))  # padding is index 0
        assert im.tobytes()[0] == 1 and im.tobytes()[1] == 2 and im.tobytes()[5] == 0x11  # A, C, then soft-masked a
    kinds = [k for k, _ in chunks(p)]
    assert kinds.index("iTXt") < kinds.index("IDAT")  # metadata before the pixels
    ihdr = dict(chunks(p))["IHDR"]
    assert ihdr[8] == 8 and ihdr[9] == 3  # 8-bit indexed colour


def test_protein_and_other_alphabets(tmp_path):
    prot = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD*"
    p = tmp_path / "TP53-201_ENSP00000269305.protein.pxil.png"
    pxil.write(p, prot, kind="protein", id="ENSP00000269305", genetic_code="standard")
    got, m = pxil.read(p)
    assert got == prot and m["alphabet"] == pxil.PROTEIN and m["kind"] == "protein"
    assert m["sha256"] == __import__("hashlib").sha256(prot.encode()).hexdigest()
    rna = tmp_path / "r.pxil.png"
    pxil.write(rna, "ACGUUGCA", alphabet="-ACGUN", width=3)
    assert pxil.read(rna)[0] == "ACGUUGCA"
    assert dict(chunks(rna))["IHDR"][8] == 8  # still 8-bit with a small alphabet


def test_n_looks_like_no_base_or_mixture():
    # the classic colours of gpxiltools: N near-black and never paler, since the mean of C (blue) and T (yellow) is grey
    cols = pxil.default_colours(pxil.DNA)
    n = pxil.DNA.index("N")
    assert cols[n] == cols[pxil.DNA.index("n")] == (32, 32, 32)
    assert cols[pxil.DNA.index("Y")] == (128, 128, 128)  # C or T
    assert all(cols[i] != cols[n] for i in range(len(pxil.DNA)) if pxil.DNA[i] not in "Nn")
    assert cols[pxil.DNA.index("a")] == (255, 115, 115)  # repeats are paler


def test_greyscale_files_are_read(tmp_path):
    seq = "ACGT" * 5
    codes = bytes(pxil.DNA.index(c) for c in seq)
    im = Image.frombytes("L", (8, 3), codes + bytes(4))
    info = PngImagePlugin.PngInfo()
    info.add_itxt("pxil", json.dumps({"pxil": "0.1", "kind": "sequence", "alphabet": pxil.DNA, "length": 20, "width": 8}), zip=True)
    p = tmp_path / "g.pxil.png"
    im.save(p, format="PNG", pnginfo=info)
    assert pxil.read(p)[0] == seq


def test_errors(tmp_path):
    with pytest.raises(pxil.PxilError):
        pxil.write(tmp_path / "a.pxil.png", "ACGTX", kind="gene")  # X is not DNA
    with pytest.raises(pxil.PxilError):
        pxil.write(tmp_path / "a.pxil.png", "ACGT", kind="genome")
    plain = tmp_path / "plain.png"
    Image.new("P", (4, 1)).save(plain)
    with pytest.raises(pxil.PxilError, match="no pxil metadata"):
        pxil.read(plain)

    def fake(meta, size=(4, 1), pixels=b"\x01\x02\x04\x08"):
        p = tmp_path / "f.pxil.png"
        im = Image.frombytes("P", size, pixels)
        im.putpalette([v for i in range(256) for v in (i, i, i)])  # 256 colours: an 8-bit PNG
        info = PngImagePlugin.PngInfo()
        info.add_itxt("pxil", json.dumps(meta))
        im.save(p, format="PNG", pnginfo=info)
        return p

    ok = {"pxil": "0.1", "kind": "gene", "alphabet": pxil.DNA, "length": 4, "width": 4}
    assert pxil.read(fake(ok))[0] == "ACGT"
    assert pxil.read(fake({**ok, "pxil": "0.7", "future": True}))[0] == "ACGT"  # same major version
    with pytest.raises(pxil.PxilError, match="version"):
        pxil.read(fake({**ok, "pxil": "1.0"}))
    with pytest.raises(pxil.PxilError, match="lacks"):
        pxil.read(fake({k: v for k, v in ok.items() if k != "width"}))
    with pytest.raises(pxil.PxilError, match="cannot hold"):
        pxil.read(fake({**ok, "length": 9}))
    with pytest.raises(pxil.PxilError, match="outside the alphabet"):
        pxil.read(fake({**ok, "alphabet": "-AC"}))


def test_command_line(tmp_path):
    p = tmp_path / "t.transcript.pxil.png"
    pxil.write(p, "ATG" * 30, kind="transcript", id="ENST00000000011", name="PIXL1-201")
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = subprocess.run([sys.executable, "-m", "pxil", "--line", "50", str(p)], capture_output=True, text=True, cwd=root, check=True).stdout
    assert out.splitlines() == [">ENST00000000011 PIXL1-201 transcript", "ATG" * 16 + "AT", "G" + "ATG" * 13]
    # to a file, with -o; and the metadata with --info
    fa = tmp_path / "t.fa"
    run = subprocess.run([sys.executable, "-m", "pxil", str(p), "-o", str(fa)], capture_output=True, text=True, cwd=root, check=True)
    assert run.stdout == "" and "wrote" in run.stderr
    assert fa.read_text().splitlines()[0] == ">ENST00000000011 PIXL1-201 transcript"
    assert pxil.main(["--info", str(p), "-o", str(tmp_path / "t.json")]) == 0
    assert '"kind": "transcript"' in (tmp_path / "t.json").read_text()


def test_fasta_names_of_a_persons_sequences(tmp_path):
    # the reference's protein and each haplotype's, as gpxil2pxil writes them: one name each
    meta = {"kind": "protein", "id": "ENSP00000269305", "name": "TP53-201", "sample": "HG002"}
    pxil.write(tmp_path / "ref.protein.pxil.png", "MEEPQSDPSVEPPLSQETFS", **{k: v for k, v in meta.items() if k != "sample"})
    pxil.write(tmp_path / "h1.protein.pxil.png", "MEEPQSDPSVEPPLSQETFS", haplotype=1, notes=["differs from the reference: P72R"], **meta)
    pxil.write(tmp_path / "h2.protein.pxil.png", "MEEPQSDPSVEPPLSQETFS", haplotype=2, **meta)
    fa = tmp_path / "all.fa"
    assert pxil.main([str(tmp_path / f) for f in ("ref.protein.pxil.png", "h1.protein.pxil.png", "h2.protein.pxil.png")] + ["-o", str(fa)]) == 0
    assert [x for x in fa.read_text().splitlines() if x.startswith(">")] == [
        ">ENSP00000269305 TP53-201 protein",
        ">ENSP00000269305_HG002-hap1 TP53-201 protein",
        ">ENSP00000269305_HG002-hap2 TP53-201 protein",
    ]
    assert pxil.fasta_name({"haplotype": 2}, "x.pxil.png") == "hap2"
    assert pxil.fasta_name({}, "x.pxil.png") == "x.pxil.png"


def test_a_reader_that_stops_early(tmp_path):
    # pxil … | head: no traceback when the pipe closes
    p = tmp_path / "long.sequence.pxil.png"
    pxil.write(p, "ACGT" * 50000, kind="sequence", id="long")
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    run = subprocess.run(f'"{sys.executable}" -m pxil --line 10 "{p}" | head -n 2', shell=True, capture_output=True, text=True, cwd=root)
    assert run.stdout.splitlines() == [">long sequence", "ACGTACGTAC"]
    assert "Traceback" not in run.stderr and "BrokenPipe" not in run.stderr and "Exception ignored" not in run.stderr, run.stderr
    # (true reads nothing and is gone before anything is written: the write fails, every time)
    run = subprocess.run(f'"{sys.executable}" -m pxil "{p}" | true', shell=True, capture_output=True, text=True, cwd=root)
    assert "Traceback" not in run.stderr and "BrokenPipe" not in run.stderr and "Exception ignored" not in run.stderr, run.stderr


@pytest.mark.skipif(not os.path.isdir(EXAMPLES), reason="no examples")
def test_examples_from_the_gpxil_browser():
    files = sorted(f for f in os.listdir(EXAMPLES) if f.endswith(".pxil.png"))
    assert files
    for f in files:
        seq, meta = pxil.read(os.path.join(EXAMPLES, f))
        assert len(seq) == meta["length"] and meta["kind"] in f
        if "sha256" in meta:
            assert meta["sha256"] == __import__("hashlib").sha256(seq.encode()).hexdigest()
        if meta["kind"] == "protein":
            assert seq.startswith("M")
