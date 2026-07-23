#!/usr/bin/env python3
"""
make_spectra_golden.py — ground-truth fixtures for the spectrum-reading layer (Rust crate
msgf-io). For each input file we freeze per-spectrum parse facts plus a canonical hash of the
peak list, so the Rust MGF/mzML readers are checked to parse byte-for-byte the same values.

Canonical peak hash: sha1 over lines "%.5f %.5f" % (mz, intensity), one per peak, in file order,
joined by "\n". The Rust reader must format identically to reproduce the hash.

Big files are bounded: full per-spectrum records for the first N, plus a whole-file aggregate
(count, total peaks, and a rolling hash over every spectrum's peak hash) so nothing is unchecked.
"""
import base64, hashlib, json, os, struct, sys, xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "..", "golden", "spectra")
FULL_CAP = 150  # per-spectrum records to store in full for large files


def peak_hash(peaks):
    h = hashlib.sha1()
    h.update("\n".join(f"{mz:.5f} {it:.5f}" for mz, it in peaks).encode())
    return h.hexdigest()


def summarize(idx, title, scan, charge, prec_mz, peaks):
    mzs = [p[0] for p in peaks]
    its = [p[1] for p in peaks]
    base_i = max(range(len(peaks)), key=lambda i: its[i]) if peaks else -1
    return {
        "index": idx, "title": title, "scan": scan, "charge": charge,
        "precursor_mz": prec_mz, "n_peaks": len(peaks),
        "mz_min": min(mzs) if mzs else None, "mz_max": max(mzs) if mzs else None,
        "base_peak_mz": mzs[base_i] if peaks else None,
        "tic": round(sum(its), 4),
        "first3_peaks": [[round(m, 5), round(i, 5)] for m, i in peaks[:3]],
        "peaks_sha1": peak_hash(peaks),
    }


def parse_mgf(path):
    idx = 0
    title = scan = charge = prec = None
    peaks = []
    in_ions = False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                in_ions = True; title = scan = charge = prec = None; peaks = []
            elif line == "END IONS":
                yield summarize(idx, title, scan, charge, prec, peaks)
                idx += 1; in_ions = False
            elif not in_ions or not line:
                continue
            elif line.startswith("TITLE="):
                title = line[6:]
            elif line.startswith("SCANS="):
                scan = line[6:]
            elif line.startswith("PEPMASS="):
                prec = float(line[8:].split()[0])
            elif line.startswith("CHARGE="):
                charge = int(line[7:].rstrip("+"))
            elif line[0].isdigit():
                a = line.split()
                peaks.append((float(a[0]), float(a[1]) if len(a) > 1 else 0.0))


def _lt(tag):  # strip xml namespace
    return tag.rsplit("}", 1)[-1]


def _decode_binary(bda):
    """Decode one mzML <binaryDataArray> -> (kind, list[float])."""
    accs = {c.get("accession") for c in bda if _lt(c.tag) == "cvParam"}
    is64 = "MS:1000523" in accs
    zlib_c = "MS:1000574" in accs
    kind = "mz" if "MS:1000514" in accs else ("intensity" if "MS:1000515" in accs else "other")
    b64 = next((c.text for c in bda if _lt(c.tag) == "binary"), "") or ""
    raw = base64.b64decode(b64)
    if zlib_c:
        import zlib
        raw = zlib.decompress(raw)
    fmt = "<%dd" % (len(raw) // 8) if is64 else "<%df" % (len(raw) // 4)
    return kind, list(struct.unpack(fmt, raw)) if raw else []


def parse_mzml(path):
    for idx, sp in enumerate(e for _, e in ET.iterparse(path) if _lt(e.tag) == "spectrum"):
        cv = {c.get("accession"): c.get("value") for c in sp.iter() if _lt(c.tag) == "cvParam"}
        title = sp.get("id")
        arrays = {}
        for bda in (b for b in sp.iter() if _lt(b.tag) == "binaryDataArray"):
            kind, vals = _decode_binary(bda)
            arrays[kind] = vals
        peaks = list(zip(arrays.get("mz", []), arrays.get("intensity", [])))
        # precursor m/z (MS:1000744) and charge (MS:1000041) if present
        prec = float(cv["MS:1000744"]) if "MS:1000744" in cv else None
        charge = int(cv["MS:1000041"]) if "MS:1000041" in cv else None
        scan = title.split("scan=")[-1] if "scan=" in title else None
        yield summarize(idx, title, scan, charge, prec, peaks)
        sp.clear()


def build(path, parser, cap=None):
    records, count, total_peaks = [], 0, 0
    roll = hashlib.sha1()
    for rec in parser(path):
        count += 1
        total_peaks += rec["n_peaks"]
        roll.update(rec["peaks_sha1"].encode())
        if cap is None or count <= cap:
            records.append(rec)
    return {
        "file": os.path.basename(path),
        "format": "mgf" if path.endswith(".mgf") else "mzML",
        "n_spectra": count,
        "total_peaks": total_peaks,
        "rolling_peak_sha1": roll.hexdigest(),   # covers EVERY spectrum, not just the stored ones
        "stored_full": len(records),
        "compare": {"per_spectrum": "exact", "note": "n_peaks/charge/precursor_mz exact; peaks_sha1 must match canonical format"},
        "spectra": records,
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    jobs = [
        ("spectra/F13_subset25.mgf", parse_mgf, None,     "F13_subset25.golden.json"),
        ("spectra/F13.mgf",          parse_mgf, None,     "F13.golden.json"),
        ("spectra/test.mgf",         parse_mgf, FULL_CAP, "test_mgf.golden.json"),
        ("spectra/tiny.pwiz.mzML",   parse_mzml, None,    "tiny_mzML.golden.json"),
    ]
    for rel, parser, cap, out in jobs:
        p = os.path.join(DATA, rel)
        if not os.path.exists(p):
            print(f"  skip (missing): {rel}"); continue
        g = build(p, parser, cap)
        with open(os.path.join(OUT, out), "w") as fh:
            json.dump(g, fh, indent=2)
        print(f"  {rel:28s} -> {out}  ({g['n_spectra']} spectra, {g['total_peaks']} peaks, "
              f"{g['stored_full']} stored)")


if __name__ == "__main__":
    main()
