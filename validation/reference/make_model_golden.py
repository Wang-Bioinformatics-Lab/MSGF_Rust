#!/usr/bin/env python3
"""
make_model_golden.py — turn the AUTHORITATIVE MS-GF+ plain-text model dump (produced by
NewRankScorer.writeParametersPlainText via jshell) into a structured golden for the Rust
.param reader (crate msgf-scorer).

The dump is the ground truth: Java parsed the binary and printed every value. We record scalars,
the charge histogram, partitions, precursor/fragment offset frequencies, and rank/error
distributions — as full lists where small, and as grand sums + counts + exact first-partition
samples where large. The Rust reader must reproduce these; combined with the binary's
`0x7FFFFFFF` terminator (an alignment proof), matching these aggregates means the parse is right.

Usage: make_model_golden.py <model.model.txt> <model.param> -o <out.golden.json>
"""
import argparse, hashlib, json


def parse_dump(path):
    lines = [l.rstrip("\n") for l in open(path)]
    i = 0
    d = {}

    def val(line, pfx):
        return line.split(pfx, 1)[1].strip()

    # scalar header
    while i < len(lines):
        l = lines[i]
        if l.startswith("#Activation Method:"): d["activation"] = val(l, ":")
        elif l.startswith("#Instrument type:"): d["instrument"] = val(l, ":")
        elif l.startswith("#Enzyme:"): d["enzyme"] = val(l, ":")
        elif l.startswith("#Protocol:"): d["protocol"] = val(l, ":")
        elif l.startswith("#Maximum mass error:"):
            t = val(l, ":").split()
            d["mme"] = {"value": float(t[0]), "ppm": t[1].lower() == "ppm"}
        elif l.startswith("Apply deconvolution:"): d["apply_deconvolution"] = val(l, ":") == "true"
        elif l.startswith("Deconvolution error tolerance:"): d["deconvolution_error_tolerance"] = float(val(l, ":"))
        elif l.startswith("#ChargeHistogram"): break
        i += 1

    # charge histogram
    n = int(lines[i].split("\t")[1]); i += 1
    d["charge_histogram"] = []
    for _ in range(n):
        c, cnt = lines[i].split("\t"); d["charge_histogram"].append([int(c), int(cnt)]); i += 1

    # partitions
    assert lines[i].startswith("#Partitions"); n = int(lines[i].split("\t")[1]); i += 1
    parts = []
    for _ in range(n):
        c, seg, m = lines[i].split("\t"); parts.append([int(c), int(seg), float(m)]); i += 1
    d["partitions"] = parts
    d["num_segments"] = max(p[1] for p in parts) + 1

    # precursor offset
    assert lines[i].startswith("#PrecursorOffsetFrequencyFunction"); n = int(lines[i].split("\t")[1]); i += 1
    poff = []
    for _ in range(n):
        c, rc, off, _tol, freq = lines[i].split("\t"); poff.append([int(c), int(rc), float(off), float(freq)]); i += 1
    d["precursor_off"] = poff

    # fragment offset frequencies
    assert lines[i].startswith("#FragmentOffsetFrequencyFunction"); nparts = int(lines[i].split("\t")[1]); i += 1
    frag_freq_sum = frag_off_sum = 0.0; frag_entries = 0; frag_sample = None
    for pidx in range(nparts):
        hdr = lines[i].split("\t"); k = int(hdr[4]); i += 1
        block = {}
        for _ in range(k):
            name, freq, off = lines[i].split("\t")
            block[name] = [float(freq), float(off)]
            frag_freq_sum += float(freq); frag_off_sum += float(off); frag_entries += 1; i += 1
        if pidx == 0: frag_sample = block
    d.update(frag_off_total_entries=frag_entries, frag_off_freq_sum=frag_freq_sum,
             frag_off_offset_sum=frag_off_sum, frag_off_sample=frag_sample)

    # rank distributions
    assert lines[i].startswith("#RankDistributions"); nparts = int(lines[i].split("\t")[1]); i += 1
    rank_sum = 0.0; rank_floats = 0; rank_sample = None; max_rank = None
    for pidx in range(nparts):
        hdr = lines[i].split("\t"); n_ions = int(hdr[4]); max_rank = int(hdr[5]); i += 1
        block = {}
        for _ in range(n_ions):
            parts_ = lines[i].split("\t"); name = parts_[0]
            freqs = [float(x) for x in parts_[1:]]
            block[name] = freqs; rank_sum += sum(freqs); rank_floats += len(freqs); i += 1
        if pidx == 0: rank_sample = block
    d.update(max_rank=max_rank, rank_dist_total_floats=rank_floats,
             rank_dist_freq_sum=rank_sum, rank_dist_sample=rank_sample)

    # error distributions (optional)
    d["error_scaling_factor"] = 0
    if i < len(lines) and lines[i].startswith("#ErrorDistributions"):
        esf = int(lines[i].split("\t")[1]); i += 1
        sig_sum = noise_sum = ionex_sum = 0.0; err_sample = None; nparts = len(parts); pidx = 0
        while i < len(lines) and lines[i].startswith("Partition"):
            main_ion = lines[i].split("\t")[4]; i += 1
            sig = [float(x) for x in lines[i].split("\t")[1:]]; i += 1
            noi = [float(x) for x in lines[i].split("\t")[1:]]; i += 1
            iex = [float(x) for x in lines[i].split("\t")[1:]]; i += 1
            sig_sum += sum(sig); noise_sum += sum(noi); ionex_sum += sum(iex)
            if pidx == 0:
                err_sample = {"main_ion": main_ion, "signal": sig, "noise": noi, "ion_existence": iex}
            pidx += 1
        d.update(error_scaling_factor=esf, error_signal_sum=sig_sum, error_noise_sum=noise_sum,
                 ion_existence_sum=ionex_sum, error_sample=err_sample)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump"); ap.add_argument("param"); ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    d = parse_dump(a.dump)
    d["file"] = a.param.split("/")[-1]
    d["sha256"] = hashlib.sha256(open(a.param, "rb").read()).hexdigest()
    d["compare"] = {
        "scalars": "exact/str; floats abs 1e-6",
        "partitions.parentMass": "abs 0.05 (f32)",
        "sums": "abs tol 1e-2 (f32 summation order)",
        "samples": "abs 1e-6 (exact f32 values)",
        "terminator": "reader must consume to 0x7FFFFFFF",
    }
    json.dump(d, open(a.out, "w"), indent=1)
    print(f"wrote {a.out}: {len(d['partitions'])} partitions, max_rank={d['max_rank']}, "
          f"rank_floats={d['rank_dist_total_floats']}, esf={d['error_scaling_factor']}")


if __name__ == "__main__":
    main()
