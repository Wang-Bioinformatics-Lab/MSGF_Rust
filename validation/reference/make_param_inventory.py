#!/usr/bin/env python3
"""
make_param_inventory.py — ground-truth inventory of the MS-GF+ scoring model (.param) files
(Rust crate msgf-scorer, model loader). For each model we freeze:
  - size + sha256      : file integrity / correct download
  - header enums       : (activation, resolution/instrument, enzyme, protocol) — the model's identity
  - header_prefix_hex  : the leading bytes before the enum block, for format documentation

The .param header format (reverse-engineered): a 4-byte prefix, then a run of length-prefixed
UTF-16BE strings — 1 length byte, then N*2 bytes of UTF-16BE chars. The Rust loader must read
the same identity tuple; the sha256 guards that the model bytes themselves never silently change.
"""
import hashlib, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "..", "data", "models")
OUT = os.path.join(HERE, "..", "golden", "models")
PREFIX = 4       # bytes before the enum block
N_ENUMS = 4      # activation, resolution/instrument, enzyme, protocol


def parse_header(b):
    off = PREFIX
    enums = []
    for _ in range(N_ENUMS):
        n = b[off]; off += 1
        enums.append(b[off:off + 2 * n].decode("utf-16-be"))
        off += 2 * n
    return enums, off


def main():
    os.makedirs(OUT, exist_ok=True)
    entries = []
    for name in sorted(os.listdir(MODELS)):
        if not name.endswith(".param"):
            continue
        b = open(os.path.join(MODELS, name), "rb").read()
        enums, data_off = parse_header(b)
        entries.append({
            "file": name,
            "size": len(b),
            "sha256": hashlib.sha256(b).hexdigest(),
            "header_prefix_hex": b[:PREFIX].hex(),
            "identity": {
                "activation": enums[0],
                "resolution_or_instrument": enums[1],
                "enzyme": enums[2],
                "protocol": enums[3] or "(default)",
            },
            "binary_data_offset": data_off,
        })
        print(f"  {name:32s} {len(b):>8} B  {enums}")
    out = {
        "note": "MS-GF+ scoring-model inventory; Rust msgf-scorer must read the same identity tuple",
        "header_format": "4-byte prefix, then 4 length-prefixed UTF-16BE strings",
        "compare": {"size": "exact", "sha256": "exact", "identity": "exact"},
        "models": entries,
    }
    with open(os.path.join(OUT, "param_inventory.golden.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote param_inventory.golden.json ({len(entries)} models)")


if __name__ == "__main__":
    main()
