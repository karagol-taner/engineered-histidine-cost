#!/usr/bin/env python3
import argparse, csv, gzip, os, sys, time, urllib.error, urllib.request

URL = "https://files.rcsb.org/download/%s.pdb.gz"

def identifiers(path):
    with open(path, encoding="utf-8") as fh:
        head = fh.readline()
    if "pdb" in [c.strip().lower() for c in head.split(",")]:
        with open(path, newline="", encoding="utf-8") as fh:
            return [r["pdb"].strip().upper() for r in csv.DictReader(fh)
                    if r.get("pdb", "").strip()]
    return [l.strip().upper() for l in open(path, encoding="utf-8") if l.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", default="../../data/ladder_geometry_structures.csv")
    ap.add_argument("--out", default="../../data/pdb")
    a = ap.parse_args()
    ids = identifiers(a.ids)
    os.makedirs(a.out, exist_ok=True)
    have, got, failed = 0, 0, []
    for i in ids:
        p = os.path.join(a.out, "%s.pdb.gz" % i)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            have += 1
            continue
        try:
            with urllib.request.urlopen(URL % i, timeout=30) as fh:
                blob = fh.read()
            with open(p, "wb") as out:
                out.write(blob)
            with gzip.open(p, "rt", errors="replace") as fh:
                if not any(l.startswith("ATOM") for l in fh):
                    os.remove(p)
                    failed.append((i, "no ATOM records"))
                    continue
            got += 1
            time.sleep(0.1)
        except urllib.error.HTTPError as e:
            failed.append((i, "HTTP %s" % e.code))
        except Exception as e:
            failed.append((i, str(e)[:60]))
    print("identifiers %d, already present %d, downloaded %d, failed %d"
          % (len(ids), have, got, len(failed)))
    for i, why in failed:
        print("   FAILED %s: %s" % (i, why))
    if failed:
        sys.exit("structures missing, the analysis cannot be run on a partial set")


if __name__ == "__main__":
    main()
