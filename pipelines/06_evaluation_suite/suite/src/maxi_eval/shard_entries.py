"""Print comma-separated entry lists, one line per shard (round-robin, Boltz entries only unless --all).  Used to split CPU stages."""
import argparse, json
ap = argparse.ArgumentParser(); ap.add_argument("--manifest", default="work/inputs/manifest.json"); ap.add_argument("-k", type=int, default=2); ap.add_argument("--all", action="store_true")
a = ap.parse_args()
m = json.load(open(a.manifest))
ents = [e for e, v in m["entries"].items() if (a.all or e.endswith("_boltz")) and any(r["arm"] == "tsa" for r in v["replicates"])]
shards = [ents[i::a.k] for i in range(a.k)]
for s in shards:
    if s:
        print(",".join(s))
