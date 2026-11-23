"""Write Boltz-2 inputs (no MSA, protein A + ligand B = CCD 6NT) for every sequence in sequences.csv and split into chunks of CHUNK predictions.
usage: python3 prep_boltz.py [CHUNK=15]; skips sequences that already have a prediction in boltz_out/.   -> boltz_in/chunk_NN/<seq_id>.yaml, chunks.txt"""
import glob, os, sys, yaml, pandas as pd
H = os.path.expanduser("~/HACKATHON_ETH"); D = H + "/runs/20261121_r5_downstream"; CH = int(sys.argv[1]) if len(sys.argv) > 1 else 15
df = pd.read_csv(D + "/sequences.csv"); done = {os.path.basename(p) for p in glob.glob(D + "/boltz_out/chunk_*/boltz_results_*/predictions/*")}
df = df[~df.seq_id.isin(done) & df.sequence.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+")]
# existing chunk dirs keep their content; new chunks continue the numbering
start = len(glob.glob(D + "/boltz_in/chunk_*")); n = 0
for k in range(0, len(df), CH):
    cd = f"{D}/boltz_in/chunk_{start + k // CH:02d}"; os.makedirs(cd, exist_ok=True)
    for _, r in df.iloc[k:k + CH].iterrows():
        open(f"{cd}/{r.seq_id}.yaml", "w").write(yaml.safe_dump({"version": 1, "sequences": [{"protein": {"id": "A", "sequence": r.sequence, "msa": "empty"}}, {"ligand": {"id": "B", "ccd": "6NT"}}]}, sort_keys=False)); n += 1
print(n, "inputs in", -(-len(df) // CH), "new chunks (first new chunk index", start, ")")
