import pandas as pd
au2ev = 27.2114
res = pd.read_csv("ref_cs.csv")
for w_ev in res["w"]:
    print w_ev/au2ev
