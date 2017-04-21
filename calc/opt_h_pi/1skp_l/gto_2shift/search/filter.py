import pandas as pd

df = pd.read_csv("search.csv").sort_values(by=["re_0"])

eps = pow(10.0, -3.0)
first = True
res = []
val_old = None
for (k, val) in df.iterrows():
    z0 = val['re_0']+1.0j*val['im_0']
    z1 = val['re_1']+1.0j*val['im_1']

    if (first):
        res.append([z0, z1])
        first = False
    else:
        z0_old = val_old['re_0']+1.0j*val_old['im_0']
        z1_old = val_old['re_1']+1.0j*val_old['im_1']
        print k, z0, z1, abs((z0_old-z0)/z0)+ abs((z1_old-z1)/z1)
        if(abs(z0_old-z0) + abs(z1_old-z1) >eps):
            res.append([z0, z1])
        
    val_old = val

f = open('filtered.csv', 'w')
for [z0, z1] in res:
    print >>f, "{0},{1},{2},{3}".format(z0.real, z0.imag, z1.real, z1.imag)
