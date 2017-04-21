import pandas as pd

df = pd.read_csv("search.csv")


df0 = df[df['re_0']<df['re_1']]
df1 = df[df['re_0']>df['re_1']].ix[:, ['re_1', 'im_1', 're_0', 'im_0']]
df1.columns = ['re_0', 'im_0', 're_1', 'im_1']
df =  pd.concat([df0, df1]).sort_values(by=["re_0"])

eps = pow(10.0, -4.0)
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
        print k, z0, z1, abs(z0_old-z0)+ abs(z1_old-z1)
        if(abs(z0_old-z0) + abs(z1_old-z1) >eps):
            if(z0.real < 1.0 and z0.imag > -1.0 and z1.real < 1.0 and z0.imag > -1.0):
                res.append([z0, z1])
        
    val_old = val

f = open('filtered.csv', 'w')
for [z0, z1] in res:
    print >>f, "{0},{1},{2},{3}".format(z0.real, z0.imag, z1.real, z1.imag)
    
"""
for i in range(len(df)-1):
    print i
    z0  = df.ix[i,:]['re_0'] + 1.0j * df.ix[i,:]['im_0']
    z1  = df.ix[i,:]['re_1'] + 1.0j * df.ix[i,:]['im_1']
    z0p = df.ix[i+1,:]['re_0'] + 1.0j * df.ix[i+1,:]['im_0']
    z1p = df.ix[i+1,:]['re_1'] + 1.0j * df.ix[i+1,:]['im_1']
    if(abs(z0-z0p)>eps and abs(z1-z1p)>eps):
        res.append([z0p, z1p])

print res
print len(df)

"""
