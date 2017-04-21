to_au = 1.0/27.2114
e0 = -1.10263
for ene in [10.0, 20.0, 50.0, 250.0]:
    print ene*to_au - e0
