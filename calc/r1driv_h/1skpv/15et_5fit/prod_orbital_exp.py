z1 = 0.01
zN = 30.0
num = 15
r = (zN/z1)**(1.0/(num-1))
zs = [z1*r**n for n in range(num)]
for z in zs:
    print z.real, ",", z.imag
