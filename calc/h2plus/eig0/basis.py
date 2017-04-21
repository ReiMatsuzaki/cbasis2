zs_s_huzinaga = [0.0285649, 0.0812406, 0.190537, 0.463925, 1.202518,
               3.379649, 10.60720, 38.65163, 173.5822, 1170.498]
print [3.0 * z for z in zs_s_huzinaga]

zs_p_huzinaga = [0.015442, 0.035652, 0.085676, 0.227763, 0.710128, 3.009711]
print [3.0 * z for z in zs_p_huzinaga]

zs_d_h = [15.0, 5.0, 5.0/3.0, 5.0/9.0]
print zs_d_h
zs_s_cen = [14.4/(3.0**n) for n in range(12)]
print zs_s_cen
zeta_d_cen = [10.0/(3.0**n) for n in range(9)]
print zeta_d_cen
