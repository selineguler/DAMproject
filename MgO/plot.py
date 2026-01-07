from eos import volumes, vfine

outdir = Path(sys.argv[1])
energies=[]
for subdir in sorted(outdir.iterdir()):
  if subdir.is_dir():
    energy=[]
    for dat in sorted(subdir.glob("*.dat")):
      with open(dat, 'r') as f:
          lines = f.readlines()
      for line in lines:
          data = line.split()
          E = float(data[1])
          energy.append(E)
    energies.append(energy)


bulks = []
primed_bulks = []
v0_list = []
e0_list = []


for subdir in sorted(outdir.iterdir()):
    bulk_file = subdir / "bulk.out"
    if bulk_file.exists():
        with open(bulk_file, 'r') as f:
            for line in f:
                if line.startswith("Bulk"):
                    bulks.append(float(line.split()[1]))
                elif line.startswith("Primed bulk"):
                    primed_bulks.append(float(line.split()[2]))
                elif line.startswith("v0"):
                    v0_list.append(float(line.split()[1]))
                elif line.startswith("e0"):
                    e0_list.append(float(line.split()[1]))


def pressure(V, V0, K0, K0_prime):
    P = (3 * K0 / 2) * (np.power(V0 / V, 7/3.0) - np.power(V0 / V, 5/3.0)) * \
        (1 + (3/4) * (K0_prime - 4) * (np.power((V0 / V), 2/3) - 1))
    return P

def pressure5th(V,V0, K0, K0_prime,K0_prime2, K0_prime3):
  f = (1/2)*(np.power(V0 / V, 2/3) - 1)
  a1 = (3/2)*(K0_prime-4)
  a2 = (3/2)*(K0*K0_prime2 + K0_prime*(K0_prime - 7) + 143/7)
  a3 = (1/8)*(9*np.power(K0_prime,2)*K0_prime3 + 12*(3*K0_prime - 8)*K0*K0_prime2 + K0_prime*(118 + np.power(3*K0_prime-16,2)) - 1888/3)
  P=3*K0*f*np.power(1+2*f, 5/2)*(1+a1*f+a2*np.power(f,2)*a3*np.power(f,3))
  return P

def P_2003(K0, K0_prime, V, V0_oganov):
  Pres = - ((3/2)*K0*(xi-1)*np.power(x,-5/3)) - ((3/2)*K0*(1-2*xi)*np.power(x,-7/3)) - ((3/2)*K0*xi*np.power(x,-10/3))
  return Pres


def E(V, V0_oganov, E0, K0, K0_prime):

  #2003 oganov
  K0 *= GPa
  x = V0_oganov / V
  xi = (3/4) * (K0_prime - 4)

  term1 = (3/2) * (xi - 1) * np.power(x, (2/3))
  term2 = (3/4) * (1 - 2 * xi) * np.power(x,(4/3))
  term3 = (1/2) * xi * np.power(x,(6/3))
  constant_term = (2 * xi - 3) / 4
  E_x = E0 + (3/2) * K0 * V0_oganov * (term1 + term2 + term3 - constant_term)

  return E_x


pressures = []
pressuresd3=[]
pressuresd4=[]
pressuresPAW=[]
pressuresECP=[]
pressure5=[]
pressure0K=[]
pressurestat=[]
pressuresC=[]
pressuresCd3 = []



for vol in Vfine:
  PC = pressure(vol, v0_list[0], bulks[0], primed_bulks[0])
  PCd3 = pressure(vol, v0_list[1], bulks[1], primed_bulks[1])
  P = pressure(vol, v0_list[2], bulks[2], primed_bulks[2])
  Pd3 = pressure(vol, v0_list[3], bulks[3], primed_bulks[3])
  Pd4 = pressure(vol, v0_list[4], bulks[4], primed_bulks[4])
  P_Paw_Large = pressure(vol, 76.049/4,  154.183, 4.141)
  P_Ecp_Large = pressure(vol, 77.629/4, 151.707, 4.212)
  P_5th = pressure5th(vol, 18.76,  161.5, 4.0, -0.026, 0.0013)
  pressures.append(P)
  pressuresd3.append(Pd3)
  pressuresd4.append(Pd4)
  pressuresPAW.append(P_Paw_Large)
  pressuresECP.append(P_Ecp_Large)
  pressure5.append(P_5th)
  pressuresC.append(PC)
  pressuresCd3.append(PCd3)



E_static = E(Vfine, 73.425/4, e0_list[2], 181.240, 3.997)
E_ECP_large = E(Vfine, 77.629/4, e0_list[2], 151.707, 4.212)
E_PAW_large = E(Vfine, 76.049/4, e0_list[2], 154.183, 4.141)
E_0k = E(Vfine, 74.439/4, e0_list[2], 173.480, 4.014)
ECfit = E(Vfine, v0_list[0], e0_list[0], bulks[0], primed_bulks[0])
ECfitD3 = E(Vfine, v0_list[1], e0_list[1], bulks[1], primed_bulks[1])
Efit = E(Vfine, v0_list[2], e0_list[2], bulks[2], primed_bulks[2])
EfitD3 = E(Vfine, v0_list[3], e0_list[3], bulks[3], primed_bulks[3])
EfitD4 = E(Vfine, v0_list[4], e0_list[4], bulks[4], primed_bulks[4])





# EXPERIMENTAL VALUES

exp_pressures=[2.89, 6.73, 11.02, 15.82, 21.17, 27.15, 33.84, 41.32, 49.67, 59.02,
             69.49, 81.22, 87.61, 94.38, 101.55, 109.15, 117.22, 125.77, 134.85, 144.49,
             154.73, 165.60, 177.16, 189.46, 202.54, 216.47, 231.31, 247.12, 263.98, 281.97,
             301.18, 321.70]

exp_volumes=[0.98, 0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.84, 0.82, 0.80,
        0.78, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.70, 0.69, 0.68,
        0.67, 0.66, 0.65, 0.64, 0.63, 0.62, 0.61, 0.60, 0.59, 0.58,
        0.57, 0.56]



fig, (axi, axi_E) = plt.subplots(1, 2, figsize=(14,8))

axi.plot(pressuresC, Vfine, label='T=0K CHGNet', linestyle='--', linewidth=2)
axi.plot(pressuresCd3, Vfine, label='T=0K CHGNet+d3', linestyle='--', linewidth=2)
axi.plot(pressures, Vfine, label='T=0K M3GNet', color='green',linestyle='--', linewidth=2)
#axi.plot(np.array(exp_pressures), np.array(exp_volumes)*v0_list[2], label='T=300K', color='red')
axi.plot(pressuresd3, Vfine, label='T=0K M3GNet+D3', color='red',linestyle='--', linewidth=2)
axi.plot(pressuresd4, Vfine, label='T=0K M3GNet+D4', color='blue',linestyle='--', linewidth=2)
axi.plot(pressuresPAW, Vfine, label='T=0K PAW Large-core MgO', color='gray')
axi.plot(pressuresECP, Vfine, label='T=0K ECP Large-core MgO', color='purple')
axi.plot(pressure5, Vfine, label='5th order DFT-GGA', color='pink')


axi_E.plot(volumes, energies[0], 'o', label='CHGNet')
axi_E.plot(volumes, energies[1], 'o', label='CHGNet+D3')
axi_E.plot(volumes, energies[2], 'o', label='M3GNet', color='green')
axi_E.plot(volumes, energies[3], 'o', label='M3GNet+D3', color='brown')
axi_E.plot(volumes, energies[4], 'o', label='M3GNet+D4', color='gray')
axi_E.plot(Vfine, Efit)
axi_E.plot(Vfine, EfitD3)
axi_E.plot(Vfine, EfitD4)
axi_E.plot(Vfine, ECfit)
axi_E.plot(Vfine, ECfitD3)
axi_E.plot(Vfine, E_ECP_large, label='ECP Large Core MgO', color='pink',linestyle='--', linewidth=2)
axi_E.plot(Vfine, E_PAW_large, label='PAW Large Core MgO', color='purple',linestyle='--', linewidth=2)
axi_E.plot(Vfine, E_static, label='2003 Oganov Static', color='red',linestyle='--', linewidth=2)
axi_E.plot(Vfine, E_0k, label='2003 Oganov 0K', color='blue',linestyle='--', linewidth=2)

axi.set_xlabel(r'$P$ (GPa)', fontsize=16)
axi.set_ylabel(r'$V(A^3)$', fontsize=16)
axi.tick_params(axis='y', labelsize=12)
axi.tick_params(axis='x', labelsize=12)
axi_E.set_ylabel(r'$E$ (eV)', fontsize=16)
axi_E.set_xlabel(r'$V(A^3)$', fontsize=16)
axi_E.tick_params(axis='y', labelsize=12)
axi_E.tick_params(axis='x', labelsize=12)
axi.legend()
axi_E.legend()
axi.set_title('Pressure vs Volume', fontsize=14)
axi_E.set_title('Energy vs Volume', fontsize=14)

plt.subplots_adjust(hspace=0.5)
plt.tight_layout()
#plt.savefig(pdffil)
plt.show()
print("----------------------------")