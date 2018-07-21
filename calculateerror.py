from lammps import lammps
datafilename="data.rp3.2250new"
bondfilename="bonds.reaxc"
ffieldfilename="ffield.reax.cho"
lmp=lammps()
lmp.command("units real")
lmp.command("atom_style    charge")
lmp.command("read_data "+datafilename)
lmp.command("pair_style reax/c lmp_control")
lmp.command("pair_coeff * * "+ffieldfilename+" C H O")
lmp.command("compute reax all pair reax/c")
lmp.command("fix 3 all qeq/reax 1 0.0 10 1.0e-6 reax/c")
lmp.command("fix 5 all reax/c/bonds 100 "+bondfilename)
lmp.command("timestep 0.1")
lmp.command("run 0")
N=lmp.get_natoms()
number=[0 for x in range(N+1)]
bond=[[] for x in range(N+1)]
done=[False for x in range(N+1)]
qa=[0.0 for x in range(N+1)]
q=lmp.gather_atoms("q",1,1)
RMSD=[]
molecular=[]
m=0
lmp.close()
with open(bondfilename) as f:
    for line in f:
        if line[0]!="#":
            s=line.split()
            bid=int(s[0])
            number[bid]=int(s[2])
            for i in range(0,number[bid]):
                bond[bid].append(int(s[i+3]))
def mo(i,n):
    molecular[n].append(i)
    done[i]=True
    for b in bond[i]:
        if not done[b]:
            mo(b,n)
for i in range(1,N+1):
    if not done[i]:
        molecular.append([])
        mo(i,m)
        m+=1
for i in range(5,11):
    print(str(i))
    for m in molecular:
        lmp=lammps()
        lmp.command("units real")
        lmp.command("atom_style    charge")
        lmp.command("atom_modify map array")
        lmp.command("read_data "+datafilename)
        lmp.command("pair_style reax/c lmp_control")
        lmp.command("pair_coeff * * "+ffieldfilename+" C H O")
        lmp.command("compute reax all pair reax/c")
        lmp.command("variable i equal "+str(i))
        lmp.command("variable j equal 1")
        lmp.command("variable a1 equal x[v_j]")
        lmp.command("variable a2 equal y[v_j]")
        lmp.command("variable a3 equal z[v_j]")
        for atom in m:
            lmp.command("variable j equal "+str(atom))
            lmp.command("region nearby sphere ${a1} ${a2} ${a3} ${i}")
            lmp.command("group nearby region nearby")
            lmp.command("region nearby delete")
        lmp.command("fix 3 nearby qeq/reax 1 0.0 10 1.0e-6 reax/c")
        lmp.command("timestep 0.1")
        lmp.command("run 0")
        qb = lmp.gather_atoms("q",1,1)
        lmp.close()
        for atom in m:
            qa[atom-1]=qb[atom-1]
    qtotal=0
    for j in range(0,N):
        qtotal+=(q[j]-qa[j])**2
    RMSDnow=str(i)+"A, RMSD="+str((qtotal/N)**0.5)
    RMSD.append(RMSDnow)
    print(qa)
    print(RMSDnow)
with open("error.txt",'w') as f:
    for RMSDeach in RMSD:
        print(RMSDeach,file=f)
