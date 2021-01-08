import csv, random

f = open("gemAl.csv", 'w')
w = csv.writer(f)

for r in [-1,1]:
  for s in [1,2]:
    if s == 1: nsec = 36
    if s == 2: nsec = 18
    if s == 2: continue
    for c in range(nsec):
      detNum = r*(s*100+c+1)
      dx = float("%1.3f"%random.gauss(0,0.2))
      dy = float("%1.3f"%random.gauss(0,0.2))
      dz = 0.0
      dpx = 0.0
      dpy = 0.0
      dpz = float("%1.5f"%random.gauss(0,0.004)) #chambers ~ 100cm tall: dpz = dx/50cm
      w.writerow([detNum, dx, dy, dz, dpx, dpy, dpz])
      
