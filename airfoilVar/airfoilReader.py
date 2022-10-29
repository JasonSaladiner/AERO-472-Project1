from csv import reader

def readAirfoil(fileName):
    with open(fileName, 'r') as f:
        ready = False
        foilReader = reader(f)
        Cl_curve = []
        Cd_curve = []
        for row in foilReader:
            if len(row) == 0:
                continue
            if row[0] == "Alpha":
                ready = True
                continue
            if not ready:
                continue
            Cl_curve.append(float(row[1]))
            Cd_curve.append(float(row[2]))
            #print(row)
    return Cl_curve, Cd_curve
#print(readAirfoil("n0009.csv"))