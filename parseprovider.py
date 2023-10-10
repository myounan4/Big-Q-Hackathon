Mbase = 100.0
def parseprovider(input_file):
    with open(input_file) as file:
        lines = [line.rstrip() for line in file]
    noofproviders = int(lines[0])
    demand = float(lines[1])/Mbase
    providerinfo = lines[2:]
    cost = []
    xmin = []
    xmax = []
    for ps in providerinfo:
        providerdetails = ps.split(',')
        cost.append(float(providerdetails[0])*Mbase)
        xmin.append(float(providerdetails[1])/Mbase)
        xmax.append(float(providerdetails[2])/Mbase)

    return noofproviders,demand,cost,xmin,xmax
