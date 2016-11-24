

datafile = open('SingleElectron.dat')
content = datafile.read()

for line in content.split('\n'):
    if 'name' not in line:
        print line

datafile.close()





