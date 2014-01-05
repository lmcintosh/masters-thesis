#makeHostFile.py
import numpy as np


numCornServers = 30
cornServers    = range(1,numCornServers+1)
taus           = np.linspace(0.01,40,numCornServers)
versions       = [1,2,3]


hostFile = open('host_file.txt','a')
for f in xrange(numCornServers):
    if cornServers[f] < 10:
        hostFile.write('lanemc@corn0' + str(cornServers[f]) + '.stanford.edu python ~/masters/python/mastersRun_chunking10e4_shotnoise.py ' + str(versions[0]) + ' ' + str(taus[f]) + ' </dev/null >nohup.out 2>&1 &' + '\n')
    else:
        hostFile.write('lanemc@corn' + str(cornServers[f]) + '.stanford.edu python ~/masters/python/mastersRun_chunking10e4_shotnoise.py ' + str(versions[0]) + ' ' + str(taus[f]) + ' </dev/null >nohup.out 2>&1 &' + '\n')

hostFile.close()
