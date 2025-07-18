file = open("Clusters.txt",'r')
file.readline()
clust1 = []
for row in file:
    clusts = row.split(',')
    if int(clusts[14]) == 1:
        clust1.append(clusts[1])

print(clust1)


