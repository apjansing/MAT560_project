import matplotlib.pyplot as plt
import numpy as np

points = []
headerSkipped = False
with open('plate/verticesToEigenvalues.csv', 'r') as F:
	for line in F:
		if not headerSkipped:
			headerSkipped = True
		else:
			point = line.split(',')
			points += [[int(point[0]), int(point[1])]]
points = np.array(points).T

# fig, ax = plt.subplots(1, 1)

plt.plot(points[0], points[1], 'bo' )
plt.plot(points[0], points[1])

plt.xticks(np.arange(0, 2500.1, step=250))
plt.xlabel('Number of Veritces')

plt.yticks(np.arange(0, 5000.1, step=500))
plt.ylabel('Number of Eignenvalues')


plt.grid(True)
plt.show()