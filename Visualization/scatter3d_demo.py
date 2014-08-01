import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def plot_points(anglepath, target):
	#print anglepath
	xs = []
	ys = []
	zs = []
	x2 = []
	y2  = []
	z2 = []
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1, projection='3d')
	#ay = fig.add_subplot(2,1,2, projection = '3d')
	counter = 0
	for i, val in enumerate(anglepath):
		if counter == 0:
			xs.append(val[0])
			ys.append(val[1])
			zs.append(val[2])

	#		ax.plot(xs, ys, zs, c= 'r', marker='o')
			counter = counter + 1
		elif counter == 1:
			x2.append(val[0])
			y2.append(val[1])
			z2.append(val[2])
	#		ay.plot(xs, ys, zs, c='b', marker='o')
			counter = counter -1
	ax.scatter(xs, ys, zs, c = 'r', marker = 'o')
	ax.scatter(x2, y2, z2, c = 'r', marker = 'o')
#	ax.scatter([anglepath[0][0]], [anglepath[0][1]], [anglepath[0][2]], c = 'y', s = 100)
#	ax.scatter([anglepath[1][0]], [anglepath[1][1]], [anglepath[1][2]], c = 'y', s = 100)
#	ax.scatter([anglepath[-1][0]], [anglepath[-1][1]], [anglepath[-1][2]], c = 'g', s = 150)
#	ax.scatter([anglepath[-2][0]], [anglepath[-2][1]], [anglepath[-2][2]], c = 'g', s = 150)
#	for i, val in enumerate(target):
#		ax.scatter([val[0]], [val[1]], [val[2]], c = 'g', s = 60)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	plt.show()

