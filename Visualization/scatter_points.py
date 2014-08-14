import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def scatter_points(anglepath):
	#print anglepath
	xs = []
	ys = []
	zs = []
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1, projection='3d')
	for i, val in enumerate(anglepath):
		xs.append(val[0])
		ys.append(val[1])
		zs.append(val[2])
			
	ax.scatter(xs, ys, zs, c = 'r', marker = 'o')
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	plt.show()

