import pygmsh
import random
import time

geom = pygmsh.opencascade.Geometry(
  characteristic_length_min=0.1,
  characteristic_length_max=0.1,
  )

random.seed(int(time.time()))

outerDisk = geom.add_disk([ 0.0,  0.0, 0.0], 2.0)
for i in range(random.randint(10,20)):
	x = random.random()*2.
	y = random.random()*2.
	quadrant = random.randint(1,4)
	if quadrant == 1:
		quadrant = [1,1]
	elif quadrant == 2:
		quadrant = [1,-1]
	elif quadrant == 3:
		quadrant = [-1,-1]
	elif quadrant == 4:
		quadrant = [-1,1]	
	size = random.random()
	innerDisk = geom.add_disk([ x*quadrant[0],  y*quadrant[1], 0.0], size*.75)
	outerDisk = geom.boolean_difference([outerDisk], [innerDisk])
# geom.extrude(outerDisk, [0, 0, 0.3])

points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom, geo_filename = "swiss_cheese.geo")