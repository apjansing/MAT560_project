import pygmsh

geom = pygmsh.opencascade.Geometry(
  characteristic_length_min=0.1,
  characteristic_length_max=0.1,
  )

outerDisk = geom.add_disk([ 2.0,  0.0, 0.0], 1.0)
innerDisk = geom.add_disk([ 2.0,  0.0, 0.0], 0.5)
flat = geom.boolean_difference([outerDisk], [innerDisk])
# geom.extrude(flat, [0, 0, 0.3])

points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom, geo_filename = "donut.geo")

print(geom.get_code())