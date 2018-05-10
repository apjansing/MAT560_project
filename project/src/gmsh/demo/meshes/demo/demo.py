import pygmsh
import numpy as np

geom = pygmsh.built_in.Geometry()

# Draw a cross.
poly = geom.add_polygon([
    [ 0.0,  0.0, 0.0],
    [ 1.0,  0.0, 0.0],
    [ 1.0,  1.0, 0.0],
    [ 0.0,  1.0, 0.0]
    ],
    lcar=1
    )


p2 = np.array([
    [ 0.0,  0.0,  0.0],
    [ 0.5,  0.5, -0.5],
    [ 1.0,  1.0,  0.0],
    [ 0.0,  1.0,  0.0],
    [ 0.5,  0.5, -0.5],
    [ 1.0,  0.0,  0.0],
    ])

# Draw a cross.
poly2 = geom.add_polygon(p2,
    lcar=1
    )

flip = lambda x : np.multiply(x, np.array([1,1,-1]))
p3 = [flip(p) for p in p2 ]
raize = lambda x : np.add(x, np.array([0,0,1]))
p3 = [raize(p) for p in p3 ]

# Draw a cross.
poly3 = geom.add_polygon(p3,
    lcar=1
    )


axis = [0, 0, 1]

geom.extrude(
    poly,
    translation_axis=axis,
    # rotation_axis=axis,
    point_on_axis=[0, 0, 0],
    angle=2.0 / 6.0 * np.pi
    )

points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom, geo_filename = "demo.geo")