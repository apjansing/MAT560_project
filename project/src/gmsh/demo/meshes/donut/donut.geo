// This code was created by pygmsh v4.1.5.
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.1;
s0 = news;
Disk(s0) = {2.0, 0.0, 0.0, 1.0};
s1 = news;
Disk(s1) = {2.0, 0.0, 0.0, 0.5};
bo1[] = BooleanDifference{ Surface{s0}; Delete; } { Surface{s1}; Delete;};