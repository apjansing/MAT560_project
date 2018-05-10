// This code was created by pygmsh v4.1.5.
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.1;
s0 = news;
Disk(s0) = {0.0, 0.0, 0.0, 0.5};