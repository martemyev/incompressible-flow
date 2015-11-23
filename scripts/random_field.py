#!/usr/bin/python

import random, struct

nx = 100
ny = 100
nz = 1

values=[]
for iz in range(nz):
  for iy in range(ny):
    for ix in range(nx):
	values.append(random.uniform(0.1, 0.5))

# write the properties into a binary file with a single precision
out = open('properties.bin', 'wb')
s = struct.pack('f'*len(values), *values)
out.write(s)
out.close()

