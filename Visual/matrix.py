#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   matrix.py
#   library of operations on vectors and matrices
#   written by Georgy Samsonidze (October 2008)
#
#-------------------------------------------------------------------------------

import math
from common import eps9

def main(argv = None):
   if argv is None:
      argv = sys.argv
   sys.exit("\n   Module matrix.py contains the following names:\n" +
            "   align_vector       - align vector along z axis\n" +
            "   vector_dot_vector  - compute vector vector dot product\n" +
            "   vector_dot_matrix  - compute vector matrix dot product\n" +
            "   matrix_dot_vector  - compute matrix vector dot product\n" +
            "   matrix_dot_matrix  - compute matrix matrix dot product\n" +
            "   rotation_matrix    - generate rotation matrix\n" +
            "   transpose_matrix   - transpose matrix\n" +
            "   determinant_matrix - compute determinant of 3x3 matrix\n" +
            "   invert_matrix      - invert 3x3 matrix\n")

def align_vector(vector):
  number = 0
  axes = []
  angles = []
  axis = 2
  angle = math.atan2(vector[1], vector[0])
  if abs(angle) > eps9 and angle < math.pi - eps9 and angle > eps9 - math.pi:
    number += 1
    axes.append(axis)
    angles.append(angle)
    vector = matrix_dot_vector(rotation_matrix(axis, angle), vector)
  axis = 1
  angle = math.atan2(vector[0], vector[2])
  if abs(angle) > eps9 and angle < math.pi - eps9 and angle > eps9 - math.pi:
    number += 1
    axes.append(axis)
    angles.append(angle)
  return number, axes, angles

def vector_dot_vector(vector1, vector2):
   scalar3 = 0.0
   column1 = len(vector1)
   row2 = len(vector2)
   ierr = 0
   if row2 != column1:
      ierr = 1
   if ierr == 0:
      xx = 0.0
      for ii in range(column1):
         xx += vector1[ii] * vector2[ii]
      scalar3 = xx
   return scalar3

def vector_dot_matrix(vector1, matrix2):
   vector3 = []
   column1 = len(vector1)
   row2 = len(matrix2)
   column2 = len(matrix2[0])
   ierr = 0
   for ii in range(1, row2):
      if len(matrix2[ii]) != column2:
         ierr = 1
   if row2 != column1:
      ierr = 1
   if ierr == 0:
      for ii in range(column2):
         xx = 0.0
         for jj in range(column1):
            xx += vector1[jj] * matrix2[jj][ii]
         vector3.append(xx)
   return vector3

def matrix_dot_vector(matrix1, vector2):
   vector3 = []
   row1 = len(matrix1)
   column1 = len(matrix1[0])
   row2 = len(vector2)
   ierr = 0
   for ii in range(1, row1):
      if len(matrix1[ii]) != column1:
         ierr = 1
   if row2 != column1:
      ierr = 1
   if ierr == 0:
      for ii in range(row1):
         xx = 0.0
         for jj in range(column1):
            xx += matrix1[ii][jj] * vector2[jj]
         vector3.append(xx)
   return vector3

def matrix_dot_matrix(matrix1, matrix2):
   matrix3 = []
   row1 = len(matrix1)
   column1 = len(matrix1[0])
   row2 = len(matrix2)
   column2 = len(matrix2[0])
   ierr = 0
   for ii in range(1, row1):
      if len(matrix1[ii]) != column1:
         ierr = 1
   for ii in range(1, row2):
      if len(matrix2[ii]) != column2:
         ierr = 1
   if row2 != column1:
      ierr = 1
   if ierr == 0:
      for ii in range(row1):
         matrix3.append([])
         for jj in range(column2):
            xx = 0.0
            for kk in range(column1):
               xx += matrix1[ii][kk] * matrix2[kk][jj]
            matrix3[ii].append(xx)
   return matrix3

def rotation_matrix(axis, angle):
   if axis == 0:
      matrix = [[1.0, 0.0, 0.0], [0.0, math.cos(angle), math.sin(angle)], [0.0, -math.sin(angle), math.cos(angle)]]
   elif axis == 1:
      matrix = [[math.cos(angle), 0.0, -math.sin(angle)], [0.0, 1.0, 0.0], [math.sin(angle), 0.0, math.cos(angle)]]
   elif axis == 2:
      matrix = [[math.cos(angle), math.sin(angle), 0.0], [-math.sin(angle), math.cos(angle), 0.0], [0.0, 0.0, 1.0]]
   else:
      matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
   return matrix

def transpose_matrix(matrix1):
   matrix2 = []
   row1 = len(matrix1)
   column1 = len(matrix1[0])
   ierr = 0
   for ii in range(1, row1):
      if len(matrix1[ii]) != column1:
         ierr = 1
   if ierr == 0:
      for ii in range(column1):
         matrix2.append([])
         for jj in range(row1):
            matrix2[ii].append(matrix1[jj][ii])
   return matrix2

def determinant_matrix(matrix1):
   ierr = 0
   det1 = 0.0
   row1 = len(matrix1)
   column1 = len(matrix1[0])
   for ii in range(1, row1):
      if len(matrix1[ii]) != column1:
         ierr = 1
   if ierr == 0:
      if row1 == 3 and column1 == 3:
         det1 = matrix1[0][0] * (matrix1[1][1] * matrix1[2][2] - matrix1[1][2] * matrix1[2][1]) + matrix1[0][1] * (matrix1[1][2] * matrix1[2][0] - matrix1[1][0] * matrix1[2][2]) + matrix1[0][2] * (matrix1[1][0] * matrix1[2][1] - matrix1[1][1] * matrix1[2][0])
      else:
         ierr = 1
   return ierr, det1

def invert_matrix(matrix1):
   ierr = 0
   matrix2 = []
   row1 = len(matrix1)
   column1 = len(matrix1[0])
   for ii in range(1, row1):
      if len(matrix1[ii]) != column1:
         ierr = 1
   if ierr == 0:
      if row1 == 3 and column1 == 3:
         matrix2.append([])
         matrix2[0].append(matrix1[1][1] * matrix1[2][2] - matrix1[1][2] * matrix1[2][1])
         matrix2[0].append(matrix1[0][2] * matrix1[2][1] - matrix1[0][1] * matrix1[2][2])
         matrix2[0].append(matrix1[0][1] * matrix1[1][2] - matrix1[0][2] * matrix1[1][1])
         matrix2.append([])
         matrix2[1].append(matrix1[1][2] * matrix1[2][0] - matrix1[1][0] * matrix1[2][2])
         matrix2[1].append(matrix1[0][0] * matrix1[2][2] - matrix1[0][2] * matrix1[2][0])
         matrix2[1].append(matrix1[0][2] * matrix1[1][0] - matrix1[0][0] * matrix1[1][2])
         matrix2.append([])
         matrix2[2].append(matrix1[1][0] * matrix1[2][1] - matrix1[1][1] * matrix1[2][0])
         matrix2[2].append(matrix1[0][1] * matrix1[2][0] - matrix1[0][0] * matrix1[2][1])
         matrix2[2].append(matrix1[0][0] * matrix1[1][1] - matrix1[0][1] * matrix1[1][0])
         ierr, det1 = determinant_matrix(matrix1)
         if abs(det1) > eps9:
            for ii in range(row1):
               for jj in range(column1):
                  matrix2[ii][jj] = matrix2[ii][jj] / det1
         else:
            ierr = 1
      else:
         ierr = 1
   return ierr, matrix2

if __name__ == "__main__":
   import sys
   sys.exit(main())

