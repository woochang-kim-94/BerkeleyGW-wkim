/*------------------------------------------------------------------------------

   surface.cpp
   generates an isosurface of the scalar field
   written by Georgy Samsonidze (October 2008)
   Objects used by ICM refactored out into Common/wfn_utils.cpp by DAS

   based on marching cubes and marching tetrahedra codes by P. Bourke
   http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
   See surface.inp for details on input file.

------------------------------------------------------------------------------*/

// FIXME: classes, c++ standard library

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "wfn_utils.h"

const int MAXBOND = 16;
const double ATOMSIZE = 0.20;
const double BONDRADIUS = 0.15;
const double BONDTOLERANCE = 0.45;
const double MINIMUMBONDINGDISTANCE = 0.40;
const double STRICTVALENCE = false;

typedef struct {
   CARTESIAN point[8];
   CARTESIAN vector[8];
   double scalar[8];
} CUBE;

typedef struct {
   CARTESIAN point[4];
   CARTESIAN vector[4];
   double scalar[4];
} TETRAHEDRON;

typedef struct {
   CARTESIAN point[3];
} TRIANGLE;

struct TRIANGLELIST {
   TRIANGLE triangle;
   TRIANGLELIST *node;
};

typedef struct {
   CARTESIAN vector[3];
} NORMAL;

struct NORMALLIST {
   NORMAL normal;
   NORMALLIST *node;
};

const int edge[12][2] = {
   {0, 1}, {1, 2}, {2, 3}, {3, 0},
   {4, 5}, {5, 6}, {6, 7}, {7, 4},
   {0, 4}, {1, 5}, {2, 6}, {3, 7}};

const int face[12][3] = {
   {0, 1, 2}, {0, 3, 2}, {4, 5, 6}, {4, 7, 6},
   {0, 1, 5}, {0, 4, 5}, {1, 2, 6}, {1, 5, 6},
   {2, 3, 7}, {2, 6, 7}, {3, 0, 4}, {3, 7, 4}};

const int vertex_tetrahedron[6][4] = {
   {0, 5, 1, 6}, {0, 1, 2, 6}, {0, 2, 3, 6},
   {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}};

const int edge_tetrahedron[6][2] = {
   {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

const int edge_table_tetrahedron[16] = {
   0x0000, 0x000d, 0x0013, 0x001e, 0x0026, 0x002b, 0x0035, 0x0038,
   0x0038, 0x0035, 0x002b, 0x0026, 0x001e, 0x0013, 0x000d, 0x0000};

const int triangle_table_tetrahedron[16][7] = {
   {-1, -1, -1, -1, -1, -1, -1},
   { 0,  3,  2, -1, -1, -1, -1},
   { 0,  1,  4, -1, -1, -1, -1},
   { 1,  4,  2,  2,  4,  3, -1},
   { 1,  2,  5, -1, -1, -1, -1},
   { 0,  3,  5,  0,  5,  1, -1},
   { 0,  2,  5,  0,  5,  4, -1},
   { 5,  4,  3, -1, -1, -1, -1},
   { 3,  4,  5, -1, -1, -1, -1},
   { 4,  5,  0,  5,  2,  0, -1},
   { 1,  5,  0,  5,  3,  0, -1},
   { 5,  2,  1, -1, -1, -1, -1},
   { 3,  4,  2,  2,  4,  1, -1},
   { 4,  1,  0, -1, -1, -1, -1},
   { 2,  3,  0, -1, -1, -1, -1},
   {-1, -1, -1, -1, -1, -1, -1}};

const int edge_table_cube[256] = {
   0x0000, 0x0109, 0x0203, 0x030a, 0x0406, 0x050f, 0x0605, 0x070c,
   0x080c, 0x0905, 0x0a0f, 0x0b06, 0x0c0a, 0x0d03, 0x0e09, 0x0f00,
   0x0190, 0x0099, 0x0393, 0x029a, 0x0596, 0x049f, 0x0795, 0x069c,
   0x099c, 0x0895, 0x0b9f, 0x0a96, 0x0d9a, 0x0c93, 0x0f99, 0x0e90,
   0x0230, 0x0339, 0x0033, 0x013a, 0x0636, 0x073f, 0x0435, 0x053c,
   0x0a3c, 0x0b35, 0x083f, 0x0936, 0x0e3a, 0x0f33, 0x0c39, 0x0d30,
   0x03a0, 0x02a9, 0x01a3, 0x00aa, 0x07a6, 0x06af, 0x05a5, 0x04ac,
   0x0bac, 0x0aa5, 0x09af, 0x08a6, 0x0faa, 0x0ea3, 0x0da9, 0x0ca0,
   0x0460, 0x0569, 0x0663, 0x076a, 0x0066, 0x016f, 0x0265, 0x036c,
   0x0c6c, 0x0d65, 0x0e6f, 0x0f66, 0x086a, 0x0963, 0x0a69, 0x0b60,
   0x05f0, 0x04f9, 0x07f3, 0x06fa, 0x01f6, 0x00ff, 0x03f5, 0x02fc,
   0x0dfc, 0x0cf5, 0x0fff, 0x0ef6, 0x09fa, 0x08f3, 0x0bf9, 0x0af0,
   0x0650, 0x0759, 0x0453, 0x055a, 0x0256, 0x035f, 0x0055, 0x015c,
   0x0e5c, 0x0f55, 0x0c5f, 0x0d56, 0x0a5a, 0x0b53, 0x0859, 0x0950,
   0x07c0, 0x06c9, 0x05c3, 0x04ca, 0x03c6, 0x02cf, 0x01c5, 0x00cc,
   0x0fcc, 0x0ec5, 0x0dcf, 0x0cc6, 0x0bca, 0x0ac3, 0x09c9, 0x08c0,
   0x08c0, 0x09c9, 0x0ac3, 0x0bca, 0x0cc6, 0x0dcf, 0x0ec5, 0x0fcc,
   0x00cc, 0x01c5, 0x02cf, 0x03c6, 0x04ca, 0x05c3, 0x06c9, 0x07c0,
   0x0950, 0x0859, 0x0b53, 0x0a5a, 0x0d56, 0x0c5f, 0x0f55, 0x0e5c,
   0x015c, 0x0055, 0x035f, 0x0256, 0x055a, 0x0453, 0x0759, 0x0650,
   0x0af0, 0x0bf9, 0x08f3, 0x09fa, 0x0ef6, 0x0fff, 0x0cf5, 0x0dfc,
   0x02fc, 0x03f5, 0x00ff, 0x01f6, 0x06fa, 0x07f3, 0x04f9, 0x05f0,
   0x0b60, 0x0a69, 0x0963, 0x086a, 0x0f66, 0x0e6f, 0x0d65, 0x0c6c,
   0x036c, 0x0265, 0x016f, 0x0066, 0x076a, 0x0663, 0x0569, 0x0460,
   0x0ca0, 0x0da9, 0x0ea3, 0x0faa, 0x08a6, 0x09af, 0x0aa5, 0x0bac,
   0x04ac, 0x05a5, 0x06af, 0x07a6, 0x00aa, 0x01a3, 0x02a9, 0x03a0,
   0x0d30, 0x0c39, 0x0f33, 0x0e3a, 0x0936, 0x083f, 0x0b35, 0x0a3c,
   0x053c, 0x0435, 0x073f, 0x0636, 0x013a, 0x0033, 0x0339, 0x0230,
   0x0e90, 0x0f99, 0x0c93, 0x0d9a, 0x0a96, 0x0b9f, 0x0895, 0x099c,
   0x069c, 0x0795, 0x049f, 0x0596, 0x029a, 0x0393, 0x0099, 0x0190,
   0x0f00, 0x0e09, 0x0d03, 0x0c0a, 0x0b06, 0x0a0f, 0x0905, 0x080c,
   0x070c, 0x0605, 0x050f, 0x0406, 0x030a, 0x0203, 0x0109, 0x0000};

const int triangle_table_cube[256][16] = {
   {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  8,  3,  9,  8,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  3,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 9,  2, 10,  0,  2,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 2,  8,  3,  2, 10,  8, 10,  9,  8, -1, -1, -1, -1, -1, -1, -1},
   { 3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0, 11,  2,  8, 11,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  9,  0,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1, 11,  2,  1,  9, 11,  9,  8, 11, -1, -1, -1, -1, -1, -1, -1},
   { 3, 10,  1, 11, 10,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0, 10,  1,  0,  8, 10,  8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
   { 3,  9,  0,  3, 11,  9, 11, 10,  9, -1, -1, -1, -1, -1, -1, -1},
   { 9,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  3,  0,  7,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  1,  9,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  1,  9,  4,  7,  1,  7,  3,  1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 10,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 3,  4,  7,  3,  0,  4,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1},
   { 9,  2, 10,  9,  0,  2,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1},
   { 2, 10,  9,  2,  9,  7,  2,  7,  3,  7,  9,  4, -1, -1, -1, -1},
   { 8,  4,  7,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11,  4,  7, 11,  2,  4,  2,  0,  4, -1, -1, -1, -1, -1, -1, -1},
   { 9,  0,  1,  8,  4,  7,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1},
   { 4,  7, 11,  9,  4, 11,  9, 11,  2,  9,  2,  1, -1, -1, -1, -1},
   { 3, 10,  1,  3, 11, 10,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1},
   { 1, 11, 10,  1,  4, 11,  1,  0,  4,  7, 11,  4, -1, -1, -1, -1},
   { 4,  7,  8,  9,  0, 11,  9, 11, 10, 11,  0,  3, -1, -1, -1, -1},
   { 4,  7, 11,  4, 11,  9,  9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
   { 9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 9,  5,  4,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  5,  4,  1,  5,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 8,  5,  4,  8,  3,  5,  3,  1,  5, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 10,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 3,  0,  8,  1,  2, 10,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1},
   { 5,  2, 10,  5,  4,  2,  4,  0,  2, -1, -1, -1, -1, -1, -1, -1},
   { 2, 10,  5,  3,  2,  5,  3,  5,  4,  3,  4,  8, -1, -1, -1, -1},
   { 9,  5,  4,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0, 11,  2,  0,  8, 11,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1},
   { 0,  5,  4,  0,  1,  5,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1},
   { 2,  1,  5,  2,  5,  8,  2,  8, 11,  4,  8,  5, -1, -1, -1, -1},
   {10,  3, 11, 10,  1,  3,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1},
   { 4,  9,  5,  0,  8,  1,  8, 10,  1,  8, 11, 10, -1, -1, -1, -1},
   { 5,  4,  0,  5,  0, 11,  5, 11, 10, 11,  0,  3, -1, -1, -1, -1},
   { 5,  4,  8,  5,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1},
   { 9,  7,  8,  5,  7,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 9,  3,  0,  9,  5,  3,  5,  7,  3, -1, -1, -1, -1, -1, -1, -1},
   { 0,  7,  8,  0,  1,  7,  1,  5,  7, -1, -1, -1, -1, -1, -1, -1},
   { 1,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 9,  7,  8,  9,  5,  7, 10,  1,  2, -1, -1, -1, -1, -1, -1, -1},
   {10,  1,  2,  9,  5,  0,  5,  3,  0,  5,  7,  3, -1, -1, -1, -1},
   { 8,  0,  2,  8,  2,  5,  8,  5,  7, 10,  5,  2, -1, -1, -1, -1},
   { 2, 10,  5,  2,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1},
   { 7,  9,  5,  7,  8,  9,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1},
   { 9,  5,  7,  9,  7,  2,  9,  2,  0,  2,  7, 11, -1, -1, -1, -1},
   { 2,  3, 11,  0,  1,  8,  1,  7,  8,  1,  5,  7, -1, -1, -1, -1},
   {11,  2,  1, 11,  1,  7,  7,  1,  5, -1, -1, -1, -1, -1, -1, -1},
   { 9,  5,  8,  8,  5,  7, 10,  1,  3, 10,  3, 11, -1, -1, -1, -1},
   { 5,  7,  0,  5,  0,  9,  7, 11,  0,  1,  0, 10, 11, 10,  0, -1},
   {11, 10,  0, 11,  0,  3, 10,  5,  0,  8,  0,  7,  5,  7,  0, -1},
   {11, 10,  5,  7, 11,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  3,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 9,  0,  1,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  8,  3,  1,  9,  8,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1},
   { 1,  6,  5,  2,  6,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  6,  5,  1,  2,  6,  3,  0,  8, -1, -1, -1, -1, -1, -1, -1},
   { 9,  6,  5,  9,  0,  6,  0,  2,  6, -1, -1, -1, -1, -1, -1, -1},
   { 5,  9,  8,  5,  8,  2,  5,  2,  6,  3,  2,  8, -1, -1, -1, -1},
   { 2,  3, 11, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11,  0,  8, 11,  2,  0, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1},
   { 0,  1,  9,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1},
   { 5, 10,  6,  1,  9,  2,  9, 11,  2,  9,  8, 11, -1, -1, -1, -1},
   { 6,  3, 11,  6,  5,  3,  5,  1,  3, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8, 11,  0, 11,  5,  0,  5,  1,  5, 11,  6, -1, -1, -1, -1},
   { 3, 11,  6,  0,  3,  6,  0,  6,  5,  0,  5,  9, -1, -1, -1, -1},
   { 6,  5,  9,  6,  9, 11, 11,  9,  8, -1, -1, -1, -1, -1, -1, -1},
   { 5, 10,  6,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  3,  0,  4,  7,  3,  6,  5, 10, -1, -1, -1, -1, -1, -1, -1},
   { 1,  9,  0,  5, 10,  6,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1},
   {10,  6,  5,  1,  9,  7,  1,  7,  3,  7,  9,  4, -1, -1, -1, -1},
   { 6,  1,  2,  6,  5,  1,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2,  5,  5,  2,  6,  3,  0,  4,  3,  4,  7, -1, -1, -1, -1},
   { 8,  4,  7,  9,  0,  5,  0,  6,  5,  0,  2,  6, -1, -1, -1, -1},
   { 7,  3,  9,  7,  9,  4,  3,  2,  9,  5,  9,  6,  2,  6,  9, -1},
   { 3, 11,  2,  7,  8,  4, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1},
   { 5, 10,  6,  4,  7,  2,  4,  2,  0,  2,  7, 11, -1, -1, -1, -1},
   { 0,  1,  9,  4,  7,  8,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1},
   { 9,  2,  1,  9, 11,  2,  9,  4, 11,  7, 11,  4,  5, 10,  6, -1},
   { 8,  4,  7,  3, 11,  5,  3,  5,  1,  5, 11,  6, -1, -1, -1, -1},
   { 5,  1, 11,  5, 11,  6,  1,  0, 11,  7, 11,  4,  0,  4, 11, -1},
   { 0,  5,  9,  0,  6,  5,  0,  3,  6, 11,  6,  3,  8,  4,  7, -1},
   { 6,  5,  9,  6,  9, 11,  4,  7,  9,  7, 11,  9, -1, -1, -1, -1},
   {10,  4,  9,  6,  4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4, 10,  6,  4,  9, 10,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1},
   {10,  0,  1, 10,  6,  0,  6,  4,  0, -1, -1, -1, -1, -1, -1, -1},
   { 8,  3,  1,  8,  1,  6,  8,  6,  4,  6,  1, 10, -1, -1, -1, -1},
   { 1,  4,  9,  1,  2,  4,  2,  6,  4, -1, -1, -1, -1, -1, -1, -1},
   { 3,  0,  8,  1,  2,  9,  2,  4,  9,  2,  6,  4, -1, -1, -1, -1},
   { 0,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 8,  3,  2,  8,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1},
   {10,  4,  9, 10,  6,  4, 11,  2,  3, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  2,  2,  8, 11,  4,  9, 10,  4, 10,  6, -1, -1, -1, -1},
   { 3, 11,  2,  0,  1,  6,  0,  6,  4,  6,  1, 10, -1, -1, -1, -1},
   { 6,  4,  1,  6,  1, 10,  4,  8,  1,  2,  1, 11,  8, 11,  1, -1},
   { 9,  6,  4,  9,  3,  6,  9,  1,  3, 11,  6,  3, -1, -1, -1, -1},
   { 8, 11,  1,  8,  1,  0, 11,  6,  1,  9,  1,  4,  6,  4,  1, -1},
   { 3, 11,  6,  3,  6,  0,  0,  6,  4, -1, -1, -1, -1, -1, -1, -1},
   { 6,  4,  8, 11,  6,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 7, 10,  6,  7,  8, 10,  8,  9, 10, -1, -1, -1, -1, -1, -1, -1},
   { 0,  7,  3,  0, 10,  7,  0,  9, 10,  6,  7, 10, -1, -1, -1, -1},
   {10,  6,  7,  1, 10,  7,  1,  7,  8,  1,  8,  0, -1, -1, -1, -1},
   {10,  6,  7, 10,  7,  1,  1,  7,  3, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2,  6,  1,  6,  8,  1,  8,  9,  8,  6,  7, -1, -1, -1, -1},
   { 2,  6,  9,  2,  9,  1,  6,  7,  9,  0,  9,  3,  7,  3,  9, -1},
   { 7,  8,  0,  7,  0,  6,  6,  0,  2, -1, -1, -1, -1, -1, -1, -1},
   { 7,  3,  2,  6,  7,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 2,  3, 11, 10,  6,  8, 10,  8,  9,  8,  6,  7, -1, -1, -1, -1},
   { 2,  0,  7,  2,  7, 11,  0,  9,  7,  6,  7, 10,  9, 10,  7, -1},
   { 1,  8,  0,  1,  7,  8,  1, 10,  7,  6,  7, 10,  2,  3, 11, -1},
   {11,  2,  1, 11,  1,  7, 10,  6,  1,  6,  7,  1, -1, -1, -1, -1},
   { 8,  9,  6,  8,  6,  7,  9,  1,  6, 11,  6,  3,  1,  3,  6, -1},
   { 0,  9,  1, 11,  6,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 7,  8,  0,  7,  0,  6,  3, 11,  0, 11,  6,  0, -1, -1, -1, -1},
   { 7, 11,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 3,  0,  8, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  1,  9, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 8,  1,  9,  8,  3,  1, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1},
   {10,  1,  2,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 10,  3,  0,  8,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1},
   { 2,  9,  0,  2, 10,  9,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1},
   { 6, 11,  7,  2, 10,  3, 10,  8,  3, 10,  9,  8, -1, -1, -1, -1},
   { 7,  2,  3,  6,  2,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 7,  0,  8,  7,  6,  0,  6,  2,  0, -1, -1, -1, -1, -1, -1, -1},
   { 2,  7,  6,  2,  3,  7,  0,  1,  9, -1, -1, -1, -1, -1, -1, -1},
   { 1,  6,  2,  1,  8,  6,  1,  9,  8,  8,  7,  6, -1, -1, -1, -1},
   {10,  7,  6, 10,  1,  7,  1,  3,  7, -1, -1, -1, -1, -1, -1, -1},
   {10,  7,  6,  1,  7, 10,  1,  8,  7,  1,  0,  8, -1, -1, -1, -1},
   { 0,  3,  7,  0,  7, 10,  0, 10,  9,  6, 10,  7, -1, -1, -1, -1},
   { 7,  6, 10,  7, 10,  8,  8, 10,  9, -1, -1, -1, -1, -1, -1, -1},
   { 6,  8,  4, 11,  8,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 3,  6, 11,  3,  0,  6,  0,  4,  6, -1, -1, -1, -1, -1, -1, -1},
   { 8,  6, 11,  8,  4,  6,  9,  0,  1, -1, -1, -1, -1, -1, -1, -1},
   { 9,  4,  6,  9,  6,  3,  9,  3,  1, 11,  3,  6, -1, -1, -1, -1},
   { 6,  8,  4,  6, 11,  8,  2, 10,  1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 10,  3,  0, 11,  0,  6, 11,  0,  4,  6, -1, -1, -1, -1},
   { 4, 11,  8,  4,  6, 11,  0,  2,  9,  2, 10,  9, -1, -1, -1, -1},
   {10,  9,  3, 10,  3,  2,  9,  4,  3, 11,  3,  6,  4,  6,  3, -1},
   { 8,  2,  3,  8,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1},
   { 0,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  9,  0,  2,  3,  4,  2,  4,  6,  4,  3,  8, -1, -1, -1, -1},
   { 1,  9,  4,  1,  4,  2,  2,  4,  6, -1, -1, -1, -1, -1, -1, -1},
   { 8,  1,  3,  8,  6,  1,  8,  4,  6,  6, 10,  1, -1, -1, -1, -1},
   {10,  1,  0, 10,  0,  6,  6,  0,  4, -1, -1, -1, -1, -1, -1, -1},
   { 4,  6,  3,  4,  3,  8,  6, 10,  3,  0,  3,  9, 10,  9,  3, -1},
   {10,  9,  4,  6, 10,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  9,  5,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  3,  4,  9,  5, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1},
   { 5,  0,  1,  5,  4,  0,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1},
   {11,  7,  6,  8,  3,  4,  3,  5,  4,  3,  1,  5, -1, -1, -1, -1},
   { 9,  5,  4, 10,  1,  2,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1},
   { 6, 11,  7,  1,  2, 10,  0,  8,  3,  4,  9,  5, -1, -1, -1, -1},
   { 7,  6, 11,  5,  4, 10,  4,  2, 10,  4,  0,  2, -1, -1, -1, -1},
   { 3,  4,  8,  3,  5,  4,  3,  2,  5, 10,  5,  2, 11,  7,  6, -1},
   { 7,  2,  3,  7,  6,  2,  5,  4,  9, -1, -1, -1, -1, -1, -1, -1},
   { 9,  5,  4,  0,  8,  6,  0,  6,  2,  6,  8,  7, -1, -1, -1, -1},
   { 3,  6,  2,  3,  7,  6,  1,  5,  0,  5,  4,  0, -1, -1, -1, -1},
   { 6,  2,  8,  6,  8,  7,  2,  1,  8,  4,  8,  5,  1,  5,  8, -1},
   { 9,  5,  4, 10,  1,  6,  1,  7,  6,  1,  3,  7, -1, -1, -1, -1},
   { 1,  6, 10,  1,  7,  6,  1,  0,  7,  8,  7,  0,  9,  5,  4, -1},
   { 4,  0, 10,  4, 10,  5,  0,  3, 10,  6, 10,  7,  3,  7, 10, -1},
   { 7,  6, 10,  7, 10,  8,  5,  4, 10,  4,  8, 10, -1, -1, -1, -1},
   { 6,  9,  5,  6, 11,  9, 11,  8,  9, -1, -1, -1, -1, -1, -1, -1},
   { 3,  6, 11,  0,  6,  3,  0,  5,  6,  0,  9,  5, -1, -1, -1, -1},
   { 0, 11,  8,  0,  5, 11,  0,  1,  5,  5,  6, 11, -1, -1, -1, -1},
   { 6, 11,  3,  6,  3,  5,  5,  3,  1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 10,  9,  5, 11,  9, 11,  8, 11,  5,  6, -1, -1, -1, -1},
   { 0, 11,  3,  0,  6, 11,  0,  9,  6,  5,  6,  9,  1,  2, 10, -1},
   {11,  8,  5, 11,  5,  6,  8,  0,  5, 10,  5,  2,  0,  2,  5, -1},
   { 6, 11,  3,  6,  3,  5,  2, 10,  3, 10,  5,  3, -1, -1, -1, -1},
   { 5,  8,  9,  5,  2,  8,  5,  6,  2,  3,  8,  2, -1, -1, -1, -1},
   { 9,  5,  6,  9,  6,  0,  0,  6,  2, -1, -1, -1, -1, -1, -1, -1},
   { 1,  5,  8,  1,  8,  0,  5,  6,  8,  3,  8,  2,  6,  2,  8, -1},
   { 1,  5,  6,  2,  1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  3,  6,  1,  6, 10,  3,  8,  6,  5,  6,  9,  8,  9,  6, -1},
   {10,  1,  0, 10,  0,  6,  9,  5,  0,  5,  6,  0, -1, -1, -1, -1},
   { 0,  3,  8,  5,  6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {10,  5,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11,  5, 10,  7,  5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11,  5, 10, 11,  7,  5,  8,  3,  0, -1, -1, -1, -1, -1, -1, -1},
   { 5, 11,  7,  5, 10, 11,  1,  9,  0, -1, -1, -1, -1, -1, -1, -1},
   {10,  7,  5, 10, 11,  7,  9,  8,  1,  8,  3,  1, -1, -1, -1, -1},
   {11,  1,  2, 11,  7,  1,  7,  5,  1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  3,  1,  2,  7,  1,  7,  5,  7,  2, 11, -1, -1, -1, -1},
   { 9,  7,  5,  9,  2,  7,  9,  0,  2,  2, 11,  7, -1, -1, -1, -1},
   { 7,  5,  2,  7,  2, 11,  5,  9,  2,  3,  2,  8,  9,  8,  2, -1},
   { 2,  5, 10,  2,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1},
   { 8,  2,  0,  8,  5,  2,  8,  7,  5, 10,  2,  5, -1, -1, -1, -1},
   { 9,  0,  1,  5, 10,  3,  5,  3,  7,  3, 10,  2, -1, -1, -1, -1},
   { 9,  8,  2,  9,  2,  1,  8,  7,  2, 10,  2,  5,  7,  5,  2, -1},
   { 1,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  7,  0,  7,  1,  1,  7,  5, -1, -1, -1, -1, -1, -1, -1},
   { 9,  0,  3,  9,  3,  5,  5,  3,  7, -1, -1, -1, -1, -1, -1, -1},
   { 9,  8,  7,  5,  9,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 5,  8,  4,  5, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1},
   { 5,  0,  4,  5, 11,  0,  5, 10, 11, 11,  3,  0, -1, -1, -1, -1},
   { 0,  1,  9,  8,  4, 10,  8, 10, 11, 10,  4,  5, -1, -1, -1, -1},
   {10, 11,  4, 10,  4,  5, 11,  3,  4,  9,  4,  1,  3,  1,  4, -1},
   { 2,  5,  1,  2,  8,  5,  2, 11,  8,  4,  5,  8, -1, -1, -1, -1},
   { 0,  4, 11,  0, 11,  3,  4,  5, 11,  2, 11,  1,  5,  1, 11, -1},
   { 0,  2,  5,  0,  5,  9,  2, 11,  5,  4,  5,  8, 11,  8,  5, -1},
   { 9,  4,  5,  2, 11,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 2,  5, 10,  3,  5,  2,  3,  4,  5,  3,  8,  4, -1, -1, -1, -1},
   { 5, 10,  2,  5,  2,  4,  4,  2,  0, -1, -1, -1, -1, -1, -1, -1},
   { 3, 10,  2,  3,  5, 10,  3,  8,  5,  4,  5,  8,  0,  1,  9, -1},
   { 5, 10,  2,  5,  2,  4,  1,  9,  2,  9,  4,  2, -1, -1, -1, -1},
   { 8,  4,  5,  8,  5,  3,  3,  5,  1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  4,  5,  1,  0,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 8,  4,  5,  8,  5,  3,  9,  0,  5,  0,  3,  5, -1, -1, -1, -1},
   { 9,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4, 11,  7,  4,  9, 11,  9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
   { 0,  8,  3,  4,  9,  7,  9, 11,  7,  9, 10, 11, -1, -1, -1, -1},
   { 1, 10, 11,  1, 11,  4,  1,  4,  0,  7,  4, 11, -1, -1, -1, -1},
   { 3,  1,  4,  3,  4,  8,  1, 10,  4,  7,  4, 11, 10, 11,  4, -1},
   { 4, 11,  7,  9, 11,  4,  9,  2, 11,  9,  1,  2, -1, -1, -1, -1},
   { 9,  7,  4,  9, 11,  7,  9,  1, 11,  2, 11,  1,  0,  8,  3, -1},
   {11,  7,  4, 11,  4,  2,  2,  4,  0, -1, -1, -1, -1, -1, -1, -1},
   {11,  7,  4, 11,  4,  2,  8,  3,  4,  3,  2,  4, -1, -1, -1, -1},
   { 2,  9, 10,  2,  7,  9,  2,  3,  7,  7,  4,  9, -1, -1, -1, -1},
   { 9, 10,  7,  9,  7,  4, 10,  2,  7,  8,  7,  0,  2,  0,  7, -1},
   { 3,  7, 10,  3, 10,  2,  7,  4, 10,  1, 10,  0,  4,  0, 10, -1},
   { 1, 10,  2,  8,  7,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  9,  1,  4,  1,  7,  7,  1,  3, -1, -1, -1, -1, -1, -1, -1},
   { 4,  9,  1,  4,  1,  7,  0,  8,  1,  8,  7,  1, -1, -1, -1, -1},
   { 4,  0,  3,  7,  4,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 4,  8,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 9, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 3,  0,  9,  3,  9, 11, 11,  9, 10, -1, -1, -1, -1, -1, -1, -1},
   { 0,  1, 10,  0, 10,  8,  8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
   { 3,  1, 10, 11,  3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  2, 11,  1, 11,  9,  9, 11,  8, -1, -1, -1, -1, -1, -1, -1},
   { 3,  0,  9,  3,  9, 11,  1,  2,  9,  2, 11,  9, -1, -1, -1, -1},
   { 0,  2, 11,  8,  0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 3,  2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 2,  3,  8,  2,  8, 10, 10,  8,  9, -1, -1, -1, -1, -1, -1, -1},
   { 9, 10,  2,  0,  9,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 2,  3,  8,  2,  8, 10,  0,  1,  8,  1, 10,  8, -1, -1, -1, -1},
   { 1, 10,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 1,  3,  8,  9,  1,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  9,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   { 0,  3,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

int *uas = NULL;
CARTESIAN *uap = NULL;
int *cs = NULL;
int *nb = NULL;
CARTESIAN **bonds = NULL;
int **sci = NULL;
CARTESIAN ***vector = NULL;
TRIANGLELIST *triangle[2] = {NULL, NULL};
NORMALLIST *normal[2] = {NULL, NULL};

int par_read(char *pfn, char *header, char *ifn, char *iff, char *ofn, char *off, double *isovalue, int *sign, int *power, int *algorithm, bool *smooth, bool *box, bool *basis, bool *uc, CARTESIAN *uco, CARTESIAN *ucv, int *ucf, bool *sc, CARTESIAN *sco, CARTESIAN *scv, int *scf, bool *sct)
{
   int ii, ierr, icount = 0, icheck = 0;
   char s1[MAXCHAR], s2[MAXCHAR], s3[MAXCHAR];
   char* trash;
   FILE *hh;

   hh = fopen(pfn, "r");
   if (hh == NULL)
      return -1;

   while (!feof(hh))
   {
      strncpy(s1, "\0", 1);
      trash = fgets(s1, MAXCHAR, hh);
      parse(s1, s2, s3);
      lowercase(s2);
      if (strcmp(s2, "header") == 0)
      {
         icount++;
         strncpy(header, s3, MAXCHAR);
         header[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "inputfilename") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(ifn, s3, MAXCHAR);
         ifn[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "inputfileformat") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(iff, s3, MAXCHAR);
         iff[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "outputfilename") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(ofn, s3, MAXCHAR);
         ofn[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "outputfileformat") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(off, s3, MAXCHAR);
         off[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "isovalue") == 0)
      {
         icount++;
         *isovalue = atof(s3);
         if (*isovalue < EPS9 || *isovalue > 1.0 - EPS9)
	 {
	    fprintf(stderr, "Value for 'isovalue' is out of bounds.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "sign") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "positive") == 0)
            *sign = 0;
         else if (strcmp(s3, "negative") == 0)
            *sign = 1;
         else if (strcmp(s3, "both") == 0)
            *sign = 2;
         else
	 {
	    fprintf(stderr, "Value for 'sign' is not legal.\n");
	    icheck--;
	 }
      }
      else if (strcmp(s2, "power") == 0)
      {
         icount++;
         *power = atoi(s3);
         if (*power < 0 || *power > 2)
	 {
	    fprintf(stderr, "Value for 'power' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "algorithm") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "cube") == 0)
            *algorithm = 0;
         else if (strcmp(s3, "tetrahedron") == 0)
            *algorithm = 1;
         else
	 {
	    fprintf(stderr, "Value for 'algorithm' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "smooth") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            *smooth = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            *smooth = true;
         else
	 {
	    fprintf(stderr, "Value for 'smooth' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "box") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            *box = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            *box = true;
         else
	 {
	    fprintf(stderr, "Value for 'box' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "basis") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            *basis = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            *basis = true;
         else
	 {
	    fprintf(stderr, "Value for 'basis' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "uc") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            *uc = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            *uc = true;
         else
	 {
	    fprintf(stderr, "Value for 'uc' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "uco") == 0)
      {
         icount++;
         ierr = fscanf(hh, "%le%le%le\n", &uco->xx, &uco->yy, &uco->zz);
      }
      else if (strcmp(s2, "ucv") == 0)
      {
         icount++;
         for (ii = 0; ii < 3; ii++)
            ierr = fscanf(hh, "%le%le%le\n", &ucv[ii].xx, &ucv[ii].yy, &ucv[ii].zz);
      }
      else if (strcmp(s2, "ucu") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "bohr") == 0)
            *ucf = 0;
         else if (strcmp(s3, "angstrom") == 0)
            *ucf = 1;
         else if (strcmp(s3, "latvec") == 0)
            *ucf = 2;
         else
	 {
	    fprintf(stderr, "Value for 'ucu' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "sc") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            *sc = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            *sc = true;
         else
	 {
	    fprintf(stderr, "Value for 'sc' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "sco") == 0)
      {
         icount++;
         ierr = fscanf(hh, "%le%le%le\n", &sco->xx, &sco->yy, &sco->zz);
      }
      else if (strcmp(s2, "scv") == 0)
      {
         icount++;
         for (ii = 0; ii < 3; ii++)
            ierr = fscanf(hh, "%le%le%le\n", &scv[ii].xx, &scv[ii].yy, &scv[ii].zz);
      }
      else if (strcmp(s2, "scu") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "bohr") == 0)
            *scf = 0;
         else if (strcmp(s3, "angstrom") == 0)
            *scf = 1;
         else if (strcmp(s3, "latvec") == 0)
            *scf = 2;
         else
	 {
	    fprintf(stderr, "Value for 'scu' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "sct") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            *sct = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            *sct = true;
         else
	 {
	    fprintf(stderr, "Value for 'sct' is not legal.\n");
            icheck--;
	 }
      }
   }

   ierr = fclose(hh);
   if (ierr != 0)
      return -1;

   if(icount != 21)
   {
      fprintf(stderr, "Not all parameters were set.\n");
   }

   if (icount != 21 || icheck < 0)
      return -1;

   return 0;
}

int atom_clone(int na, CARTESIAN *ucv, CARTESIAN sco, CARTESIAN *scv, int *scmin, int *scmax)
{
   bool flag;
   int ii, jj, kk, ll, mm, nn, ss, tt, apmin[3], apmax[3];
   double scvw[3], ucvw[3], pp, psccenter[3], pscmin[3], pscmax[3], pa[3];
   CARTESIAN scvn[3], ucvn[3], sccenter, sccorner[8], sa;

   normal_make(scv, scvn, scvw);
   normal_make(ucv, ucvn, ucvw);

   sccenter.xx = sco.xx + scv[0].xx * 0.5 + scv[1].xx * 0.5 + scv[2].xx * 0.5;
   sccenter.yy = sco.yy + scv[0].yy * 0.5 + scv[1].yy * 0.5 + scv[2].yy * 0.5;
   sccenter.zz = sco.zz + scv[0].zz * 0.5 + scv[1].zz * 0.5 + scv[2].zz * 0.5;

   for (ii = 0; ii < 8; ii++)
   {
      sccorner[ii].xx = sccenter.xx;
      sccorner[ii].yy = sccenter.yy;
      sccorner[ii].zz = sccenter.zz;
      for (jj = 0; jj < 3; jj++)
      {
         sccorner[ii].xx += scv[jj].xx * (double(vertex[ii][jj]) - 0.5);
         sccorner[ii].yy += scv[jj].yy * (double(vertex[ii][jj]) - 0.5);
         sccorner[ii].zz += scv[jj].zz * (double(vertex[ii][jj]) - 0.5);
      }
   }

   for (ii = 0; ii < 3; ii++)
   {
      scmin[ii] = INT_MAX;
      scmax[ii] = INT_MIN;
      for (jj = 0; jj < 8; jj++)
      {
         pp = (sccorner[jj].xx * ucvn[ii].xx + sccorner[jj].yy * ucvn[ii].yy + sccorner[jj].zz * ucvn[ii].zz) / ucvw[ii];
         // this is similar to double_to_int but we want to convert
         // the interval (-2,-1] to -2, (-1,0] to -1, (0,1] to 1,
         // (1,2] to 2, that is, skipping the integer value of 0.
         if (pp < 0.0)
            kk = int(pp - 1.0);
         else
            kk = int(pp + 1.0);
         if (kk < scmin[ii])
            scmin[ii] = kk;
         if (kk > scmax[ii])
            scmax[ii] = kk;
      }
   }

   for (ii = 0; ii < 3; ii++)
   {
      apmin[ii] = INT_MAX;
      apmax[ii] = INT_MIN;
      for (jj = 0; jj < na; jj++)
      {
         pp = (ap[jj].xx * ucvn[ii].xx + ap[jj].yy * ucvn[ii].yy + ap[jj].zz * ucvn[ii].zz) / ucvw[ii];
         // this is similar to double_to_int but we want to convert
         // the interval (-2,-1] to -2, (-1,0] to -1, (0,1] to 1,
         // (1,2] to 2, that is, skipping the integer value of 0.
         if (pp < 0.0)
            kk = int(pp - 1.0);
         else
            kk = int(pp + 1.0);
         if (kk < apmin[ii])
            apmin[ii] = kk;
         if (kk > apmax[ii])
            apmax[ii] = kk;
      }
   }

   jj = 0;
   for (ii = 0; ii < 3; ii++)
   {
      if (abs(apmin[ii]) > jj)
         jj = abs(apmin[ii]);
      if (abs(apmax[ii]) > jj)
         jj = abs(apmax[ii]);
   }
   for (ii = 0; ii < 3; ii++)
   {
      scmin[ii] -= jj;
      scmax[ii] += jj;
   }

   for (ii = 0; ii < 3; ii++)
      psccenter[ii] = sccenter.xx * scvn[ii].xx + sccenter.yy * scvn[ii].yy + sccenter.zz * scvn[ii].zz;
   for (ii = 0; ii < 3; ii++)
      pscmin[ii] = psccenter[ii] - 0.5 * scvw[ii] - EPS9;
   for (ii = 0; ii < 3; ii++)
      pscmax[ii] = psccenter[ii] + 0.5 * scvw[ii] - EPS9;

   // loop over tt is used to construct a supercell twice
   // tt = 0 is a dry run intended to identify number of atoms in the supercell nn
   // tt = 1 is for allocating and filling arrays of atomic species and positions
   uas = as;
   uap = ap;
   nn = 0;
   ss = 0;
   for (tt = 0; tt < 2; tt++)
   {
      if (tt != 0)
      {
         as = new int[MAX(1,nn)];
         ap = new CARTESIAN[MAX(1,nn)];
         sci = new int*[MAX(1,nn)];
         for (ii = 0; ii < nn; ii++)
            sci[ii] = new int[4];
      }
      for (ii = scmin[0]; ii <= scmax[0]; ii++)
      for (jj = scmin[1]; jj <= scmax[1]; jj++)
      for (kk = scmin[2]; kk <= scmax[2]; kk++)
      for (ll = 0; ll < na; ll++)
      {
         sa.xx = ucv[0].xx * double(ii) + ucv[1].xx * double(jj) + ucv[2].xx * double(kk) + uap[ll].xx;
         sa.yy = ucv[0].yy * double(ii) + ucv[1].yy * double(jj) + ucv[2].yy * double(kk) + uap[ll].yy;
         sa.zz = ucv[0].zz * double(ii) + ucv[1].zz * double(jj) + ucv[2].zz * double(kk) + uap[ll].zz;
         for (mm = 0; mm < 3; mm++)
            pa[mm] = sa.xx * scvn[mm].xx + sa.yy * scvn[mm].yy + sa.zz * scvn[mm].zz;
         flag = true;
         for (mm = 0; mm < 3; mm++)
            if (pa[mm] < pscmin[mm] || pa[mm] > pscmax[mm])
               flag = false;
         if (flag)
         {
            if (tt == 0)
               nn++;
            else
            {
               as[ss] = uas[ll];
               ap[ss].xx = sa.xx;
               ap[ss].yy = sa.yy;
               ap[ss].zz = sa.zz;
               sci[ss][0] = ll;
               sci[ss][1] = ii;
               sci[ss][2] = jj;
               sci[ss][3] = kk;
               ss++;
            }
         }
      }
   }

   return nn;
}

void scalar_normal(int ni, int nj, int nk, const scalar_field & sf)
{
   int ii, jj, kk, i1, j1, k1, i2, j2, k2;
   double weight, gradient[3];
   CARTESIAN normal;

   vector = new CARTESIAN**[ni];
   for (ii = 0; ii < ni; ii++)
      vector[ii] = new CARTESIAN*[nj];
   for (ii = 0; ii < ni; ii++)
   for (jj = 0; jj < nj; jj++)
      vector[ii][jj] = new CARTESIAN[nk];

   for (ii = 0; ii < ni; ii++)
   {
      if (ii == 0)
         i1 = ii;
      else
         i1 = ii - 1;
      if (ii == ni - 1)
         i2 = ii;
      else
         i2 = ii + 1;
      for (jj = 0; jj < nj; jj++)
      {
         if (jj == 0)
            j1 = jj;
         else
            j1 = jj - 1;
         if (jj == nj - 1)
            j2 = jj;
         else
            j2 = jj + 1;
         for (kk = 0; kk < nk; kk++)
         {
            if (kk == 0)
               k1 = kk;
            else
               k1 = kk - 1;
            if (kk == nk - 1)
               k2 = kk;
            else
               k2 = kk + 1;
            gradient[0] = sf.scalar[i2][jj][kk] - sf.scalar[i1][jj][kk];
            gradient[1] = sf.scalar[ii][j2][kk] - sf.scalar[ii][j1][kk];
            gradient[2] = sf.scalar[ii][jj][k2] - sf.scalar[ii][jj][k1];
            normal.xx = gradient[0] * sf.stepv[0].xx + gradient[1] * sf.stepv[1].xx + gradient[2] * sf.stepv[2].xx;
            normal.yy = gradient[0] * sf.stepv[0].yy + gradient[1] * sf.stepv[1].yy + gradient[2] * sf.stepv[2].yy;
            normal.zz = gradient[0] * sf.stepv[0].zz + gradient[1] * sf.stepv[1].zz + gradient[2] * sf.stepv[2].zz;
            weight = sqrt(normal.xx * normal.xx + normal.yy * normal.yy + normal.zz * normal.zz);
            if (weight < EPS9)
               weight = 1.0;
            normal.xx = normal.xx / weight;
            normal.yy = normal.yy / weight;
            normal.zz = normal.zz / weight;
            vector[ii][jj][kk] = normal;
         }
      }
   }
}

int chemical_species(int na)
{
   bool flag;
   int ii, jj, ns = 0;
   int *species;

   species = new int[MAX(1,na)];
   for (ii = 0; ii < na; ii++)
      species[ii] = 0;
   for (ii = 0; ii < na; ii++)
   {
      flag = true;
      for (jj = 0; jj < ns; jj++)
         if (species[jj] == as[ii])
            flag = false;
      if (flag)
      {
         species[ns] = as[ii];
         ns++;
      }
   }

   cs = new int[MAX(1,ns)];
   for (ii = 0; ii < ns; ii++)
      cs[ii] = species[ii];
   delete [] species;

   return ns;
}

int bonds_st(int na, CARTESIAN *scv)
{
   int ii, i1, i2, j1, j2, j3, v1;
   double rr, r1, r2;
   CARTESIAN pp;

   nb = new int[MAX(1,na)];
   bonds = new CARTESIAN*[MAX(1,na)];
   for (ii = 0; ii < na; ii++)
      bonds[ii] = new CARTESIAN[MAXBOND];

   for (i1 = 0; i1 < na; i1++)
   {
      v1 = periodic_table[as[i1]].valence;
      r1 = periodic_table[as[i1]].rcov;
      ii = 0;
      for (i2 = 0; i2 < na; i2++)
      {
         r2 = periodic_table[as[i2]].rcov;
         for (j1 = -1; j1 < 2; j1++)
         for (j2 = -1; j2 < 2; j2++)
         for (j3 = -1; j3 < 2; j3++)
         if ((j3 != 0 || j2 != 0 || j1 != 0 || i2 != i1) && as[i2] != 0 && as[i1] != 0)
         {
            pp.xx = scv[2].xx * double(j3) + scv[1].xx * double(j2) + scv[0].xx * double(j1) + ap[i2].xx - ap[i1].xx;
            pp.yy = scv[2].yy * double(j3) + scv[1].yy * double(j2) + scv[0].yy * double(j1) + ap[i2].yy - ap[i1].yy;
            pp.zz = scv[2].zz * double(j3) + scv[1].zz * double(j2) + scv[0].zz * double(j1) + ap[i2].zz - ap[i1].zz;
            rr = sqrt(pp.xx * pp.xx + pp.yy * pp.yy + pp.zz * pp.zz);
            if (rr < MINIMUMBONDINGDISTANCE)
               return -1;
            if (rr < r1 + r2 + BONDTOLERANCE)
            {
               if (ii >= MAXBOND)
                  return -1;
               if (STRICTVALENCE && ii >= v1)
                  return -1;
               bonds[i1][ii].xx = ap[i1].xx + pp.xx * r1 / (r1 + r2);
               bonds[i1][ii].yy = ap[i1].yy + pp.yy * r1 / (r1 + r2);
               bonds[i1][ii].zz = ap[i1].zz + pp.zz * r1 / (r1 + r2);
               ii++;
            }
         }
      }
      nb[i1] = ii;
   }

   return 0;
}

int bonds_ut(int na, int una, CARTESIAN *ucv, int *scmin, int *scmax)
{
   int ii, i1, i2, j1, j2, j3, v1;
   double rr, r1, r2;
   CARTESIAN pp;

   nb = new int[MAX(1,na)];
   bonds = new CARTESIAN*[MAX(1,na)];
   for (ii = 0; ii < na; ii++)
      bonds[ii] = new CARTESIAN[MAXBOND];

   for (i1 = 0; i1 < na; i1++)
   {
      v1 = periodic_table[as[i1]].valence;
      r1 = periodic_table[as[i1]].rcov;
      ii = 0;
      for (i2 = 0; i2 < una; i2++)
      {
         r2 = periodic_table[uas[i2]].rcov;
         for (j1 = scmin[0] - 1; j1 < scmax[0] + 2; j1++)
         for (j2 = scmin[1] - 1; j2 < scmax[1] + 2; j2++)
         for (j3 = scmin[2] - 1; j3 < scmax[2] + 2; j3++)
         if ((j3 != sci[i1][3] || j2 != sci[i1][2] || j1 != sci[i1][1] || i2 != sci[i1][0]) && uas[i2] != 0 && as[i1] != 0)
         {
            pp.xx = ucv[2].xx * double(j3) + ucv[1].xx * double(j2) + ucv[0].xx * double(j1) + uap[i2].xx - ap[i1].xx;
            pp.yy = ucv[2].yy * double(j3) + ucv[1].yy * double(j2) + ucv[0].yy * double(j1) + uap[i2].yy - ap[i1].yy;
            pp.zz = ucv[2].zz * double(j3) + ucv[1].zz * double(j2) + ucv[0].zz * double(j1) + uap[i2].zz - ap[i1].zz;
            rr = sqrt(pp.xx * pp.xx + pp.yy * pp.yy + pp.zz * pp.zz);
            if (rr < MINIMUMBONDINGDISTANCE)
               return -1;
            if (rr < r1 + r2 + BONDTOLERANCE)
            {
               if (ii >= MAXBOND)
                  return -1;
               if (STRICTVALENCE && ii >= v1)
                  return -1;
               bonds[i1][ii].xx = ap[i1].xx + pp.xx * r1 / (r1 + r2);
               bonds[i1][ii].yy = ap[i1].yy + pp.yy * r1 / (r1 + r2);
               bonds[i1][ii].zz = ap[i1].zz + pp.zz * r1 / (r1 + r2);
               ii++;
            }
         }
      }
      nb[i1] = ii;
   }

   return 0;
}

CARTESIAN interpolate_vertex(double scalar, CARTESIAN cartesian1, CARTESIAN cartesian2, double scalar1, double scalar2)
{
   double weight1, weight2;
   CARTESIAN cartesian;

   if (fabs(scalar1 - scalar) < fabs(scalar2 - scalar1) * EPS9)
      return cartesian1;
   if (fabs(scalar2 - scalar) < fabs(scalar2 - scalar1) * EPS9)
      return cartesian2;

   if (fabs(scalar2 - scalar1) < EPS9)
   {
      weight1 = 0.5;
      weight2 = 0.5;
   }
   else
   {
      weight1 = (scalar2 - scalar) / (scalar2 - scalar1);
      weight2 = (scalar - scalar1) / (scalar2 - scalar1);
   }

   cartesian.xx = weight1 * cartesian1.xx + weight2 * cartesian2.xx;
   cartesian.yy = weight1 * cartesian1.yy + weight2 * cartesian2.yy;
   cartesian.zz = weight1 * cartesian1.zz + weight2 * cartesian2.zz;

   return cartesian;
}

int polygonise_cube(CUBE cube, double isovalue, TRIANGLE *triangle, NORMAL *normal, int smooth)
{
   int ii, jj, kk = 0, ll = 1, mm = 1, tt = 0;
   CARTESIAN point[12], vector[12];

   for (ii = 0; ii < 8; ii++)
   {
      if (cube.scalar[ii] < isovalue)
         kk |= ll;
      ll <<= 1;
   }

   if (edge_table_cube[kk] == 0)
      return 0;

   for (ii = 0; ii < 12; ii++)
   {
      if (edge_table_cube[kk] & mm)
      {
         point[ii] = interpolate_vertex(isovalue, cube.point[edge[ii][0]], cube.point[edge[ii][1]], cube.scalar[edge[ii][0]], cube.scalar[edge[ii][1]]);
         if (smooth)
            vector[ii] = interpolate_vertex(isovalue, cube.vector[edge[ii][0]], cube.vector[edge[ii][1]], cube.scalar[edge[ii][0]], cube.scalar[edge[ii][1]]);
      }
      mm <<= 1;
   }

   for (ii = 0; triangle_table_cube[kk][ii] != -1; ii += 3)
   {
      for (jj = 0; jj < 3; jj++)
         triangle[tt].point[jj] = point[triangle_table_cube[kk][ii + jj]];
      if (smooth)
         for (jj = 0; jj < 3; jj++)
            normal[tt].vector[jj] = vector[triangle_table_cube[kk][ii + jj]];
      tt++;
   }

   return tt;
}

int polygonise_tetrahedron(TETRAHEDRON tetrahedron, double isovalue, TRIANGLE *triangle, NORMAL *normal, int smooth)
{
   bool flag[6];
   int ii, jj, kk = 0, ll = 1, mm = 1, tt = 0;
   CARTESIAN point[6], vector[6];

   for (ii = 0; ii < 4; ii++)
   {
      if (tetrahedron.scalar[ii] < isovalue)
         kk |= ll;
      ll <<= 1;
   }

   if (edge_table_tetrahedron[kk] == 0)
      return 0;

   for (ii = 0; ii < 6; ii++)
   {
      if (edge_table_tetrahedron[kk] & mm)
         flag[ii] = true;
      else
         flag[ii] = false;
      mm <<= 1;
   }

   for (ii = 0; ii < 6; ii++)
   {
   if (flag[ii])
      point[ii] = interpolate_vertex(isovalue, tetrahedron.point[edge_tetrahedron[ii][0]], tetrahedron.point[edge_tetrahedron[ii][1]], tetrahedron.scalar[edge_tetrahedron[ii][0]], tetrahedron.scalar[edge_tetrahedron[ii][1]]);
   if (smooth && flag[ii])
      vector[ii] = interpolate_vertex(isovalue, tetrahedron.vector[edge_tetrahedron[ii][0]], tetrahedron.vector[edge_tetrahedron[ii][1]], tetrahedron.scalar[edge_tetrahedron[ii][0]], tetrahedron.scalar[edge_tetrahedron[ii][1]]);
   }

   for (ii = 0; triangle_table_tetrahedron[kk][ii] != -1; ii += 3)
   {
      for (jj = 0; jj < 3; jj++)
         triangle[tt].point[jj] = point[triangle_table_tetrahedron[kk][ii + jj]];
      if (smooth)
      for (jj = 0; jj < 3; jj++)
         normal[tt].vector[jj] = vector[triangle_table_tetrahedron[kk][ii + jj]];
      tt++;
   }

   return tt;
}

void polygonise(double isovalue, int algorithm, int smooth, int sign, const scalar_field & sf)
{
   int ii, jj, kk, cc, aa, tt, na, nt, lo, hi;
   CUBE cube;
   TETRAHEDRON tetrahedron;
   TRIANGLE polygon_triangle[5];
   NORMAL polygon_normal[5];
   TRIANGLELIST *previous_triangle, *current_triangle, *first_triangle = NULL;
   NORMALLIST *previous_normal, *current_normal, *first_normal = NULL;

   if (algorithm == 0)
      na = 1;
   else
      na = 6;

   for (ii = 0; ii < sf.ni - 1; ii++)
   for (jj = 0; jj < sf.nj - 1; jj++)
   for (kk = 0; kk < sf.nk - 1; kk++)
   {
      lo = 0;
      hi = 0;
      for (cc = 0; cc < 8; cc++)
      {
         if (sf.scalar[ii + vertex[cc][0]][jj + vertex[cc][1]][kk + vertex[cc][2]] <= isovalue)
            lo++;
         if (sf.scalar[ii + vertex[cc][0]][jj + vertex[cc][1]][kk + vertex[cc][2]] >= isovalue)
            hi++;
      }
      if (lo != 0 && hi != 0)
      {
         for (cc = 0; cc < 8; cc++)
         {
            cube.point[cc].xx = sf.origin.xx + double(ii + vertex[cc][0]) * sf.stepv[0].xx + double(jj + vertex[cc][1]) * sf.stepv[1].xx + double(kk + vertex[cc][2]) * sf.stepv[2].xx;
            cube.point[cc].yy = sf.origin.yy + double(ii + vertex[cc][0]) * sf.stepv[0].yy + double(jj + vertex[cc][1]) * sf.stepv[1].yy + double(kk + vertex[cc][2]) * sf.stepv[2].yy;
            cube.point[cc].zz = sf.origin.zz + double(ii + vertex[cc][0]) * sf.stepv[0].zz + double(jj + vertex[cc][1]) * sf.stepv[1].zz + double(kk + vertex[cc][2]) * sf.stepv[2].zz;
            if (smooth)
               cube.vector[cc] = vector[ii + vertex[cc][0]][jj + vertex[cc][1]][kk + vertex[cc][2]];
            cube.scalar[cc] = sf.scalar[ii + vertex[cc][0]][jj + vertex[cc][1]][kk + vertex[cc][2]];
         }

         for (aa = 0; aa < na; aa++)
         {
            if (algorithm == 0)
               nt = polygonise_cube(cube, isovalue, polygon_triangle, polygon_normal, smooth);
            else
            {
               for (tt = 0; tt < 4; tt++)
               {
                  tetrahedron.point[tt] = cube.point[vertex_tetrahedron[aa][tt]];
                  if (smooth)
                     tetrahedron.vector[tt] = cube.vector[vertex_tetrahedron[aa][tt]];
                  tetrahedron.scalar[tt] = cube.scalar[vertex_tetrahedron[aa][tt]];
               }
               nt = polygonise_tetrahedron(tetrahedron, isovalue, polygon_triangle, polygon_normal, smooth);
            }
            for (tt = 0; tt < nt; tt++)
            {
               current_triangle = new TRIANGLELIST;
               current_triangle->triangle = polygon_triangle[tt];
               current_triangle->node = NULL;
               if (first_triangle == NULL)
                  first_triangle = current_triangle;
               else
                  previous_triangle->node = current_triangle;
               previous_triangle = current_triangle;
               if (smooth)
               {
                  current_normal = new NORMALLIST;
                  current_normal->normal = polygon_normal[tt];
                  current_normal->node = NULL;
                  if (first_normal == NULL)
                     first_normal = current_normal;
                  else
                     previous_normal->node = current_normal;
                  previous_normal = current_normal;
               }
            }
         }
      }
   }

   triangle[sign] = first_triangle;
   if (smooth)
      normal[sign] = first_normal;
}

int pov_write(char *ofn, char *header, int sign, bool smooth, bool box, bool basis, int na, int ns, CARTESIAN sco, CARTESIAN *scv)
{
   const CARTESIAN cameralocation = {-0.1, -0.5, 0.9};
   const CARTESIAN cameralookat = {0.0, 0.0, 0.0};
   const double camerafactor = 2.0;

   bool flag;
   int ii, jj, ss, ierr;
   double scscale, cameradistance, camerascale, axislength, axisradius, ratom, rbond, color[3];
   double xmin = INF9, ymin = INF9, zmin = INF9, xmax = -INF9, ymax = -INF9, zmax = -INF9, dx, dy, dz;
   char symbol[4];
   CARTESIAN sccenter, sccorner[8], cameralocationsc, cameralookatsc;
   TRIANGLELIST *current_triangle;
   NORMALLIST *current_normal;
   FILE *hh;

   hh = fopen(ofn, "w");
   if (hh == NULL)
      return -1;

   sccenter.xx = sco.xx + scv[0].xx * 0.5 + scv[1].xx * 0.5 + scv[2].xx * 0.5;
   sccenter.yy = sco.yy + scv[0].yy * 0.5 + scv[1].yy * 0.5 + scv[2].yy * 0.5;
   sccenter.zz = sco.zz + scv[0].zz * 0.5 + scv[1].zz * 0.5 + scv[2].zz * 0.5;
   for (ii = 0; ii < 8; ii++)
   {
      sccorner[ii].xx = sccenter.xx;
      sccorner[ii].yy = sccenter.yy;
      sccorner[ii].zz = sccenter.zz;
      for (jj = 0; jj < 3; jj++)
      {
         sccorner[ii].xx += scv[jj].xx * (double(vertex[ii][jj]) - 0.5);
         sccorner[ii].yy += scv[jj].yy * (double(vertex[ii][jj]) - 0.5);
         sccorner[ii].zz += scv[jj].zz * (double(vertex[ii][jj]) - 0.5);
      }
   }

   scscale = 0.0;
   for (ii = 0; ii < 3; ii++)
      scscale += sqrt(scv[ii].xx * scv[ii].xx + scv[ii].yy * scv[ii].yy + scv[ii].zz * scv[ii].zz);
   scscale = scscale / 3.0;
   cameradistance = sqrt((cameralocation.xx - cameralookat.xx) * (cameralocation.xx - cameralookat.xx) + (cameralocation.yy - cameralookat.yy) * (cameralocation.yy - cameralookat.yy) + (cameralocation.zz - cameralookat.zz) * (cameralocation.zz - cameralookat.zz));
   camerascale = camerafactor * scscale / cameradistance;
   cameralocationsc.xx = sccenter.xx / camerascale + cameralocation.xx;
   cameralocationsc.yy = sccenter.yy / camerascale + cameralocation.yy;
   cameralocationsc.zz = sccenter.zz / camerascale + cameralocation.zz;
   cameralookatsc.xx = sccenter.xx / camerascale + cameralookat.xx;
   cameralookatsc.yy = sccenter.yy / camerascale + cameralookat.yy;
   cameralookatsc.zz = sccenter.zz / camerascale + cameralookat.zz;

   fprintf(hh, "\n");
   fprintf(hh, "// %s\n",header);
   fprintf(hh, "\n");
   fprintf(hh, "#declare camera_location = <%.2f, %.2f, %.2f>;\n", cameralocationsc.xx, cameralocationsc.yy, cameralocationsc.zz);
   fprintf(hh, "#declare camera_look_at = <%.2f, %.2f, %.2f>;\n", cameralookatsc.xx, cameralookatsc.yy, cameralookatsc.zz);
   fprintf(hh, "#declare camera_scale = %.2f;\n", camerascale);
   fprintf(hh, "#declare light_location = camera_location - camera_look_at;\n");
   fprintf(hh, "#declare light_scale = 1e6;\n");
   fprintf(hh, "#declare color_light = rgb <2.00, 2.00, 2.00>;\n");
   fprintf(hh, "#declare color_background = rgb <0.00, 0.00, 0.00>;\n");
   if (box)
   {
      fprintf(hh, "#declare radius_frame = 0.01;\n");
      fprintf(hh, "#declare color_frame = rgb <0.75, 0.75, 0.75>;\n");
      fprintf(hh, "#declare color_box = rgbf <1.00, 1.00, 1.00, 0.75>;\n");
   }
   if (basis)
   {
      axislength = scscale / 3.0;
      axisradius = axislength / 50.0;
      fprintf(hh, "#declare length_axis = %.2f;\n", axislength);
      fprintf(hh, "#declare radius_axis = %.2f;\n", axisradius);
      fprintf(hh, "#declare color_axis_x = rgb <1.00, 0.00, 0.00>;\n");
      fprintf(hh, "#declare color_axis_y = rgb <0.00, 1.00, 0.00>;\n");
      fprintf(hh, "#declare color_axis_z = rgb <0.00, 0.00, 1.00>;\n");
      fprintf(hh, "#declare length_arrow = 0.2 * length_axis;\n");
      fprintf(hh, "#declare radius_arrow = 2.0 * radius_axis;\n");
      fprintf(hh, "#declare color_arrow_x = color_axis_x;\n");
      fprintf(hh, "#declare color_arrow_y = color_axis_y;\n");
      fprintf(hh, "#declare color_arrow_z = color_axis_z;\n");
   }
   for (ii = 0; ii < ns; ii++)
   {
      strncpy(symbol, periodic_table[cs[ii]].symbol, 4);
      symbol[3] = '\0';
      lowercase(symbol);
      ratom = periodic_table[cs[ii]].rvdw * ATOMSIZE;
      rbond = BONDRADIUS;
      for (jj = 0; jj < 3; jj++)
         color[jj] = periodic_table[cs[ii]].color[jj];
      fprintf(hh, "#declare radius_atom_%s = %.2f;\n", symbol, ratom);
      if (cs[ii] != 0)
         fprintf(hh, "#declare radius_bond_%s = %.2f;\n", symbol, rbond);
      fprintf(hh, "#declare color_atom_%s = rgb <%.2f, %.2f, %.2f>;\n", symbol, color[0], color[1], color[2]);
      if (cs[ii] != 0)
         fprintf(hh, "#declare color_bond_%s = rgb <%.2f, %.2f, %.2f>;\n", symbol, color[0], color[1], color[2]);
   }
   if (sign == 2)
   {
      fprintf(hh, "#declare color_isosurface_positive = rgbf <0.40, 0.40, 0.80, 0.75>;\n");
      fprintf(hh, "#declare color_isosurface_negative = rgbf <0.80, 0.40, 0.40, 0.75>;\n");
   }
   else
      fprintf(hh, "#declare color_isosurface = rgbf <0.40, 0.80, 0.40, 0.75>;\n");
   fprintf(hh, "\n");
   fprintf(hh, "camera { location camera_location sky <0.00, 0.00, 1.00> up <0.00, 0.00, 1.00> right <-1.33, 0.00, 0.00> direction <0.00, -1.00, 0.00> look_at camera_look_at scale camera_scale }\n");
   fprintf(hh, "light_source { light_location color color_light shadowless scale light_scale }\n");
   fprintf(hh, "background { color color_background }\n");
   fprintf(hh, "\n");
   if (box)
   {
      fprintf(hh, "union {\n");
      for (ii = 0; ii < 8; ii++)
         fprintf(hh, "sphere { <%.9f, %.9f, %.9f>, radius_frame }\n", sccorner[ii].xx, sccorner[ii].yy, sccorner[ii].zz);
      for (ii = 0; ii < 12; ii++)
         fprintf(hh, "cylinder { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, radius_frame }\n", sccorner[edge[ii][0]].xx, sccorner[edge[ii][0]].yy, sccorner[edge[ii][0]].zz, sccorner[edge[ii][1]].xx, sccorner[edge[ii][1]].yy, sccorner[edge[ii][1]].zz);
      fprintf(hh, "texture { pigment { color color_frame } }\n");
      fprintf(hh, "}\n");
      fprintf(hh, "\n");
      fprintf(hh, "union {\n");
      for (ii = 0; ii < 12; ii++)
         fprintf(hh, "triangle { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f> }\n", sccorner[face[ii][0]].xx, sccorner[face[ii][0]].yy, sccorner[face[ii][0]].zz, sccorner[face[ii][1]].xx, sccorner[face[ii][1]].yy, sccorner[face[ii][1]].zz, sccorner[face[ii][2]].xx, sccorner[face[ii][2]].yy, sccorner[face[ii][2]].zz);
      fprintf(hh, "texture { pigment { color color_box } }\n");
      fprintf(hh, "}\n");
      fprintf(hh, "\n");
   }
   if (basis)
   {
      fprintf(hh, "union {\n");
      fprintf(hh, "cylinder { <0.00, 0.00, 0.00>, <length_axis, 0.00, 0.00>, radius_axis texture { pigment { color color_axis_x } } }\n");
      fprintf(hh, "cylinder { <0.00, 0.00, 0.00>, <0.00, length_axis, 0.00>, radius_axis texture { pigment { color color_axis_y } } }\n");
      fprintf(hh, "cylinder { <0.00, 0.00, 0.00>, <0.00, 0.00, length_axis>, radius_axis texture { pigment { color color_axis_z } } }\n");
      fprintf(hh, "cone { <length_axis, 0.00, 0.00>, radius_arrow <length_axis + length_arrow, 0.00, 0.00>, 0.00 texture { pigment { color color_arrow_x } } }\n");
      fprintf(hh, "cone { <0.00, length_axis, 0.00>, radius_arrow <0.00, length_axis + length_arrow, 0.00>, 0.00 texture { pigment { color color_arrow_y } } }\n");
      fprintf(hh, "cone { <0.00, 0.00, length_axis>, radius_arrow <0.00, 0.00, length_axis + length_arrow>, 0.00 texture { pigment { color color_arrow_z } } }\n");
      fprintf(hh, "}\n");
      fprintf(hh, "\n");
   }
   fprintf(hh, "union {\n");
   for (ii = 0; ii < na; ii++)
   {
      strncpy(symbol, periodic_table[as[ii]].symbol, 4);
      symbol[3] = '\0';
      lowercase(symbol);
      fprintf(hh, "sphere { <%.9f, %.9f, %.9f>, radius_atom_%s texture { pigment { color color_atom_%s } } }\n", ap[ii].xx, ap[ii].yy, ap[ii].zz, symbol, symbol);
   }
   fprintf(hh, "}\n");
   fprintf(hh, "\n");
   fprintf(hh, "union {\n");
   for (ii = 0; ii < na; ii++)
   {
      strncpy(symbol, periodic_table[as[ii]].symbol, 4);
      symbol[3] = '\0';
      lowercase(symbol);
      for (jj = 0; jj < nb[ii]; jj++)
         fprintf(hh, "cylinder { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, radius_bond_%s texture { pigment { color color_bond_%s } } }\n", ap[ii].xx, ap[ii].yy, ap[ii].zz, bonds[ii][jj].xx, bonds[ii][jj].yy, bonds[ii][jj].zz, symbol, symbol);
   }
   fprintf(hh, "}\n");
   fprintf(hh, "\n");
   for (ss = 0; ss < 2; ss++)
   {
      current_triangle = triangle[ss];
      if (smooth)
         current_normal = normal[ss];
      if (current_triangle != NULL)
      {
         fprintf(hh, "mesh {\n");
         do
         {
            for (ii = 0; ii < 3; ii++)
            {
               if (current_triangle->triangle.point[ii].xx < xmin) xmin = current_triangle->triangle.point[ii].xx;
               if (current_triangle->triangle.point[ii].xx > xmax) xmax = current_triangle->triangle.point[ii].xx;
               if (current_triangle->triangle.point[ii].yy < ymin) ymin = current_triangle->triangle.point[ii].yy;
               if (current_triangle->triangle.point[ii].yy > ymax) ymax = current_triangle->triangle.point[ii].yy;
               if (current_triangle->triangle.point[ii].zz < zmin) zmin = current_triangle->triangle.point[ii].zz;
               if (current_triangle->triangle.point[ii].zz > zmax) zmax = current_triangle->triangle.point[ii].zz;
            }
            if (smooth)
               flag = fabs(current_normal->normal.vector[0].xx) + fabs(current_normal->normal.vector[0].yy) + fabs(current_normal->normal.vector[0].zz) < EPS9 || fabs(current_normal->normal.vector[1].xx) + fabs(current_normal->normal.vector[1].yy) + fabs(current_normal->normal.vector[1].zz) < EPS9 || fabs(current_normal->normal.vector[2].xx) + fabs(current_normal->normal.vector[2].yy) + fabs(current_normal->normal.vector[2].zz) < EPS9;
            else
               flag = true;
            if (flag)
               fprintf(hh, "triangle { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f> }\n", current_triangle->triangle.point[0].xx, current_triangle->triangle.point[0].yy, current_triangle->triangle.point[0].zz, current_triangle->triangle.point[1].xx, current_triangle->triangle.point[1].yy, current_triangle->triangle.point[1].zz, current_triangle->triangle.point[2].xx, current_triangle->triangle.point[2].yy, current_triangle->triangle.point[2].zz);
            else
               fprintf(hh, "smooth_triangle { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f> }\n", current_triangle->triangle.point[0].xx, current_triangle->triangle.point[0].yy, current_triangle->triangle.point[0].zz, current_normal->normal.vector[0].xx, current_normal->normal.vector[0].yy, current_normal->normal.vector[0].zz, current_triangle->triangle.point[1].xx, current_triangle->triangle.point[1].yy, current_triangle->triangle.point[1].zz, current_normal->normal.vector[1].xx, current_normal->normal.vector[1].yy, current_normal->normal.vector[1].zz, current_triangle->triangle.point[2].xx, current_triangle->triangle.point[2].yy, current_triangle->triangle.point[2].zz, current_normal->normal.vector[2].xx, current_normal->normal.vector[2].yy, current_normal->normal.vector[2].zz);
            current_triangle = current_triangle->node;
            if (smooth)
               current_normal = current_normal->node;
         }
         while (current_triangle != NULL);
         if (sign == 2)
         {
            if (ss == 0)
               fprintf(hh, "texture { pigment { color color_isosurface_positive } }\n");
            else
               fprintf(hh, "texture { pigment { color color_isosurface_negative } }\n");
         }
         else
            fprintf(hh, "texture { pigment { color color_isosurface } }\n");
         fprintf(hh, "}\n");
         fprintf(hh, "\n");
      }
   }

   ierr = fclose(hh);
   if (ierr != 0)
      return -1;

   xmin = xmin / BOHR;
   xmax = xmax / BOHR;
   ymin = ymin / BOHR;
   ymax = ymax / BOHR;
   zmin = zmin / BOHR;
   zmax = zmax / BOHR;
   dx = xmax - xmin;
   dy = ymax - ymin;
   dz = zmax - zmin;
   printf("\n the minimal bounding box for the isosurface\n\n xmin = %.9f xmax = %.9f dx = %.9f bohr\n ymin = %.9f ymax = %.9f dy = %.9f bohr\n zmin = %.9f zmax = %.9f dz = %.9f bohr\n\n", xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz);

   return 0;
}

int terminate(int ierr, int na, int ni, int nj, scalar_field & sf)
{
   int ii, jj, ss;
   TRIANGLELIST *current_triangle, *next_triangle;
   NORMALLIST *current_normal, *next_normal;

   if (as != NULL) delete [] as;
   if (ap != NULL) delete [] ap;
   if (uas != NULL) delete [] uas;
   if (uap != NULL) delete [] uap;
   if (cs != NULL) delete [] cs;
   if (nb != NULL) delete [] nb;

   if (bonds != NULL)
   {
      for (ii = 0; ii < na; ii++)
         delete [] bonds[ii];
      delete [] bonds;
   }

   if (sci != NULL)
   {
      for (ii = 0; ii < na; ii++)
         delete [] sci[ii];
      delete [] sci;
   }

   if (sf.scalar != NULL)
   {
      deallocate_scalar(sf.scalar, ni, nj);
   }

   if (vector != NULL)
   {
      for (ii = 0; ii < ni; ii++)
      for (jj = 0; jj < nj; jj++)
         delete [] vector[ii][jj];
      for (ii = 0; ii < ni; ii++)
         delete [] vector[ii];
      delete [] vector;
   }

   for (ss = 0; ss < 2; ss++)
   {
      current_triangle = triangle[ss];
      while (current_triangle != NULL)
      {
         next_triangle = current_triangle->node;
         delete current_triangle;
         current_triangle = next_triangle;
      }
   }

   for (ss = 0; ss < 2; ss++)
   {
      current_normal = normal[ss];
      while (current_normal != NULL)
      {
         next_normal = current_normal->node;
         delete current_normal;
         current_normal = next_normal;
      }
   }

   return ierr;
}

int main(int argc, char *argv[])
{
   bool smooth, box, basis, uc, sc, sct;
   int sign, power, algorithm, ucf, scf, ns, una;
   int scmin[3], scmax[3];
   int ierr = 0, na = 0, ni = 0, nj = 0;
   double isovalue;
   char pfn[MAXCHAR], ifn[MAXCHAR], iff[MAXCHAR], ofn[MAXCHAR], off[MAXCHAR], header[MAXCHAR];
   CARTESIAN sfv[3], uco, ucv[3], sco, scv[3];
   scalar_field sf;

   sf.scalar = NULL;

   if (argc != 2)
   {
      fprintf(stderr, "\nUsage: surface.x surface.inp\n\n");
      return terminate(-1, na, ni, nj, sf);
   }
   strncpy(pfn, argv[1], MAXCHAR);
   pfn[MAXCHAR - 1] = '\0';

   ierr = par_read(pfn, header, ifn, iff, ofn, off, &isovalue, &sign, &power, &algorithm, &smooth, &box, &basis, &uc, &uco, ucv, &ucf, &sc, &sco, scv, &scf, &sct);
   if (ierr != 0)
   {
      fprintf(stderr, "\nError: failed to read %s\n\n", pfn);
      return terminate(-1, na, ni, nj, sf);
   }

   if (strcmp(iff, "cube") == 0)
      ierr = cub_read(ifn, &na, sfv, sf);
   else if (strcmp(iff, "xsf") == 0)
      ierr = xsf_read(ifn, &na, sfv, sf);
   else
   {
      fprintf(stderr, "\nError: unrecognized input format %s\n\n", iff);
      return terminate(-1, na, sf.ni, sf.nj, sf);
   }
   if (ierr != 0)
   {
      fprintf(stderr, "\nError: failed to read %s\n\n", ifn);
      return terminate(-1, na, sf.ni, sf.nj, sf);
   }

   cell_set(sf.origin, sfv, uc, &uco, ucv, ucf, sc, &sco, scv, scf);

   if (sc)
   {
      una = na;
      na = atom_clone(na, ucv, sco, scv, scmin, scmax);
      ierr = scalar_clone(uco, ucv, sco, scv, sf);
   }
   if (na < 0 || ierr != 0)
   {
      fprintf(stderr, "\nError: failed to build supercell\n\n");
      return terminate(-1, na, sf.ni, sf.nj, sf);
   }

   if (smooth)
      scalar_normal(sf.ni, sf.nj, sf.nk, sf);

   ns = chemical_species(na);

   if (sct || !sc)
      ierr = bonds_st(na, scv);
   else
      ierr = bonds_ut(na, una, ucv, scmin, scmax);
   if (ierr != 0)
   {
      fprintf(stderr, "\nError: failed to build interatomic bonds\n\n");
      return terminate(-1, na, sf.ni, sf.nj, sf);
   }

   isovalue = isovalue_scale(power, isovalue, sf);
   if (sign == 0 || sign == 2)
      polygonise(isovalue, algorithm, smooth, 0, sf);
   isovalue = -isovalue;
   if (sign == 1 || sign == 2)
      polygonise(isovalue, algorithm, smooth, 1, sf);

   if (strcmp(off, "povray") == 0)
      ierr = pov_write(ofn, header, sign, smooth, box, basis, na, ns, sco, scv);
   else
   {
      fprintf(stderr, "\nError: unrecognized output format %s\n\n", off);
      return terminate(-1, na, ni, nj, sf);
   }
   if (ierr != 0)
   {
      fprintf(stderr, "\nError: failed to write %s\n\n", ofn);
      return terminate(-1, na, ni, nj, sf);
   }

   return terminate(0, na, ni, nj, sf);
}

