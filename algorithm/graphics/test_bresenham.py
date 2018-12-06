#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np


size_x = 50

def bresenham_2d_general(at, p0, p1, max_length=np.inf):
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = dp = p1 - p0
    abs_dx, abs_dy = abs_dp = np.abs(dp)

    offset_dx = np.sign(dx)
    offset_dy = np.sign(dy) * size_x

    offset = y0 * size_x + x0

    dist = np.linalg.norm(dp)
    scale = 1.0 if (dist == 0.0) else np.min([1.0, max_length / dist]);

    if(abs_dx >= abs_dy):
        error_y = abs_dx / 2
        bresenham_2d(at, abs_dx, abs_dy, error_y, offset_dx, offset_dy, offset, int(scale*abs_dx))
    else:
        error_x = abs_dy / 2
        bresenham_2d(at, abs_dy, abs_dx, error_x, offset_dy, offset_dx, offset, int(scale*abs_dy))
        

def bresenham_2d(at, abs_da, abs_db, error_b, offset_a, offset_b, offset, max_length):
    end = np.min([max_length, abs_da])

    for i in xrange(0, end):
        at(offset)
        offset += offset_a
        error_b += abs_db
        if(int(error_b) >= abs_da):
            offset += offset_b
            error_b -= abs_da

    at(offset)


point_array = []

def draw_point(offset):
    x = offset % size_x
    y = offset / size_x
    print x,y
    global point_array
    point_array[x][y] = 1
    #plt.scatter(x,y, s=s,color='red', marker='s')


def test(p):
    p = np.array(p)
    print p
    global point_array
    point_array = np.zeros(p[1]+[1,1])
    print point_array,point_array.shape
    bresenham_2d_general(draw_point, p[0], p[1])
    a = []
    for i in xrange(point_array.shape[0]):
        for j in xrange(point_array.shape[1]):
            a.append([i,j, point_array[i][j]])

    for i in xrange(point_array.shape[0]+1):
        plt.plot([i-0.5,i-0.5],[p[0,1]-0.5,p[1,1]+0.5],color='black')
    for j in xrange(point_array.shape[1]+1):
        plt.plot([p[0,0]-0.5,p[1,0]+0.5],[j-0.5,j-0.5],color='black')

    print "="
    print np.array(a)[:,:2].T.shape
    plt.scatter(*np.array(a)[:,:2].T, s=100, color=map(lambda x:'red' if x==1 else 'orange',np.array(a)[:,2]), marker='s')
    plt.grid(True)
    plt.axes().set_aspect('equal','datalim')
    plt.plot(*p.T, color='blue')
    plt.show()

test([[0,0], 
      [40,15]])
