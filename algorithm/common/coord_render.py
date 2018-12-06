#!/usr/bin/env python

import pylab
import numpy as np
from matplotlib import pyplot as plt


def plot_show():
    plt.grid(True)
    plt.axes().set_aspect('equal','datalim')
    plt.show()

def plot_savefig(*args, **kwargs):
    plt.grid(True)
    plt.axes().set_aspect('equal','datalim')
    plt.savefig(*args, **kwargs)


def plot_arrow_xy_d(x,y,
                    dx,dy,
                    color='black',
                    width=0.01,
                    head_width=0.05):
    #print x,y," ",dx,dy," ",x+dx,y+dy," ",color
    plt.plot([x,x+dx],[y,y+dy],color=color)
    plt.arrow(x,y,
              dx,dy,
              width=width,
              shape="full",
              lw=0,
              length_includes_head=True,
              head_width=head_width,
              color=color)


def plot_arrow_xy(x0,y0,
                  x1,y1,
                  color="black",
                  width=0.01,
                  head_width=0.05):
    plot_arrow_xy_d(x0,y0,
                    x1-x0, y1-y0,
                    color,
                    width,
                    head_width)

def plot_arrow_v(v0, v1, 
                 color='black',
                 width=0.01,
                 head_width=0.05):
    plot_arrow_xy(v0[0],v0[1],
                  v1[0],v1[1],
                  color,
                  width,
                  head_width)

def plot_arrow_v_d(v, dv, 
                   color='black',
                   width=0.01,
                   head_width=0.05):
    plot_arrow_xy_d(v[0],v[1],
                    dv[0],dv[1],
                    color,
                    width,
                    head_width)
    
def get_2d_transform_mat(x,y,theta):
    #print '[ get_2d_transform_mat ]  ',x,y,theta
    t = theta
    return np.array([[  np.cos(t), -np.sin(t), x ],
                     [  np.sin(t),  np.cos(t), y ],
                     [  0        ,  0        , 1 ]])

def get_2d_rotation_mat(theta):
    t = theta
    return np.array([[  np.cos(t), -np.sin(t) ],
                     [  np.sin(t),  np.cos(t) ]])

def rotatation_mat_2d_to_theta(T):
    return np.arctan2(T[1,0],T[0,0])
    
def transform_mat_2d_to_xyt(T):
    return np.array([T[0,2], T[1,2], np.arctan2(T[1,0],T[0,0])])

def plot_2d_transform_xyt( x, y, theta, length=1.0):
    length = float(length)
    T_rot = get_2d_rotation_mat(theta)

    x_arrow = T_rot[:][0] * length
    y_arrow = T_rot[:][1] * length
    
    O = np.array([x,y])

    plot_arrow_v_d( O, x_arrow, color='red')
    plot_arrow_v_d( O, y_arrow, color='green')


def plot_2d_transform_T( T, length=1.0):
    #print T
    length = float(length)
    x_arrow = T[:2,0] * length
    y_arrow = T[:2,1] * length
    
    O = T[:2,2]

    plot_arrow_v_d( O, x_arrow, color='red')
    plot_arrow_v_d( O, y_arrow, color='green')


def mod_angle(a):
    return (a+np.pi)%(2*np.pi)-np.pi


def dot_chain(*args):
    if(len(args)==1):
        return args[0]
    T = args[0]
    for i in args[1:]:
        T = np.dot(T,i)
    return T
    
# TODO: use currying
def rotate_around_point(p,theta,T):
    theta = mod_angle(theta)
    T_o_c = get_2d_transform_mat(p[0], p[1], 0) # rotate center, relative to the Origin
    T_r = get_2d_transform_mat(0, 0, theta)   # rotation matrix
    T1 = dot_chain(
        T_o_c,
        T_r,
        np.linalg.inv(T_o_c),
        T
    )
    return T1
     
def move_arc( radius, theta, T ):
    T_c = get_2d_transform_mat(0, radius, 0) 
    T_r = get_2d_transform_mat(0, 0, theta)   # rotation matrix
    T1 = dot_chain(
        T, T_c,
        T_r,
        np.linalg.inv(T_c)
        #,  np.linalg.inv(T),  # they become I
        #T
    )
    return T1

def move_arc_param( radius, theta, T ):
    #print "[ move_arc_param ] [ radius, theta, T]",radius,theta
    #print T
    theta_T = np.arctan2(T[1][0], T[0][0])
    T_c = get_2d_transform_mat(0, radius, 0) 
    T_r = get_2d_transform_mat(0, 0, theta)   # rotation matrix
    T_o_c = dot_chain(
        T, T_c
    )
    center_x, center_y = T_o_c[:2,2]
    T_c_k0 = np.linalg.inv(T_c)
    dx, dy = T_c_k0[:2,2]
    angle_start = np.arctan2(dy, dx) + theta_T
    #angle_start = np.arctan2(dx, dy)
    #print 'T_c_k0 dx dy: ',dx,dy,np.arctan2(dy,dx),theta_T,angle_start
    T_c_k1 = dot_chain(
        T_r, T_c_k0
    )
    dx, dy = T_c_k1[:2,2]
    angle_stop = np.arctan2(dy, dx) + theta_T
    #angle_stop = np.arctan2(dx, dy)
    #print 'T_c_k1 dx dy: ',dx,dy,np.arctan2(dy,dx),theta_T, angle_stop
    if(1):
        plot_t = np.linspace(angle_start, angle_stop, 30)
        plot_x = center_x + np.abs(radius) * np.cos(plot_t)
        plot_y = center_y + np.abs(radius) * np.sin(plot_t)
        plt.plot(plot_x, plot_y,color='blue')
    #print ("["
    #       " center_x,"
    #       " center_y,"
    #       " angle_start,"
    #       " angle_stop ]"+" % .3f"*4)%(center_x, 
    #                                    center_y, 
    #                                    angle_start, 
    #                                    angle_stop)
    return center_x, center_y, angle_start, angle_stop
    
def move_straight( length, T ):
    T_t = get_2d_transform_mat( length, 0, 0 )
    T1 = dot_chain(
        T,
        T_t
    )
    return T1

if __name__ == "__main__":
    x = np.linspace(0,np.pi*2,100)
    plt.plot(np.cos(x),np.sin(x))
    T_center = get_2d_transform_mat( 0, 0, np.pi/6 )
    T_mover = get_2d_transform_mat( 1, 0, np.pi/2 ) 
    #plot_2d_transform_T( T_center, length=0.3 )
    plot_2d_transform_T( T_mover, length=0.3 )
    for i in xrange(30):
        T_mover = move_arc( 1, np.pi/15, T_mover)
        plot_2d_transform_T( T_mover, length=0.3 )

    if(False):
        T_mover = rotate_around_point( np.array([0, 0]), np.pi/6, T_mover)
        plot_2d_transform_T( T_mover, length=0.3 )
        T_mover = move_arc( 1, -np.pi/3, T_mover )
        center_x, center_y, angle_start, angle_stop = move_arc_param( 1, -np.pi/3, T_mover )
        print "[ move_arc_param ]", center_x, center_y, angle_start, angle_stop
        plot_2d_transform_T( T_mover, length=0.3 )
        T_mover = move_straight( 1, T_mover )
        plot_2d_transform_T( T_mover, length=0.3 )
    plot_show()
