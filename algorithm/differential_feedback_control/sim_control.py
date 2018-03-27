#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
from coord_render import *

import datetime
import sys

arr = np.array
pi = np.pi
norm = np.linalg.norm
inv = np.linalg.inv

#def get_2d_rot_mat(x,y,t):
#    return np.array([[ np.cos(t), -np.sin(t), x ],
#                     [ np.sin(t),  np.cos(t), y ],
#                     [ 0        ,  0        , 1 ]])

def range_enum(mode,a,b):
    assert len(mode)==2
    assert mode[0] in ["(","["]
    assert mode[1] in [")","]"]
    left_close = (mode[0] == "[" )
    right_close = (mode[1] == "]" )
    assert a<=b
    if(a==b):
        assert (left_close or right_close) 
    def range_enum__(x):
        left_check = False
        right_check = False
        if(left_close):
            left_check = (x>=a)
        else:
            left_check = (x>a)
        if not left_check:
            return -1
        if(right_close):
            right_check = (x<=b)
        else:
            right_check = (x<b)
        if not right_check:
            return 1
        return 0
        #return (left_check and right_check)
    return range_enum__
    
def is_in_range(mode,a,b):
    def is_in_range__(x):
        return not range_enum(mode,a,b)(x)
    return is_in_range__



def control_path(target_point, start_point):
    print "param: ",target_point, start_point
    pi_2 = pi/2
    sim_step = 0.01

    #target_point = arr([0,0,pi/2], dtype=float)
    #start_point = arr([50*np.cos(_t_), 50*np.sin(_t_), _angle_], dtype=float)
    #print "[ target_point ] %s"%target_point
    #print "[ start_point ] %s"%start_point

    reach_threshold_r = 0.01
    reach_threshold_b = 0.01
    
    #k_rou = 3
    #k_alpha = 8
    #k_beta = -1.5

    k_rou = 3
    k_alpha = 8
    k_beta = -1.5


    #print arr([target_point,
    #           start_point])

    T_target = get_2d_transform_mat(*target_point)
    T_start  = get_2d_transform_mat(*start_point)

    #plot_2d_transform_T(T_target, length=5)
    #plot_2d_transform_T(T_start, length=5)
    plot_arrow_v_d(T_target[:2,2], T_target[:2,0]*10, color="red", width=.5, head_width=.5)
    plot_arrow_v_d(T_start[:2,2] , T_start[:2,0]*10 , color="red", width=.5, head_width=1)

    T_start_to_target = np.dot(T_start, inv(T_target))

    s = start_point
    t = target_point
    p = s

    i = 0

    record_array = []
    control_array = []


    diff = t-p
    dx, dy = diff[:2]
    rou = norm(diff[:2])
    diff_theta = mod_angle(-diff[2])

    print diff, diff_theta

    direction = None



    diff = t-p
    dx, dy = diff[:2]
    rou = norm(diff[:2])
    diff_theta = -diff[2] #mod_angle(-diff[2])
            
    theta = mod_angle(diff_theta)
    alpha = mod_angle(-p[2]+np.arctan2(dy,dx))
    beta = mod_angle(-theta-alpha)

    is_front = is_in_range("(]",-pi_2,pi_2)(alpha)
    if(not is_front):
        direction = -1
    else:
        direction = 1

    if(direction == -1):
        p[2] = mod_angle(p[2]+np.pi)
        t[2] = mod_angle(t[2]+np.pi)

    while(True):
        #print p
        record_array.append(p)

        T_mover = get_2d_transform_mat( *p )
        #plot_2d_transform_T( T_mover, length=1 )


        diff = t-p
        dx, dy = diff[:2]
        rou = norm(diff[:2])
        diff_theta = -diff[2] #mod_angle(-diff[2])
        
        if( rou < reach_threshold_r and
            np.abs(mod_angle(diff_theta)) < reach_threshold_b ):
            #print "Goal reached."
            #print arr([t,p])
            #print "\n"*3
            #rec_arr = np.array(record_array)
            #ctrl_arr = np.array(control_array)
            #plt.plot(rec_arr[0], rec_arr[1], color="blue")
            #plot_show()
            break

        
        #beta  = mod_angle(-np.arctan2(dy,dx)+t[2])
        #theta = diff_theta
        #alpha = mod_angle(-beta-theta)

        theta = mod_angle(diff_theta)
        alpha = mod_angle(-p[2]+np.arctan2(dy,dx))
        beta = mod_angle(-theta-alpha)
        

        #print alpha

        is_front = is_in_range("(]",-pi_2,pi_2)(alpha)
        
        #theta = -diff_theta
        #alpha = mod_angle(-theta+np.arctan2(dy,dx))
        #beta  = mod_angle(-theta-alpha)

        #print "dx[ %s ]  dy[ %s ]  "%(dx,dy)
        #print "rou[ %s ]  alpha[ %s ]  beta[ %s ]  "%(rou,alpha,beta)
        #print "k_rou[ %s ]  k_alpha[ %s ]  k_beta[ %s ]  "%(k_rou,k_alpha,k_beta)

        v = k_rou * rou
        #if(direction == -1):
        #    v = -v

        #if(direction is None):
        #    direction = 1 if(is_front) else 0
        #else:
        #    if(direction):
        #        if(not is_front):
        #            alpha = mod_angle(alpha + pi_2)                    
        #    else:
        #        if(is_front):
        #            alpha = mod_angle(alpha + pi_2)


        if(not is_front):
            v = -v

        w = k_alpha*alpha + k_beta * beta


        print "% .7f  % .7f  % .7f"%(theta,v,w)
        #print v


        #if(direction is None):
        #    if(not is_front):
        #        v = -v
        #    direction = 1 if v>0 else -1
        #else:
        #    if (v * direction > 0):
        #        v = -v
            

        #print v


        sin_alpha = np.sin(alpha)
        cos_alpha = np.cos(alpha)


        if(False):
            if(is_front):
                d_rou =  -k_rou * rou * cos_alpha
                d_alpha = k_rou * sin_alpha - k_alpha * alpha - k_beta * beta
                d_beta = -k_rou * sin_alpha
            else:
                d_rou =  k_rou * rou * cos_alpha
                d_alpha = -k_rou * sin_alpha + k_alpha * alpha + k_beta * beta
                d_beta = k_rou * sin_alpha

        plot_step = 30
        if(w):
            r = v/w
            arc_theta = w*sim_step
            if(v>0 and w>0):
                # left front
                pass
            elif(v>0 and w<0):
                # right front
                pass
            elif(v<0 and w>0):
                # right back
                pass
            elif(v<0 and w<0):
                # left back
                pass
            else:
                raise AssertionError("wtf")

            #print r,v,w
            #plot_2d_transform_T(T_mover)
            center_x, center_y, angle_start, angle_stop = move_arc_param( r, arc_theta, T_mover ) # for plotting

            T_mover = move_arc( r, arc_theta, T_mover )
            #plot_2d_transform_T(T_mover)

            #plot_t = np.linspace(angle_start, angle_stop, plot_step)
            #plot_x = center_x + r * np.cos(plot_t)
            #plot_y = center_y + r * np.sin(plot_t)
            #plt.plot(plot_x, plot_y, color='blue')
            control_array.append([v,w,r,w*sim_step])
            #print "[ r arc_theta ] % .3f  %.3f "%(r, arc_theta) 
        else:
            plot_start_x = T_mover[0][2]
            plot_start_y = T_mover[1][2]
            T_mover = move_straight( v*sim_step, T_mover )
            plot_end_x = T_mover[0][2]
            plot_end_y = T_mover[1][2]
            
            plot_t = np.linspace(0,1,plot_step)
            plot_x = plot_start_x*(1-plot_t) + plot_end_x*(plot_t)
            plot_y = plot_start_y*(1-plot_t) + plot_end_y*(plot_t)
            plt.plot(plot_x, plot_y,color='blue')
            control_array.append([v,w,0,v*sim_step])
        p = transform_mat_2d_to_xyt( T_mover )
        #print "[ p from T_mover ]", p
        #print "[ T_mover ]"
        #print T_mover
        #plot_2d_transform_T(get_2d_transform_mat( *p ))
        #plot_show()

        #print "%4d [ p ] % .3f  % .3f  % .3f    [ v w ] % .3f  % .3f   [ rou a b ] % .3f  % .3f  % .3f"%((i,)+tuple(p)+(v,w,rou,alpha,beta))
        #break

        #print "\n"
        i += 1
        if(0):
            break
        if(0):
            if(i> 0):
                #print np.array(record_array)
                ra = np.array(record_array).T
                plt.plot(ra[0], ra[1], color='blue')
                plot_show()
                break



def control_path_on_circle(_t_, _angle_):
    target_point = arr([0,0,pi/2], dtype=float)
    start_point = arr([50*np.cos(_t_), 50*np.sin(_t_), _angle_], dtype=float)
    control_path(target_point, start_point)

    
if __name__ == "__main__":
    #_ts_ = np.linspace(0, 2*np.pi, 13)
    #for i in _ts_:
    #    control_path_on_circle(np.pi/4, i)

    #_ts_ = np.linspace(0, 2*np.pi, 13)
    #for i in _ts_:
    #    control_path_on_circle(i, i)

    #control_path_on_circle(0, 0)
    #_ts_ = np.linspace(0, 2*np.pi, 4)

    num = 8
    step = np.pi * 2 / num
    t__ = [i for i in xrange( num )]
    _ts_ = map(lambda x:x*step, t__)

    #for i in xrange(num):
    #    for j in xrange(num):
    #        params_ = (_ts_[i], _ts_[j])
    #        sys.stdout.write("[ %d | %d %d ] ( % .5f pi, % .5f pi )"%((num,i+1,j+1)+tuple(map(lambda x:x/np.pi,params_))))
    #        control_path_on_circle(*params_)
    #        print " [ OK ]"
            
    #for i in _ts_:
    #    params_ = (_ts_[2],i)
    #    sys.stdout.write("( % .5f pi, % .5f pi )"%tuple(map(lambda x:x/np.pi,params_)))
    #    control_path_on_circle(*params_)
    #    print " [ OK ]"

    #for k in xrange(num):
    #    i = mod_angle(_ts_[k]+np.pi*2/num/2)
    #    params_ = (i,0 if is_in_range("[]",-np.pi/2,np.pi/2)(i) else np.pi)
    #    print "%.3f"%(i/np.pi), is_in_range("[]",-np.pi/2,np.pi/2)(i), range_enum("[]",-np.pi/2,np.pi/2)(i)
    #sys.exit(0)

    

    for k in xrange(num):
        #i = mod_angle(_ts_[k]+np.pi*2/num/2)
        i = mod_angle(_ts_[k])

        #params_ = (i,0 if is_in_range("[]",-np.pi/2,np.pi/2)(i) else np.pi)
        #sys.stdout.write("[ %d | %d ] ( % .5f pi, % .5f pi )"%((num,k)+tuple(map(lambda x:x/np.pi,params_))))
        #control_path_on_circle(*params_)
        #print " [ OK ]"
        #params_ = (i,np.pi if is_in_range("[]",-np.pi/2,np.pi/2)(i) else 0)
        #sys.stdout.write("[ %d | %d ] ( % .5f pi, % .5f pi )"%((num,k)+tuple(map(lambda x:x/np.pi,params_))))
        #control_path_on_circle(*params_)
        #print " [ OK ]"

        center_point = arr([0,0,-np.pi], dtype=float)
        circle_point = arr([50*np.cos(i), 50*np.sin(i), 0], dtype=float)
        params = [center_point, circle_point]
        #lower_bound = -np.pi/4
        #upper_bound = np.pi/4*3
        lower_bound = -np.pi/2
        upper_bound =  np.pi/2
        if(False and is_in_range("(]", lower_bound, upper_bound)(i)):
            control_path(*params[::-1])
            sys.stdout.write("[ %d | %d ] ( % .5f pi, % .5f pi )\n"%((num,k)+tuple(map(lambda x:x/np.pi,circle_point[:2]))))
            sys.stdout.flush()
        else:
            control_path(*params[::1])
            sys.stdout.write("[ %d | %d ] ( % .5f pi, % .5f pi )\n"%((num,k)+tuple(map(lambda x:x/np.pi,circle_point[:2]))))
            sys.stdout.flush()
    

    plot_savefig("figure_sim_control_%s.svg"%(datetime.datetime.now().strftime("%Y%m%d-%H%M%S")),format="svg")
    plot_show()
