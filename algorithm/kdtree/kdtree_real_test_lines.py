#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

def plot_show():
    plt.grid(True)
    plt.axes().set_aspect('equal','datalim')
    plt.show()



class kd_tree(object):
    def __init__(self,p):
        assert len(p.shape) == 2
        self.p = p
        self.num = p.shape[0]
        self.dim = p.shape[1]
        self.node_type = self.gen_node_type(self.dim)
        self.tree = self.generate(self.p, 
                                  self.num,
                                  self.dim,
                                  self.node_type)

    @staticmethod
    def gen_node_type(dim):
        class node(object):
            dim_ = dim
            LEFT = 0
            RIGHT = 1
            def __init__(self, point=None, axis=None, value=None, parent=None, side=None):
                if(point is not None):
                    assert len(point.shape) == 1
                    assert point.shape[0] == self.dim_
                self.point = point
                self.value = value 
                self.axis = axis
                self.side = side
                self.left = None
                self.right = None
                self.parent = parent
        return node

    @staticmethod
    def generate(p_, num_, dim_, node_type):
        assert len(p_.shape) == 2
        assert p_.shape[0] == num_
        assert p_.shape[1] == dim_
        # split from the axis with maximum

        point_set_stack = []
        
        def execute(p, parent, side):
            axis_select = np.argmax(np.var(p,axis=0))
            p = p[np.argsort(p[:,axis_select])]

            # select the one in the middle to make the
            # binary tree (BST) balance,
            # instead of do 1-D clustering and attempting
            # to make a multi-spanning unbanlance tree
            num = p.shape[0]

            if(num == 1):
                point = p[0]
                node = node_type(point=point,
                                 axis=None, 
                                 value=None,
                                 parent=parent, 
                                 side=side)
            else:
                select_index_left = num/2+num%2-1
                select_index_right = select_index_left + 1
                split_value = ( p[select_index_left, axis_select] + 
                                p[select_index_right, axis_select] ) / 2.0

                node = node_type(point=None, 
                                 axis=axis_select,
                                 value=split_value,
                                 parent=parent,
                                 side=side)

                point_set_stack.append((p[:select_index_right],
                                        node,
                                        node_type.LEFT))
                point_set_stack.append((p[select_index_right:],
                                        node,
                                        node_type.RIGHT))            
                pass


            if(side == node_type.LEFT):
                parent.left = node
            elif(side == node_type.RIGHT):
                parent.right = node
            elif(side == None):
                pass
            else:
                raise Exception("unexpected behavior 0")
                pass

            return node
        #enddef execute(p, parent, lr)
        
        point_set_stack.append((p_, None, None))
        root = execute(*point_set_stack.pop(0))
        while(point_set_stack):
            execute(*point_set_stack.pop(0))
        
        return root
    #enddef generate(p_, num_, dim_, node_type)

    def plot_2d(self, dim_max_, dim_min_):
        dim_max_p = self.p.max(axis=0)
        dim_min_p = self.p.min(axis=0)
        
        if (dim_max_p > dim_max_).any():
            dim_max_ = dim_max_p
        if (dim_min_p < dim_min_).any():
            dim_min_ = dim_min_p

        node_type = self.node_type
        proc_stack = [(self.tree, dim_max_, dim_min_)]
        
        while(proc_stack):
            p, dim_max, dim_min = proc_stack.pop(0)
            axis = p.axis
            
            #print p.point, np.array(map(lambda x:x[0].point,proc_stack))

            if(p.point is not None):
                continue

            if(axis==0):
                plt.plot([ p.value    , p.value    ], 
                         [ dim_min[1] , dim_max[1] ], 
                         color='orange')
            elif(axis==1):
                plt.plot([ dim_min[0] , dim_max[0] ],
                         [ p.value    , p.value    ], 
                         color='orange')
            else:
                raise Exception("unexpected behavior 1")


            if(p.left):
                if(axis == 0):
                    proc_stack.append((p.left,  
                                       [p.value   , dim_max[1]],
                                       dim_min))
                if(axis == 1):
                    proc_stack.append((p.left,  
                                       [dim_max[0], p.value   ],
                                       dim_min))
                else:
                    pass

            if(p.right):
                if(axis == 0):
                    proc_stack.append((p.right,  
                                       dim_max,
                                       [p.value   , dim_min[1]]))
                if(axis == 1):
                    proc_stack.append((p.right,  
                                       dim_max,
                                       [dim_min[0], p.value   ]))
                else:
                    pass
            
    #enddef 2d_tree_plot(root, dim_max, dim_min):    

    def search_closest(self, p):
        p = np.array(p)
        assert len(p.shape) == 1
        assert p.shape[0] == self.dim

        root = self.tree
        d = [np.inf]
        c = [None]
        norm = np.linalg.norm # alias

        def recursive(n):
            if(n.point is not None):
                dist = norm(n.point-p)
                if(dist < d[0]):
                    d[0] = dist
                    c[0] = n.point
            else:
                p_axis = p[n.axis]
                dist   = p_axis - n.value
                dist_l = p_axis - d[0] - n.value
                dist_r = p_axis + d[0] - n.value

                # left side or middle, left first
                #
                # for middle, it should be OK for both sides,
                # but in my implementation, left side tend to
                # have 1 more point when number of points of 
                # the sub-region is odd number, so left first 
                # here.
                if(dist <= 0):
                    if(dist_l <= 0):
                        recursive(n.left)
                    if(dist_r >= 0):
                        recursive(n.right)
                else: # right first
                    if(dist_r >= 0):
                        recursive(n.right)                    
                    if(dist_l <= 0):
                        recursive(n.left)
                    pass
    
        recursive(root)
        return c[0]

num = 20
p = np.random.random((num,2))
plt.scatter(*p.T, color='cyan')

# plot box
plt.plot([0,1,1,0,0], [0,0,1,1,0], color='black')


kd_tree_p = kd_tree(p)
root = kd_tree_p.tree
kd_tree_p.plot_2d([1]*2, [0]*2)
tp = np.random.random(2)
cp = kd_tree_p.search_closest(tp)
plt.scatter([tp[0]],[tp[1]],color='red')
plt.scatter([cp[0]],[cp[1]],color='green')

def draw_circle(x,y,r):
    t = np.linspace(0, np.pi*2, 30)
    x_ = r * np.cos(t)
    y_ = r * np.sin(t)
    plt.plot(x_+x,y_+y,color='black')

draw_circle(tp[0], tp[1], np.linalg.norm(cp-tp))

plot_show()
