#!/usr/bin/env python
import rospy
import numpy as np
from matplotlib import pyplot as plt
from sensor_msgs.msg import LaserScan
from ransac import RANSAC_Line
import coord_render
print "library imported"


class LidarCalibrator(object):
    def __init__(self, 
                 node_name, scan_topic, 
                 angle_min, angle_max,
                 inliner_ratio_required ):
        rospy.init_node(node_name)
        self.scan_topic = scan_topic
        self.scan_received = False
        self.scan_data = None
        self.scan_cb = self.scan_cb_once
        self.loop = self.loop_once
        self.unit_vectors = None
        self.inliner_ratio_required = inliner_ratio_required
        assert( 0 < self.inliner_ratio_required < 1 )
        self.angle_min = angle_min
        self.angle_max = angle_max
        assert self.angle_max > self.angle_min

        self.sub_scan = rospy.Subscriber( self.scan_topic, 
                                          LaserScan,
                                          self.scan_cb,
                                          queue_size = 1 )
    
        self.max_rin = 0

        self.run()

    def run(self):
        rospy.loginfo("Waiting for scan topic [%s]..."%(self.scan_topic))
        while(not rospy.is_shutdown()):
            if(self.scan_received):
                rospy.loginfo("Scan received!")
                break
            else:
                rospy.sleep(0.1)
        else:
            rospy.loginfo("Shutting down")

        while(not rospy.is_shutdown()):
            self.loop() 
            rospy.sleep(0.1)
        else:
            rospy.loginfo("Shutting down")

    def scan_cb_once(self,msg):
        self.scan_cb_normal(msg)
        self.scan_received = True
    
    def scan_cb_normal(self,msg):
        self.scan_data = msg

    def loop_once(self):
        scan_data = self.scan_data
        self.range_min = scan_data.range_min
        self.range_max = scan_data.range_max
        self.gen_unit_vectors(scan_data.angle_min,
                              scan_data.angle_max,
                              scan_data.angle_increment,
                              len(scan_data.ranges))
        self.proc(scan_data.ranges)
        self.loop = self.loop_normal

    def loop_normal(self):
        self.proc(self.scan_data.ranges)
        
    def range_filter(self, r):
        r = np.array(r)
        r[r==np.inf] = 0
        r[r==np.nan] = 0
        r[r<self.range_min] = 0
        r[r>self.range_max] = 0
        return r

    def proc(self, ranges):
        ranges = ranges[self.index_min:self.index_max]
        ranges = self.range_filter(ranges)
        data = self.convert_cartesian(ranges)
        #plt.plot(data.T[0], data.T[1])
        #plt.scatter(data.T[0], data.T[1])
        #plt.axis("equal")
        #plt.show()

        b1, b2 = self.call_ransac(data)
        
        h, radians = self.get_cali_param(b1, b2)

        plt.plot([min(0,data[:,0].min()), data[:,0].max()],
                 [-h,-h],
                 color="black")

        plt.axis("equal")
        plt.show()
        rospy.signal_shutdown(0)

    def print_proc_gen(self, prefix=""):
        def print_proc(s,indent=0):
            return "\n".join(map(lambda p: prefix + " " + "  "*indent + p,
                                 s.split("\n")))
        return print_proc

    def get_cali_param(self,b1,b2):
        print_proc = self.print_proc_gen("[ cali param ]")
        print print_proc("="*30)
        x = np.array([1,2])
        y = b1 + b2 * x
        p = np.concatenate([[x],[y]]).T
        v = p[1] - p[0]
        z = np.zeros(2)

        n = np.array([-v[1],v[0]])
        n = n / np.linalg.norm(n)

        zp0 = z-p[0]

        h = np.dot(zp0, n)

        plt.plot( [p[0,0], p[0,0]+zp0[0]],
                  [p[0,1], p[0,1]+zp0[1]],
                  color = "green")

        n_x = np.array([1,0])

        mega = 100
        n_x_mega = n_x*mega
        v_mega = v*mega
        n_x_mega_norm = np.linalg.norm(n_x_mega)
        v_mega_norm = np.linalg.norm(v_mega)
        cos_a = np.dot(n_x_mega, v_mega)/(n_x_mega_norm*v_mega_norm)
        radians = np.arccos(cos_a)
        angle = radians / np.pi * 180

        print print_proc("p =")
        print print_proc("%s"%p, indent=2)
        print print_proc("v =")
        print print_proc("%s"%v, indent=2)
        print print_proc("z =")
        print print_proc("%s"%z, indent=2)
        print print_proc("n =")
        print print_proc("%s"%n, indent=2)
        print print_proc("zp0 =")
        print print_proc("%s"%zp0, indent=2)
        print print_proc("-"*30)
        print print_proc("h =")
        print print_proc("%s"%h, indent=2)
        print print_proc("radians =")
        print print_proc("%16.10lf"%radians, indent=2)
        print print_proc("angle =")
        print print_proc("%s"%angle, indent=2)
        
        coord_render.plot_2d_transform_T(np.array([[1,0,0],
                                                   [0,1,0],
                                                   [0,0,1]]),
                                         length = 0.2)
        
        coord_render.plot_arrow_v_d( p[0], n*0.2, 
                                     color='red',
                                     width=0.01,
                                     head_width=0.05)

        coord_render.plot_arrow_v_d( p[0], v*0.2, 
                                     color='green',
                                     width=0.01,
                                     head_width=0.05)

        return h,radians
        

    def gen_unit_vectors(self, amin, amax, da, n):
        print "number of original laser points = %d"%n
        amin_use = amin
        amax_use = amax

        if (amin >= self.angle_min):
            self.index_min = 0
        else:
            self.index_min = int ( ( self.angle_min - amin ) / da )
            amin_use = self.angle_min

        if (amax <= self.angle_max):
            self.index_max = n
        else:
            self.index_max = int ( ( self.angle_max - amin ) / da )
            amax_use = self.angle_max

        print "amin = %lf   amax = %lf"%(amin, amax)
        print "amin_use = %lf   amax_use = %lf"%(amin_use, amax_use)
        print "self.angle_min = %lf   self.angle_max = %lf"%(self.angle_min, self.angle_max)

        n_use = self.index_max - self.index_min
        self.point_use_ratio = n_use/float(n)
        self.inliner_ratio_required /= self.point_use_ratio
        print "laser points used = %d"%(n_use)
    
        angles = np.arange(n_use) * da + amin_use
        ys = np.sin(angles)
        xs = np.cos(angles)
        # not transposition, as multiplication will need it
        # to be an ( N x 2 ) matrix later
        self.unit_vectors = np.concatenate([[xs],[ys]], axis=0)
        
    def convert_cartesian(self, ranges):
        return ( self.unit_vectors * ranges ).T

    def stop_cond(self, it, rin):
        ret = (rin >= self.inliner_ratio_required)
        if(rin > self.max_rin):
            self.max_rin = rin
        print "[%d] %lf %lf"%(it, rin, self.max_rin)
        if(rospy.is_shutdown()):
            ret = True
        if(ret):
            print "total iteration:",it
        return ret 

    def call_ransac(self,data):

        plt.plot(data.T[0], data.T[1])
        plt.scatter(data.T[0], data.T[1])

        e = 0.05

        ymin = data[:,1].min()
        ymax = data[:,1].max()
        xmin = data[:,0].min()
        xmax = data[:,0].max()
        rmin = xmin
        rmax = xmax

        b1, b2, pin = RANSAC_Line( data = data,
                                   nmin = 2,#n * 0.5,
                                   niter = 10000, 
                                   dmax = e, 
                                   ratio_in_min = 0.2,
                                   stop_cond = self.stop_cond ).run()

        #plt.scatter(data[:,0], data[:,1], color="black")
        #plt.scatter(clean_data[:,0], clean_data[:,1], color="red", s=20, linewidth=5)

        #plt.plot(ori_line_x,ori_line_y,color="green",linewidth=2)

        #b1, b2 = LSE(data)
    
        line_x = np.array([rmin,rmax])
        line_y = b1 + b2*line_x


        #if(line_y.min()<ymin or line_y.max()>ymax):
        #    line_y = np.array([ymin,ymax])
        #    line_x = (line_y - b1)/b2
        if(line_y.min()<ymin):
            line_y = np.array([ymin,line_y.max()])
            line_x = (line_y - b1)/b2
        if(line_y.max()>ymax):
            line_y = np.array([line_y.min(),ymax])
            line_x = (line_y - b1)/b2
        
        datain = data[pin]
        plt.scatter(datain[:,0], datain[:,1], color="orange")
        
        plt.plot(line_x, line_y, color="red", linewidth=2)
        
        #plt.show()
        return b1, b2



if __name__ == "__main__":
    lc = LidarCalibrator( node_name = "lidar_calibrator", 
                          scan_topic = "scan_hokuyo",
                          angle_min = -1.57,
                          angle_max = 0,
                          inliner_ratio_required = 0.15)
