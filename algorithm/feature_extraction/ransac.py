#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

print "library imported"

def data_gen(a,b,c,e,n,nnoise,rmin,rmax):
    '''
    % a, b, c:  ( a*x + b*y + c = 0 ), descriptor of target line
    % e: radius of gaussian error 
    % n: number of normal points
    % nnoise: number of noise points
    % rmin, rmax: data range
    '''
    assert(a!=0 or b!=0)
    a = float(a)
    b = float(b)
    c = float(c)

    rd = rmax - rmin

    r = np.random.rand(n) * rd + rmin

    r[-2:] = [rmin, rmax]
    
    if(b == 0):
        y = r
        x = np.zeros(n) + c/a
    elif(a == 0):
        x = r
        y = np.zeros(n) + c/b
    else:
        x = r
        y = ( c - a*x ) / b

    #plt.plot(x[-2:],y[-2:],color="purple",linewidth=3)

    data = np.concatenate([[x],[y]], axis=0).T
    data += (np.random.random([n,2]) * 2 - 1) * e
    
    data_x_min = data[:,0].min()
    data_x_max = data[:,0].max()
    data_y_min = data[:,1].min()
    data_y_max = data[:,1].max()

    data_x_range = data_x_max - data_x_min
    data_y_range = data_y_max - data_y_min

    if(a == 0):
        data_y_min = data_y_min - (data_x_range-data_y_range)/2.0
        data_y_range = data_x_range
    elif(b == 0):
        data_x_min = data_x_min - (data_y_range-data_x_range)/2.0
        data_x_range = data_y_range
    elif(a/b < 1e-1):
        data_y_min = data_y_min - (data_x_range-data_y_range)/2.0
        data_y_range = data_x_range
    elif(b/a < 1e-1):
        data_x_min = data_x_min - (data_y_range-data_x_range)/2.0
        data_x_range = data_y_range
    else:
        pass

    data_noise = np.random.random([nnoise, 2])
    data_noise[:,0] = data_noise[:,0] * data_x_range + data_x_min
    data_noise[:,1] = data_noise[:,1] * data_y_range + data_y_min

    data_with_noise = np.concatenate([data, data_noise])

    return data_with_noise, data


def LSE(p):
    '''
    % p: 2D point array
    '''
    n = p.shape[0]

    x = p[:,0]
    y = p[:,1]

    x2 = x*x
    y2 = y*y
    xy = x*y

    sx = x.sum()
    sy = y.sum()
    sx2 = x2.sum()
    sy2 = y2.sum()
    sxy = xy.sum()

    dom = (   n * sx2 - sx * sx  ) # dominator
    if(dom == 0):
        return 0,0
    b2  = (   n * sxy - sx * sy  ) / dom
    b1  = ( sx2 * sy  - sx * sxy ) / dom

    return b1,b2


def ransac_line(data, nmin, niter, dmax, ratio_in_min):
    '''
    % data: 2D point array
    % nmin: minimum number of points in a fitting set
    % niter: number of iterations
    % dmax: maximum distance between a point to a fitting line
    % ratio_in_min: minimum ratio of inliners in a fitting set
    '''
    assert nmin > 0
    assert 0 < ratio_in_min < 1
    n = data.shape[0]
    n_in_min = round(ratio_in_min * n)
    best_in_num = 0
    d1 = 0
    d2 = 0
    pin = [False]*n
    rin = 0
    #for i in xrange(niter):
    while(rin < 0.25):
        spidx = np.random.randint( low = 0, 
                                   high = n-1, 
                                   size = nmin ) # sample
        #print spidx
        sp = data[spidx,:]

        b1, b2 = LSE(sp)
        kline_x = np.array([0,1])
        kline_y = b1 + b2 * kline_x
        #kline = sp[1] - sp[0]
        #kline = np.concatenate([[kline_x],
        #                        [kline_y]],
        #                        axis = 0)
        kline = np.array([kline_x[1]-kline_x[0],
                          kline_y[1]-kline_y[0]])
        kline_normalized = kline/np.linalg.norm(kline)
        kline_norm_v = np.array([-kline_normalized[1],
                                  kline_normalized[0]]) # 90 degrees
        d_v = data - sp[0]
        d = np.dot( kline_norm_v, d_v.T )
        nin = np.count_nonzero(np.abs(d) <= dmax)
        rin = nin/float(n)
        if ( ( nin >= n_in_min ) and
             ( nin > best_in_num ) ):
            best_in_num = nin
            # [ TODO ] May have pure vertical or pure horizontal line. Handle it.
            d1, d2 = b1, b2
            pin = np.abs(d) <= dmax
            print d1,d2,nin,rin
            pass # if
        #print i
        pass # for
        

    return d1, d2, pin

















class RANSAC_Line(object):
    def __init__(self,data, nmin, niter, dmax, ratio_in_min, 
                 stop_cond = lambda it, rin: rin >= 0.25):
        '''
        % data: 2D point array
        % nmin: minimum number of points in a fitting set
        % niter: number of iterations
        % dmax: maximum distance between a point to a fitting line
        % ratio_in_min: minimum ratio of inliners in a fitting set
        '''

        assert nmin > 0
        assert 0 < ratio_in_min < 1
        
        self.data = data
        self.nmin = nmin
        self.niter = niter
        self.dmax = dmax
        self.ratio_in_min = ratio_in_min
        self.stop_cond = stop_cond

    def run(self):
        self.n = self.data.shape[0]
        self.n_in_min = round(self.ratio_in_min * self.n)
        self.best_in_num = 0
        self.d1 = 0
        self.d2 = 0
        self.b1 = 0
        self.b2 = 0
        self.pin = [False]*self.n # for debugging
        self.rin = 0
        
        it = 0

        while(not self.stop_cond(it, self.rin)):
            it += 1
            spidx = np.random.randint( low = 0, 
                                       high = self.n-1, 
                                       size = self.nmin ) # sample
            sp = self.data[spidx,:]
            self.calc(sp)

            while ( ( self.nin >= self.n_in_min ) and
                    ( self.nin > self.best_in_num ) ):
                self.best_in_num = self.nin
                # [ TODO ] May have pure vertical or pure horizontal line. Handle it.
                self.d1, self.d2 = self.b1, self.b2
                self.pin = np.abs(self.d) <= self.dmax
                print self.d1,self.d2,self.nin,self.rin
                
                sp = self.data[self.pin]

                self.calc(sp)
                
                pass # if

            pass # while
        

        return self.d1, self.d2, self.pin



    def calc(self, sp):
        self.b1, self.b2 = LSE(sp)
        kline_x = np.array([0,1])
        kline_y = self.b1 + self.b2 * kline_x
        #kline = sp[1] - sp[0]
        #kline = np.concatenate([[kline_x],
        #                        [kline_y]],
        #                        axis = 0)
        kline = np.array([kline_x[1]-kline_x[0],
                          kline_y[1]-kline_y[0]])
        kline_normalized = kline/np.linalg.norm(kline)
        kline_norm_v = np.array([-kline_normalized[1],
                                  kline_normalized[0]]) # 90 degrees
        d_v = self.data - sp[0]
        self.d = np.dot( kline_norm_v, d_v.T )
        self.nin = np.count_nonzero(np.abs(self.d) <= self.dmax)
        self.rin = self.nin/float(self.n)













def test(a,b,c,e):
    a = float(a)
    b = float(b)
    c = float(c)
    e = float(e)
    n = 100
    nnoise = 300
    rmin = -10
    rmax = 10

    data, clean_data = data_gen( a = a,
                                 b = b,
                                 c = c,
                                 e = e,
                                 n = n,
                                 nnoise = nnoise,
                                 rmin = rmin,
                                 rmax = rmax )

    ymin = data[:,1].min()
    ymax = data[:,1].max()

    r = np.array([rmin,rmax])

    if(b == 0):
        ori_line_y = r
        ori_line_x = np.zeros(2) + c/a
    elif(a == 0):
        ori_line_x = r
        ori_line_y = np.zeros(2) + c/b
    else:
        ori_line_x = r
        ori_line_y = ( c - a*ori_line_x ) / b


    if(ori_line_y.min()<ymin or ori_line_y.max()>ymax):
        ori_line_y = np.array([ymin,ymax])
        ori_line_x = ( c - b*ori_line_y ) / a

    plt.scatter(data[:,0], data[:,1], color="black")
    #plt.scatter(clean_data[:,0], clean_data[:,1], color="red", s=20, linewidth=5)

    #plt.plot(ori_line_x,ori_line_y,color="green",linewidth=2)

    plt.show()

    def stop_cond(it, rin):
        ret = (rin >= 0.25)
        if(ret):
            print "total iteration:",it
        return ret 

    b1, b2, pin = RANSAC_Line( data = data,
                               nmin = 2,#n * 0.5,
                               niter = 10000, 
                               dmax = e, 
                               ratio_in_min = 0.2,
                               stop_cond = stop_cond ).run()

    plt.scatter(data[:,0], data[:,1], color="black")
    plt.scatter(clean_data[:,0], clean_data[:,1], color="red", s=20, linewidth=5)

    plt.plot(ori_line_x,ori_line_y,color="green",linewidth=2)

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

    plt.plot(line_x, line_y, color="blue", linewidth=2)

    plt.show()
    


if __name__ == "__main__":
    #test(-1,-2,3,0.5)
    #test(0,1,3,0.5)
    #test(1,0,3,0.5)
    #test(0.1,10,3,0.5)
    #test(10,0.1,3,0.5)
    test(0.1,1,3,0.01)
    #test(1,0.1,3,0.5)
