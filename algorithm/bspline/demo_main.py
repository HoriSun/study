#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt


class BezierTest(object):
    def __init__(self, degree=3, ctrl_points=np.array([])):
        assert isinstance(degree, int)
        assert degree > 0
        assert len(ctrl_points.shape) == 2
        assert ctrl_points.shape[1] == 2
        assert ctrl_points.shape[0] == degree+1  # 2个点一阶，3个点2阶

        self.degree = degree
        self.ctrl_points = ctrl_points

    def calc_point(self, t):
        cache_ = {}
        def calc_point_inner(t_, ipoint=0, deg=0):
            if deg == self.degree:
                return self.ctrl_points[ipoint]
            pass
            if ipoint in cache_:
                if deg in cache_[ipoint]:
                    return cache_[ipoint][deg]
                pass
            pass

            ret = (t_ * calc_point_inner(t_, ipoint, deg + 1) +
                   (1 - t_) * calc_point_inner(t_, ipoint + 1, deg + 1))

            if ipoint not in cache_:
                cache_[ipoint] = {}
            pass

            cache_[ipoint][deg] = ret

            return ret
        pass

        return calc_point_inner(t, 0, 0)




class UniformBSplineTest(object):
    def __init__(self, degree=3, ctrl_points=np.array([]), clamped=True):
        assert isinstance(clamped, bool)
        assert isinstance(degree, int)
        assert degree > 0
        assert len(ctrl_points.shape) == 2
        assert ctrl_points.shape[1] == 2
        assert ctrl_points.shape[0] > degree  # 至少2个点一阶，3个点2阶
        n_ctrl_points_ = ctrl_points.shape[0]

        if not clamped:
            n_knots_ = degree + n_ctrl_points_ + 1

            # 内部结点（Internal Knots）: n_knots[degree+1: -degree-1] 如果没有内部结点，B样条将退化为贝赛尔曲线。
            n_internal_knots_ = n_knots_ - (degree + 1) * 2
            have_internal_knots_ = n_internal_knots_ > 0

            knot_vector_ = np.linspace(0, 1, n_knots_)
        else:
            n_knots_ = degree + n_ctrl_points_ + 1 + 2

            # 内部结点（Internal Knots）: n_knots[degree+1: -degree-1] 如果没有内部结点，B样条将退化为贝赛尔曲线。
            n_internal_knots_ = n_knots_ - (degree + 1) * 2
            have_internal_knots_ = n_internal_knots_ > 0

            knot_vector_ = np.zeros(n_knots_)
            if have_internal_knots_:
                knot_vector_[degree + 1:-degree - 1] = np.linspace(0, 1, n_internal_knots_ + 2)[1:-1]
                knot_vector_[-degree - 1:] = 1
            else:
                knot_vector_[(n_knots_//2)+(n_knots_%2):] = 1
            pass
        pass

        self.degree = degree
        self.n_ctrl_points = n_ctrl_points_
        self.n_knots = n_knots_
        self.n_internal_knots = n_internal_knots_
        self.have_internal_knots = have_internal_knots_
        self.knot_vector = knot_vector_
        self.ctrl_points = ctrl_points
        self.is_clamped = clamped

        self.basis_function_table = dict()
        self.generate_basis_function_table()

        pass

    def basis_function_degree_0(self, t, iknot):
        # print ("degree[0] iknot[%d] t[%f]"
        #       # " knots=%s"
        #       "" % (iknot, t
        #             # self.knot_vector
        #             ))
        if iknot < (self.n_knots - 1):
            # print ("                            self.knot_vector[iknot]=%f, self.knot_vector[iknot+1]=%f"
            #        "" % (self.knot_vector[iknot], self.knot_vector[iknot + 1]))
            ret = 1 if ((self.knot_vector[iknot] <= t < self.knot_vector[iknot + 1]) and
                        (self.knot_vector[iknot] < self.knot_vector[iknot + 1])) else 0
        else:
            # print ("                            self.knot_vector[iknot]=%f, the last one"
            #        "" % (self.knot_vector[iknot]))
            ret = 0
        pass
        if 0:
            print ("degree[0] iknot[%d] t[%f] ret[%d]"
                   # " knots=%s"
                   "" % (iknot, t, ret
                         # self.knot_vector
                         ))
        elif 0:
            print ("                            ret[%d]"
                   # " knots=%s"
                   "" % (ret
                         # self.knot_vector
                         ))
        pass
        return ret
    
    def basis_function_degree(self, t, iknot, degree):
        if (iknot < (self.n_knots - 1 - degree)) and \
                (self.knot_vector[iknot + degree] != self.knot_vector[iknot]) and \
                (self.knot_vector[iknot + degree + 1] != self.knot_vector[iknot + 1]):
            ret = \
                ((t - self.knot_vector[iknot]) /
                 (self.knot_vector[iknot + degree] - self.knot_vector[iknot]) *
                 (self.basis_function_table[iknot][degree - 1](t))
                 ) \
                 + \
                ((self.knot_vector[iknot + degree + 1] - t) /
                 (self.knot_vector[iknot + degree + 1] - self.knot_vector[iknot + 1]) *
                 (self.basis_function_table[iknot + 1][degree - 1](t))
                 )
        else:
            ret = 0
        pass
        return ret

    def generate_basis_function_table(self):
        """
        基函数表有两个维度： 1.结点ID∈[0,self.n_knots-1]； 2.阶数∈[0,self.degree]
        :return:
        """
        for iknot_ in range(self.n_knots):
            # print ("======================= degree[0] iknot_[%d]" % iknot_)
            self.basis_function_table[iknot_] = {}
            self.basis_function_table[iknot_][0] = (lambda ik_: lambda t: self.basis_function_degree_0(t, ik_))(iknot_)
            # self.basis_function_table[iknot_][0](0.5)
            # for ik_ in self.basis_function_table:
            #    print(self.basis_function_table[ik_][0](0.5))
            # pass
        pass

        if 1:
            for d_ in range(1, self.degree+1):
                for iknot_ in range(self.n_knots):
                    self.basis_function_table[iknot_][d_] = \
                        (lambda ik_, deg_: lambda t: self.basis_function_degree(t, ik_, deg_))(iknot_, d_)
                pass
            pass
        pass

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return ("\n"
                "degree=%d\n"
                "ctrl_points[%d]=%s\n"
                "knots[%d]=%s\n"
                "n_internal_knots[%s]=%d\n"
                "is_clamped=%s\n"
                "" % (self.degree,
                      self.n_ctrl_points, self.ctrl_points,
                      self.n_knots, self.knot_vector,
                      self.have_internal_knots, self.n_internal_knots,
                      self.is_clamped))


if __name__ == '__main__':

    ctrl_points_0 = np.array([[-1, 1],
                              [-1, 0],
                              [-1, -1],
                              [1, 1],
                              [1, 0, ],
                              [1, -1]])

    ctrl_points_1 = np.array([[-2, -1],
                              [-1, 0],
                              [0, 1],
                              [0, -1],
                              [1, 0],
                              [2, 1]])

    ctrl_points_2 = np.array([[-1, 0],
                              [0, 1],
                              [0, -1],
                              [1, 0]])

    ctrl_points_3 = np.array([[-1, -2],
                              [-1, -1],
                              [-1,  0],
                              [ 1,  0],
                              [ 1,  1],
                              [ 1,  2]])

    ctrl_points_4 = np.array([[-1,  2],
                              [-1,  1],
                              [-1,  0],
                              [ 1,  0],
                              [ 1,  1],
                              [ 1,  2]])

    ctrl_points_5 = np.array([[ 1, 0],
                              [ 0, 0],
                              [-1, 0],
                              [ 0,-1],
                              [ 0, 0],
                              [ 0, 1]])

    ctrl_points_6 = np.array([[1, 0],
                              [0, 1],
                              [-1, 0],
                              [0, -1],
                              [1, 0],
                              [0, 1],
                              [-1,0]])

    ctrl_points_7 = np.array([[0, 0],
                              [1, 0],
                              [0, 1],
                              [0, 0],
                              [1, 0],
                              [0, 1]])

    ctrl_points = ctrl_points_6
    pic_name_ = "bspline_test_6"

    bs_ = UniformBSplineTest(degree=3,
                             ctrl_points=ctrl_points,
                             clamped=False)

    # bzt_ = BezierTest(degree=3,
    #                  ctrl_points=ctrl_points[1:-1])

    print(bs_)
    print(bs_.basis_function_table.keys())
    # for iknot_ in bs_.basis_function_table.keys():
    #    print("[%d] %s" % (iknot_, bs_.basis_function_table[iknot_].items()))
    # pass

    npoints_ = 1000

    t_start_ = 0
    t_end_ = 1

    t_start_ = bs_.knot_vector[bs_.degree]
    t_end_ = bs_.knot_vector[-bs_.degree-1]
    print("t_start_[%f] t_end_[%f]" % (t_start_, t_end_))
    t_ = np.linspace(t_start_, t_end_, npoints_)
    # print(t_)
    plot_table_ = {}
    color_map_ = {0:"red",
                  1:"green",
                  2:"blue",
                  3:"purple",
                  4:"black"}
    c_ = np.zeros((2, npoints_))
    p_ = np.zeros(npoints_)
    for degree_ in [3]: # range(bs_.degree+1):
        #for iknot_ in [0,5,6,7,8,9]:  # range(bs_.n_knots):
        for iknot_ in range(bs_.n_knots):
            # print("----iknot=%d" % iknot_)
            v_ = np.array(list(map(bs_.basis_function_table[iknot_][degree_], t_)))
            # plt.plot(t_, v_, color=color_map_[degree_])
            p_ += v_
            plt.plot(t_, v_)
        pass
    pass

    plt.show()

    plt.plot(t_, p_)
    plt.show()

    for ipoint_ in range(bs_.n_ctrl_points):
        v_ = np.array(list(map(bs_.basis_function_table[ipoint_][3], t_)))
        c_[0] += v_ * bs_.ctrl_points[ipoint_][0]
        c_[1] += v_ * bs_.ctrl_points[ipoint_][1]
    pass


    plt.plot(t_, c_[0], color="blue")
    plt.plot(t_, c_[1], color="green")
    # plt.scatter(bs_.ctrl_points.T[0], color="orange")
    #plt.plot(bs_.ctrl_points.T[0], color="red")
    #plt.axis("equal")
    plt.show()

    t_full_ = np.linspace(0, 1, npoints_)
    # bz_curve_ = np.array(list(map(bzt_.calc_point, t_full_))).T

    if 1:
        bspline_curve_plot_ = plt.plot(c_[0], c_[1], color="blue", label="b-spline")
        # bezier_curve_plot_ = plt.plot(bz_curve_[0], bz_curve_[1], color="green", label="bezier")
        plt.scatter(*bs_.ctrl_points.T, color="orange")
        plt.plot(*bs_.ctrl_points.T, color="red", linestyle="-.", alpha=0.5)
        plt.axis("equal")
        #plt.legend(["b-spline", "bezier"])
        plt.legend(["b-spline"])
        plt.savefig("%s.svg" % pic_name_)
        plt.savefig("%s.png" % pic_name_, dpi=300)
        plt.show()

    pass
