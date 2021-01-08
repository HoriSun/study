#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class UniformBSplineTest(object):
    def __init__(self, degree=3, ctrl_points=np.array([]),
                 knot_vector=np.array([]), clamped=True,
                 precalc_derivative=False):
        assert isinstance(clamped, bool)
        assert isinstance(degree, int)
        assert degree > 0
        assert len(ctrl_points.shape) == 2
        assert ctrl_points.shape[1] == 2
        assert ctrl_points.shape[0] > degree  # 至少2个点一阶，3个点2阶
        n_ctrl_points_ = ctrl_points.shape[0]

        if len(knot_vector) is 0:
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
                    knot_vector_[(n_knots_ // 2) + (n_knots_ % 2):] = 1
                pass
            pass
        else:
            knot_vector_ = knot_vector
            n_knots_ = knot_vector_.shape[0]
            n_internal_knots_ = n_knots_ - (degree + 1) * 2
            have_internal_knots_ = n_internal_knots_ > 0
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

        self.t_start = self.knot_vector[self.degree]
        self.t_end = self.knot_vector[-self.degree - 1]
        self.t_range = self.t_end - self.t_start

        self.__derivative = None
        if precalc_derivative:
            derv_ = self.derivative
        pass

        self.__length = None
        self.__get_length()
        len_ = self.length

        pass

    def basis_function_degree_0(self, t, iknot):
        if iknot < (self.n_knots - 1):
            ret = 1 if ((self.knot_vector[iknot] <= t < self.knot_vector[iknot + 1]) and
                        (self.knot_vector[iknot] < self.knot_vector[iknot + 1])) else 0
        else:
            ret = 0
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
            self.basis_function_table[iknot_] = {}
            self.basis_function_table[iknot_][0] = \
                (lambda ik_: lambda t: self.basis_function_degree_0(t, ik_))(iknot_)
        pass

        if 1:
            for d_ in range(1, self.degree + 1):
                for iknot_ in range(self.n_knots):
                    self.basis_function_table[iknot_][d_] = \
                        (lambda ik_, deg_: lambda t: self.basis_function_degree(t, ik_, deg_))(iknot_, d_)
                pass
            pass
        pass

    def calc_point(self, t):
        return reduce(lambda a,b:a+b,
                      map(lambda ipoint_: ( self.basis_function_table[ipoint_][self.degree](t) *
                                            self.ctrl_points[ipoint_] ),
                          range(self.n_ctrl_points)))

    def calc_path(self, npoints):
        t_ = np.linspace(self.t_start, self.t_end, npoints)
        path_ = np.array(list(map(self.calc_point, t_)))
        return path_

    def get_t_range(self):
        return self.t_start, self.t_end

    def percent_to_t(self, percent):
        return (percent * self.t_range) + self.t_start

    def get_internal_knots(self):
        return self.knot_vector[self.degree:-self.degree][:]

    def get_internal_knot_points(self):
        """
        这些就是所谓的”型值点“。
        :todo: 1. 实现用型值点反算B样条曲线
        :todo: 2. 计算曲率
        :todo: 3. 计算切线
        :todo: 4. 曲率连续/曲率平滑变化
        :return:
        """
        internal_knots_ = self.get_internal_knots()
        return np.array(list(map(self.calc_point, internal_knots_)))

    def __get_derivative(self):
        derv_ctrl_points_ = np.zeros((self.n_ctrl_points-1, 2))
        for i in range(self.n_ctrl_points-1):
            derv_ctrl_points_[i] = \
                ( self.degree /
                  (self.knot_vector[i+self.degree+1] + self.knot_vector[i+1]) *
                  (self.ctrl_points[i+1] - self.ctrl_points[i]) )
        pass

        return UniformBSplineTest(degree=self.degree-1,
                                  ctrl_points=derv_ctrl_points_,
                                  knot_vector=self.knot_vector[1:-1],
                                  clamped=self.is_clamped,
                                  precalc_derivative=False)

    def get_tangent_vector_at(self, percent):
        """
        :param percent: 总长的百分比，范围在 [self.knot_vector[self.degree], self.knot_vector[-self.degree]]
        :return:
        """
        assert 0 <= percent <= 1
        return self.derivative.calc_point(self.percent_to_t(percent))

    def get_unit_tangent_vector_at(self, percent):
        """
        :param percent: 总长的百分比，范围在 [self.knot_vector[self.degree], self.knot_vector[-self.degree]]
        :return:
        """
        tanvec_ = self.get_tangent_vector_at(percent)
        tanvec_ = tanvec_ / np.linalg.norm(tanvec_)
        return tanvec_

    def __get_length(self):
        """
        粗糙的计算方法。直接采离散点算距离相加。
        更好的方法：
        1.　转换成几条贝塞尔曲线，解析地求长度，或用高斯－勒让德数值积分
        2.　根据knot拆成几条多项式曲线，解析地求长度
        :return:
        """
        path_ = self.calc_path(npoints=1000)
        dpath_ = path_[1:] - path_[:-1]
        dlen_ = np.linalg.norm(dpath_, axis=1)
        # print(dlen_.shape)
        return dlen_.sum()

    @property
    def derivative(self):
        if not self.__derivative:
            self.__derivative = self.__get_derivative()
        pass
        return self.__derivative

    @property
    def length(self):
        if not self.__length:
            self.__length = self.__get_length()
        pass
        return self.__length

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
