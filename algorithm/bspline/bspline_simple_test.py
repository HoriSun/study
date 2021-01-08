#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import collections as mcollect
from matplotlib import pyplot as plt

from bspline_simple import UniformBSplineTest

if __name__ == "__main__":
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
                              [-1, 0],
                              [1, 0],
                              [1, 1],
                              [1, 2]])

    ctrl_points_4 = np.array([[-1, 2],
                              [-1, 1],
                              [-1, 0],
                              [1, 0],
                              [1, 1],
                              [1, 2]])

    ctrl_points_5 = np.array([[1, 0],
                              [0, 0],
                              [-1, 0],
                              [0, -1],
                              [0, 0],
                              [0, 1]])

    ctrl_points_6 = np.array([[1, 0],
                              [0, 1],
                              [-1, 0],
                              [0, -1],
                              [1, 0],
                              [0, 1],
                              [-1, 0]])

    ctrl_points_7 = np.array([[0, 0],
                              [1, 0],
                              [0, 1],
                              [0, 0],
                              [1, 0],
                              [0, 1]])

    ctrl_points_8 = np.array([[0, 2],
                              [1, 1],
                              [2, 0],
                              [3, 1],
                              [4, 0],
                              [5, 1],
                              [6, 0],
                              [7, 1],
                              [7, 6],
                              [6, 7],
                              [5, 6],
                              [4, 7],
                              [3, 6],
                              [2, 7],
                              [1, 6],
                              [0, 5]])

    ctrl_points_9 = np.array([[-1.9904146, -1.20562073],
                              [-1.7762073, -0.51281037],
                              [-1.562, 0.18],
                              [-1.858, 0.962],
                              [-1.62434735, 1.56556952],
                              [-1.39069469, 2.16913905]]
                             )

    ctrl_points = ctrl_points_9
    pic_name_ = "bspline_test_6"

    bs_ = UniformBSplineTest(degree=3,
                             ctrl_points=ctrl_points,
                             clamped=False)
    path_ = bs_.calc_path(1000).T

    bs_derv_ = bs_.derivative

    derv_point_0_ = bs_.calc_path(100).T
    derv_path_ = bs_derv_.calc_path(100).T
    derv_path_uni_ = derv_path_ / np.linalg.norm(derv_path_, axis=0)
    derv_vec_len_ = 0.1
    derv_point_1_ = derv_point_0_ + derv_path_uni_ * derv_vec_len_
    print(np.arctan2(derv_path_[1], derv_path_[0]))
    # derv_point_1_ = derv_point_0_ + derv_path_
    # print(derv_path_)
    # print(np.linalg.norm(derv_path_, axis=0))
    # print(derv_point_0_.shape)
    # print(derv_point_1_.shape)
    derv_seg_ = np.dstack([derv_point_0_.T, derv_point_1_.T]).transpose((0, 2, 1))
    # print(derv_seg_)

    percent_spec_ = 1.0
    # print(bs_.t_start)
    # print(bs_.t_end)
    # print(bs_.t_range)
    # print(percent_spec_)
    # print(bs_.percent_to_t(percent_spec_))
    derv_vec_spec_ = bs_.get_unit_tangent_vector_at(percent_spec_)
    derv_point_0_spec_ = bs_.calc_point(bs_.percent_to_t(percent_spec_))
    derv_point_1_spec_ = derv_point_0_spec_ + derv_vec_spec_

    internal_knots_ = bs_.get_internal_knot_points().T
    # print(derv_vec_spec_)
    # print(derv_point_0_spec_)
    # print(derv_point_1_spec_)

    # print(bs_.length)

    lc = mcollect.LineCollection(derv_seg_, colors="red", linewidths=2)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    # ax.autoscale()
    # ax.margins(0.1)
    plt.plot(path_[0], path_[1], color="blue")

    plt.scatter(*derv_point_0_, color="cyan")
    plt.scatter(*derv_point_1_, color="purple")

    plt.plot([derv_point_0_spec_[0],
              derv_point_1_spec_[0]],
             [derv_point_0_spec_[1],
              derv_point_1_spec_[1]],
             color="orange")

    plt.plot(*derv_path_uni_, color="green")
    plt.scatter(*internal_knots_, color="orange")
    plt.axis("equal")
    # plt.savefig("temp.png", dpi=300)
    plt.show()
