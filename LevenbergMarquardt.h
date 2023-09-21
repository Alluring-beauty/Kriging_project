#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <chrono>
#include "InitArrPointer.h"

/*
    LevenbergMarquardt具体实现优化方程类
 */
class LevenbergMarquardt
{
public:
    // 构造函数
    LevenbergMarquardt(double *a, double *b, double *c) : a_(a), b_(b), c_(c)
    {
        epsilon_1_ = 1e-6;
        epsilon_2_ = 1e-6;
        max_iter_ = 50;
        is_out_ = true;
    }
    // 设置参数函数
    void setParameters(double epsilon_1, double epsilon_2, int max_iter, bool is_out)
    {
        epsilon_1_ = epsilon_1;
        epsilon_2_ = epsilon_2;
        max_iter_ = max_iter;
        is_out_ = is_out;
    }
    // 将数据添加入观测值中
    void addObservation(const double &x, const double &y)
    {
        obs_.push_back(Eigen::Vector2d(x, y));
    }
    // 计算要优化曲线方程的雅可比矩阵
    void calcJ_fx()
    {
        J_.resize(obs_.size(), 3);
        fx_.resize(obs_.size(), 1);

        for (size_t i = 0; i < obs_.size(); i++)
        {
            const Eigen::Vector2d &ob = obs_.at(i);
            const double &x = ob(0);
            const double &y = ob(1);
            double j1 = -1;
            double j2 = exp(-1 * x * x / (*c_ * *c_)) - 1;
            double j3 = 2 * *b_ * x * x / (*c_ * *c_ * *c_) * exp(-1 * x * x / (*c_ * *c_));
            J_(i, 0) = j1;
            J_(i, 1) = j2;
            J_(i, 2) = j3;
            fx_(i, 0) = y - (*a_ + *b_ * (1 - exp(-1 * x * x / (*c_ * *c_))));
        }
    }
    // 计算黑塞矩阵
    void calcH_g()
    {
        H_ = J_.transpose() * J_;
        g_ = -J_.transpose() * fx_;
    }

    double getCost()
    {
        Eigen::MatrixXd cost = fx_.transpose() * fx_;
        return cost(0, 0);
    }

    double F(double a, double b, double c)
    {
        Eigen::MatrixXd fx;
        fx.resize(obs_.size(), 1);

        for (size_t i = 0; i < obs_.size(); i++)
        {
            const Eigen::Vector2d &ob = obs_.at(i);
            const double &x = ob(0);
            const double &y = ob(1);
            // change
            // fx(i, 0) = y - exp(a * x * x + b * x + c);
            fx(i, 0) = y - (a + b * (1 - exp(-1 * x * x / c / c)));
        }
        Eigen::MatrixXd F = 0.5 * fx.transpose() * fx;
        return F(0, 0);
    }

    double L0_L(Eigen::Vector3d &h)
    {
        Eigen::MatrixXd L = -h.transpose() * J_.transpose() * fx_ - 0.5 * h.transpose() * J_.transpose() * J_ * h;
        return L(0, 0);
    }
    // 迭代优化过程
    double **solve()
    {
        int k = 0;
        double nu = 2.0;
        calcJ_fx();
        calcH_g();
        bool found = (g_.lpNorm<Eigen::Infinity>() < epsilon_1_);

        std::vector<double> A;
        A.push_back(H_(0, 0));
        A.push_back(H_(1, 1));
        A.push_back(H_(2, 2));
        auto max_p = std::max_element(A.begin(), A.end());
        double mu = *max_p;

        double sumt = 0;

        while (!found && k < max_iter_)
        {
            k = k + 1;
            Eigen::Matrix3d G = H_ + mu * Eigen::Matrix3d::Identity();
            Eigen::Vector3d h = G.ldlt().solve(g_);

            if (h.norm() <= epsilon_2_ * (sqrt(*a_ * *a_ + *b_ * *b_ + *c_ * *c_) + epsilon_2_))
                found = true;
            else
            {
                double na = *a_ + h(0);
                double nb = *b_ + h(1);
                double nc = *c_ + h(2);

                double rho = (F(*a_, *b_, *c_) - F(na, nb, nc)) / L0_L(h);

                if (rho > 0)
                {
                    *a_ = na;
                    *b_ = nb;
                    *c_ = nc;
                    calcJ_fx();
                    calcH_g();

                    found = (g_.lpNorm<Eigen::Infinity>() < epsilon_1_);
                    mu = mu * std::max<double>(0.33, 1 - std::pow(2 * rho - 1, 3));
                    nu = 2.0;
                }
                else
                {
                    mu = mu * nu;
                    nu = 2 * nu;
                }
            }

            if (is_out_)
            {
                result[0][0] = *a_;
                result[1][0] = *b_;
                result[2][0] = *c_;
                std::cout << "Iter: " << std::left << k << " Result: "
                          << std::left << *a_ << " " << std::left << *b_ << " " << std::left
                          << *c_ << " step: " << std::left << h.norm() << " cost: " << std::left
                          << getCost() << endl;
            }
        }
        if (found == true)
        {
            std::cout << "\nFinished!\n\n";
            return result;
        }
        else
        {
            result[0][0] = 0;
            result[1][0] = 0;
            result[2][0] = 0;
            std::cout << "\nDiverged\n\n";
            return result;
        }
    }

    Eigen::MatrixXd fx_;
    Eigen::MatrixXd J_; // 雅克比矩阵
    Eigen::Matrix3d H_; // H矩阵
    Eigen::Vector3d g_;

    std::vector<Eigen::Vector2d> obs_; // 观测

    /* 要求的三个参数 */
    double *a_, *b_, *c_;

    /* parameters */
    double epsilon_1_, epsilon_2_;
    int max_iter_;
    bool is_out_;
    double **result = InitArrPointer(3, 1);
}; // class LevenbergMarquardt