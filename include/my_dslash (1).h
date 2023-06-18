#ifndef _MY_DSLASH_H
#define _MY_DSLASH_H
#include "my_vector.h"
#include <iostream>
class Dslash
{
public:
    Dslash(ComplexVector &U, int &lattice_x, int &lattice_t, int &lattice_spin, double &mass, bool &dag) : U_(U), lattice_x0(lattice_x), lattice_t_(lattice_t), lattice_spin_(lattice_spin), mass_(mass), dag_(dag)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &node_size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &node_rank_);
        lattice_x_ = lattice_x0 / node_size_;
        a_ = 2.0;
        Half_ = 0.5;
        flag_ = (dag_ == true) ? -1 : 1;
        i_ = std::complex<double>(0.0, 1.0);
        backward_rank_ = (node_rank_ + node_size_ - 1) % node_size_;
        forward_rank_ = (node_rank_ + 1) % node_size_;
        backward_send_vec_ = new std::complex<double>[lattice_t_];
        backward_recv_vec_ = new std::complex<double>[lattice_t_];
        forward_send_vec_ = new std::complex<double>[lattice_t_];
        forward_recv_vec_ = new std::complex<double>[lattice_t_];
    }

    ~Dslash()
    {
        delete[] backward_send_vec_;
        delete[] backward_recv_vec_;
        delete[] forward_send_vec_;
        delete[] forward_recv_vec_;
    }
    void dslash(ComplexVector &src, ComplexVector &dest)
    {
        dest.clean_0();
        if (0)
        {
            dslash_test(src, dest);
        }
        else if (node_size_ == 1 && 1)
        {
            dslash_old(src, dest);
        }
        else
        {
            MPI_Request backward_recv_req;
            MPI_Request forward_recv_req;
            pack_up(src, dest);

            MPI_Irecv(forward_recv_vec_, lattice_t_, MPI_DOUBLE_COMPLEX, forward_rank_, node_rank_ * 2 + 0, MPI_COMM_WORLD, &forward_recv_req);
            MPI_Irecv(backward_recv_vec_, lattice_t_, MPI_DOUBLE_COMPLEX, backward_rank_, node_rank_ * 2 + 1, MPI_COMM_WORLD, &backward_recv_req);

            for (int x = 0; x < lattice_x_; x++)
            {
                dslash_no_comm(src, dest, x);
            }

            MPI_Wait(&forward_recv_req, MPI_STATUS_IGNORE);
            // std::cout << "forward_recv:Rank# " << node_rank_ << "<-Rank# " << forward_rank_ << std::endl;

            MPI_Wait(&backward_recv_req, MPI_STATUS_IGNORE);
            // std::cout << "backward_recv:Rank# " << node_rank_ << "<-Rank# " << backward_rank_ << std::endl;

            pack_down(dest);
            std::cout << "forward_recv:Data0# " << forward_recv_vec_[0] << "-DataE# " << forward_recv_vec_[lattice_x_ - 1] << std::endl;
            std::cout << "backward_recv:Data0# " << backward_recv_vec_[0] << "-DataE# " << backward_recv_vec_[lattice_x_ - 1] << std::endl;
        }
    }

private:
    void dslash_test(ComplexVector &src, ComplexVector &dest)
    {
        dest = src * 0.2;
    }
    void dslash_no_comm(ComplexVector &src, ComplexVector &dest, int &x)
    {
        std::complex<double> tmp;
        int x0 = lattice_x_ * node_rank_ + x;
        for (int t = 0; t < lattice_t_; t++)
        {
            // mass_ term
            for (int s = 0; s < lattice_spin_; s++)
            {
                dest[(x * lattice_t_ + t) * 2 + s] += -(a_ + mass_) * src[(x * lattice_t_ + t) * 2 + s];
            }

            // backward t
            int backward_t_ = (t + lattice_t_ - 1) % lattice_t_;
            tmp = (src[(x * lattice_t_ + t) * 2 + 0] + flag_ * i_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * U_[(x0 * lattice_t_ + backward_t_) * 2 + 1];
            dest[(x * lattice_t_ + backward_t_) * 2 + 0] += tmp;
            dest[(x * lattice_t_ + backward_t_) * 2 + 1] -= flag_ * i_ * tmp;

            // forward t
            int forward_t_ = (t + 1) % lattice_t_;
            tmp = (src[(x * lattice_t_ + t) * 2 + 0] - flag_ * i_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * conj(U_[(x0 * lattice_t_ + t) * 2 + 1]);
            dest[(x * lattice_t_ + forward_t_) * 2 + 0] += tmp;
            dest[(x * lattice_t_ + forward_t_) * 2 + 1] += flag_ * i_ * tmp;

            if (x == 0 || x == lattice_x_ - 1)
            {
                continue;
            }
            // backward x
            int backward_x_ = (x + lattice_x_ - 1) % lattice_x_;
            int backward_x0 = (x0 + lattice_x0 - 1) % lattice_x0;
            tmp = (src[(x * lattice_t_ + t) * 2 + 0] + flag_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * U_[(backward_x0 * lattice_t_ + t) * 2 + 0];
            dest[(backward_x_ * lattice_t_ + t) * 2 + 0] += tmp;
            dest[(backward_x_ * lattice_t_ + t) * 2 + 1] += flag_ * tmp;

            // forward x
            int forward_x_ = (x + 1) % lattice_x_;
            tmp = (src[(x * lattice_t_ + t) * 2 + 0] - flag_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * conj(U_[(x0 * lattice_t_ + t) * 2 + 0]);
            dest[(forward_x_ * lattice_t_ + t) * 2 + 0] += tmp;
            dest[(forward_x_ * lattice_t_ + t) * 2 + 1] -= flag_ * tmp;
        }
    }
    void pack_up(ComplexVector &src, ComplexVector &dest)
    {
        MPI_Request backward_send_req;
        MPI_Request forward_send_req;
        int x = 0;
        int x0 = lattice_x_ * node_rank_ + x;
        int backward_x0 = (x0 + lattice_x0 - 1) % lattice_x0;

        for (int t = 0; t < lattice_t_; t++)
        {
            backward_send_vec_[t] = (src[(x * lattice_t_ + t) * 2 + 0] + flag_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * U_[(backward_x0 * lattice_t_ + t) * 2 + 0];
        }
        MPI_Isend(backward_send_vec_, lattice_t_, MPI_DOUBLE_COMPLEX, backward_rank_, backward_rank_ * 2 + 0, MPI_COMM_WORLD, &backward_send_req);
        MPI_Wait(&backward_send_req, MPI_STATUS_IGNORE);
        // std::cout << "backward_send:Rank# " << node_rank_ << "->Rank# " << backward_rank_ << std::endl;

        x = lattice_x_ - 1;
        x0 = lattice_x_ * node_rank_ + x;

        for (int t = 0; t < lattice_t_; t++)
        {
            forward_send_vec_[t] = (src[(x * lattice_t_ + t) * 2 + 0] - flag_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * conj(U_[(x0 * lattice_t_ + t) * 2 + 0]);
        }
        MPI_Isend(forward_send_vec_, lattice_t_, MPI_DOUBLE_COMPLEX, forward_rank_, forward_rank_ * 2 + 1, MPI_COMM_WORLD, &forward_send_req);
        MPI_Wait(&forward_send_req, MPI_STATUS_IGNORE);
        // std::cout << "forward_send:Rank# " << node_rank_ << "->Rank# " << forward_rank_ << std::endl;
    }
    void pack_down(ComplexVector &dest)
    {
        int x = 0;
        int x0 = lattice_x_ * node_rank_ + x;

        for (int t = 0; t < lattice_t_; t++)
        {
            dest[(x * lattice_t_ + t) * 2 + 0] += backward_recv_vec_[t];
            dest[(x * lattice_t_ + t) * 2 + 1] -= flag_ * backward_recv_vec_[t];
        }

        x = lattice_x_ - 1;
        x0 = lattice_x_ * node_rank_ + x;

        for (int t = 0; t < lattice_t_; t++)
        {
            dest[(x * lattice_t_ + t) * 2 + 0] += forward_recv_vec_[t];
            dest[(x * lattice_t_ + t) * 2 + 1] += flag_ * forward_recv_vec_[t];
        }
    }
    void dslash_old(ComplexVector &src, ComplexVector &dest)
    {
        std::complex<double> tmp;
        for (int x = 0; x < lattice_x_; x++)
        {
            for (int t = 0; t < lattice_t_; t++)
            {
                // mass_ term
                for (int s = 0; s < lattice_spin_; s++)
                {
                    dest[(x * lattice_t_ + t) * 2 + s] += -(a_ + mass_) * src[(x * lattice_t_ + t) * 2 + s];
                }

                // backward x
                int backward_x_ = (x + lattice_x_ - 1) % lattice_x_;
                tmp = (src[(x * lattice_t_ + t) * 2 + 0] + flag_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * U_[(backward_x_ * lattice_t_ + t) * 2 + 0];
                dest[(backward_x_ * lattice_t_ + t) * 2 + 0] += tmp;
                dest[(backward_x_ * lattice_t_ + t) * 2 + 1] += flag_ * tmp;

                // forward x
                int forward_x_ = (x + 1) % lattice_x_;
                tmp = (src[(x * lattice_t_ + t) * 2 + 0] - flag_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * conj(U_[(x * lattice_t_ + t) * 2 + 0]);
                dest[(forward_x_ * lattice_t_ + t) * 2 + 0] += tmp;
                dest[(forward_x_ * lattice_t_ + t) * 2 + 1] -= flag_ * tmp;

                // backward t
                int backward_t_ = (t + lattice_t_ - 1) % lattice_t_;
                tmp = (src[(x * lattice_t_ + t) * 2 + 0] + flag_ * i_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * U_[(x * lattice_t_ + backward_t_) * 2 + 1];
                dest[(x * lattice_t_ + backward_t_) * 2 + 0] += tmp;
                dest[(x * lattice_t_ + backward_t_) * 2 + 1] -= flag_ * i_ * tmp;

                // forward t
                int forward_t_ = (t + 1) % lattice_t_;
                tmp = (src[(x * lattice_t_ + t) * 2 + 0] - flag_ * i_ * src[(x * lattice_t_ + t) * 2 + 1]) * Half_ * conj(U_[(x * lattice_t_ + t) * 2 + 1]);
                dest[(x * lattice_t_ + forward_t_) * 2 + 0] += tmp;
                dest[(x * lattice_t_ + forward_t_) * 2 + 1] += flag_ * i_ * tmp;
            }
        }
    }
    ComplexVector U_;
    int lattice_x_, lattice_t_, lattice_spin_, lattice_x0, node_rank_, node_size_, forward_rank_, backward_rank_;
    double mass_, a_, Half_, flag_;
    bool dag_;
    std::complex<double> i_;
    std::complex<double> *backward_send_vec_;
    std::complex<double> *backward_recv_vec_;
    std::complex<double> *forward_send_vec_;
    std::complex<double> *forward_recv_vec_;
};

#endif