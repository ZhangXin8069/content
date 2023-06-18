#include <complex>
#include <iostream>
#include <mpi.h>

const int N = 100; // 数组的大小

class Example
{
public:
    Example(int rank, int size)
    {
        // 初始化 MPI 环境
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);

        // 为非阻塞通信分配缓冲区
        send_buf_ = new std::complex<double>[N];
        recv_buf_ = new std::complex<double>[N];
    }

    ~Example()
    {
        // 释放缓冲区
        delete[] send_buf_;
        delete[] recv_buf_;
    }

    void Run()
    {
        // 初始化数据
        for (int i = 0; i < N; i++)
        {
            send_buf_[i] = std::complex<double>(rank_, i);
        }

        // 发送数据
        MPI_Request send_request;
        MPI_Isend(send_buf_, N, MPI_DOUBLE_COMPLEX, (rank_ + 1) % size_, 0,
                  MPI_COMM_WORLD, &send_request);

        // 接收数据
        MPI_Request recv_request;
        MPI_Irecv(recv_buf_, N, MPI_DOUBLE_COMPLEX, (rank_ + size_ - 1) % size_, 0,
                  MPI_COMM_WORLD, &recv_request);

        // 等待通信完成
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);

        // 处理接收到的数据
        for (int i = 0; i < N; i++)
        {
            std::cout << "Rank " << rank_ << " received: " << recv_buf_[i] << std::endl;
        }
    }

private:
    int rank_; // 当前进程的编号
    int size_; // 进程总数

    std::complex<double> *send_buf_; // 发送缓冲区
    std::complex<double> *recv_buf_; // 接收缓冲区
};

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Example example(rank, size);
    example.Run();

    MPI_Finalize();
    return 0;
}
