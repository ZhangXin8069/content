#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
    // 初始化MPI环境
    MPI_Init(&argc, &argv);

    // 获取当前进程的编号和总进程数
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 创建非阻塞通信的句柄
    MPI_Request request;

    // 向所有其它线程群发数据
    for (int i = 0; i < size; i++)
    {
        if (i == rank) // 当前进程不需要发送数据给自己
        {
            continue;
        }

        int data = rank;                                              // 要发送的数据
        MPI_Isend(&data, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request); // 非阻塞发送数据
    }

    // 定义一个接收数据的缓冲区
    vector<int> recv_buf(size);

    // 接收其它线程发来的数据
    for (int i = 0; i < size; i++)
    {
        if (i == rank) // 当前进程不会接收自己的数据
        {
            continue;
        }

        int data;
        MPI_Recv(&data, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // 阻塞接收数据
        recv_buf[i] = data;                                                   // 将接收到的数据存入缓冲区
    }

    // 输出接收到的数据
    cout << "Process " << rank << " received data: ";
    for (int i = 0; i < size; i++)
    {
        cout << recv_buf[i] << " ";
    }
    cout << endl;

    // 释放MPI环境
    MPI_Finalize();

    return 0;
}