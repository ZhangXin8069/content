#include <iostream>
#include <mpi.h>
#include <thread>

void SendData(int rank, int size) {
    int send_data = rank;
    MPI_Request request;
    MPI_Isend(&send_data, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
}

void RecvData(int rank, int size) {
    int recv_data;
    MPI_Request request;
    MPI_Irecv(&recv_data, 1, MPI_INT, (rank + size - 1) % size, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    std::cout << "Process " << rank << " received " << recv_data << " from process " << (rank + size - 1) % size << std::endl;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::thread send_thread(SendData, rank, size);
    std::thread recv_thread(RecvData, rank, size);

    send_thread.join();
    recv_thread.join();

    MPI_Finalize();
    return 0;
}
