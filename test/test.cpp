//matrix_multiply.cpp

#include<mpi.h>  
#include<stdio.h>  
#include <iostream>  
#include<math.h>  
#pragma comment(lib,"mpi.lib")  
#define n 1000  
using namespace std;  
int main(int argv, char *argc[])  
{  
    int rank, p, a;  
    MPI_Init(&argv, &argc);  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    MPI_Comm_size(MPI_COMM_WORLD, &p);  
    MPI_Status status;  
    if (p!=1)  
    a = n / (p - 1);  
    if (rank == 0)  
    {  
        int* A = new int[n*n];  
        int* B = new int[n*n];  
        int* C = new int[n*n];  
    //  int * recptr = NULL;  
        for (int i = 0; i < n; i++)  
        for (int j = 0; j < n; j++)// 时间是 O nn  
        {  
            A[i*n + j] = i + j; //A[i][j]  
            B[i*n + j] = 1; //B[i][j]  
        }  
        if (p == 1)  
        {  
            double tb, te;  
            tb = MPI_Wtime();  
            for (int i = 0; i < n; i++)  
            for (int j = 0; j < n; j++)  
            {  
                C[i*n + j] = 0; //C[i][j]  
                for (int k = 0; k < n; k++)  
                {  
                    C[i*n + j] = A[i*n + k] * B[k*n + j];  
                }  


            }  
            te = MPI_Wtime();  
            cout << "time is " << te - tb;// << "s" << endl;  
        }  


        if (p != 1)  
        {  
            double tb, te;  

            tb = MPI_Wtime();  
            for (int i = 0; i < p-1; i++){//给每个寄存器发送  数组 A，B，C  
                MPI_Send(&A[0+0], n*n, MPI_INT, i+1, 1, MPI_COMM_WORLD);//每个发送 a行，a*n大小的数据   
                MPI_Send(&B[0+0], n*n, MPI_INT, i+1,2, MPI_COMM_WORLD);  

            }  
            for (int i =0; i < p-1; i++)  
                    MPI_Recv(&C[i*a+0], a*n, MPI_INT, i+1,3, MPI_COMM_WORLD, &status);//每个接受 a行，a*n大小的数据   

            te = MPI_Wtime();  
            cout << "time is " << te - tb;// << "s" << endl;  
        }  

        delete[] A;  
        delete[] B;  
        delete[] C;  
    }  

    if (p != 1)  
    if (rank != 0){  
        int* A = new int[n*n];  
        int* B = new int[n*n];  
        int* C = new int[n*n];  

        MPI_Recv(&A[0+0], n*n, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);//从A[0][0]和B[0][0]开始接受  
        MPI_Recv(&B[0+0], n*n, MPI_INT,0, 2, MPI_COMM_WORLD, &status);  

        for (int i =a*(rank-1); i < (a*(rank)); i++)//按照行间隔分，每个cpu计算自己的a行  
        for (int j = 0; j < n; j++)  
        {  
            C[i*n + j] = 0; //C[i][j]  
            for (int k = 0; k < n; k++)  
            {  
                C[i*n + j] = A[i*n + k] * B[k*n + j];  
            }  
        }  
        {//向rank=0发送自己的那a行C,大小是a*n  
            //int * sendptr = &(C[a*(rank - 1)+0]);  
            MPI_Send(&C[a*(rank - 1) + 0], a*n, MPI_INT, 0,3, MPI_COMM_WORLD);//起始地址是C[rank-1][0],大小是a*n  
        }  
    }  
    MPI_Finalize();  
    return 0;  
} 

