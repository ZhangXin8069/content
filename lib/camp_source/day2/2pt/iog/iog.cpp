//#include <lat-io.h>
//#include "iog_class.h"
#include "io_general_class.h"
#include "io_general.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

extern "C"
{
    int getsize(char *iog_name)
    {
        int size = 0;

        char *bbf = new char[1000];

        sprintf(bbf,"%s" ,iog_name);
        fstream file_exit;
        file_exit.open(bbf, ios::in);
        if (!file_exit)
        {
            printf("%s not found!\n", iog_name);
            abort();
        }
        general_data_base iog(bbf);

        iog.load();
        size = iog.get_size();

        return size;
    }

    iog_type *getdat(char *iog_name)
    {
        char *bbf = new char[1000];

        sprintf(bbf,"%s",iog_name);
        fstream file_exit;
        file_exit.open(bbf, ios::in);
        if (!file_exit)
        {
            printf("%s not found!\n", iog_name);
            abort();
        }
        general_data_base iog(bbf);

        iog.load();

        iog_type *ptr = (iog_type *)malloc(sizeof(iog_type));
        ptr = iog.get_dat();

        //for(int i = 0; i<100; i=i+4)
        // {
        // printf("%d --- %13.5e   %13.5e   %13.5e   %13.5e ---\n", i, ptr->Re[i], ptr->Re[i+1],ptr->Re[i+2],ptr->Re[i+3] );
        // }

        return ptr;
    }
}