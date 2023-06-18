#include<iostream>
#include<string.h>
using namespace std;
int main(){
    string user_name;
    cout << "please enter your name: ";
    cin >> user_name;
    cout << '\n'
         << "hello ,"
         << user_name
         << '...and  goodbye!\n';
    return 0;
}
