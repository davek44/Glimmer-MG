#include <iostream>
#include "base.h"
#include "instant.h"

using namespace::std;

int base_int = DEFAULT_INT;
char base_charstar[4];

int main() {
     base_charstar[0] = 'a';
     base_charstar[1] = 'b';
     base_charstar[2] = 'c';
     base_charstar[3] = '\0';
     
     cout << base_int << endl;
     cout << base_float << endl;
     cout << base_charstar << endl;

     f1();
     cout << base_int << endl;
     cout << base_charstar << endl;
}

static void f1() {
     base_int++;
     base_charstar[0] = 'x';
}
