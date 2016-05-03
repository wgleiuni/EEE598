/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description: 2016 Spring ASU EEE598 -- Monte Carlo Simulation of GaN MESFET 
 *
 *        Created:  04/27/2016 08:06:06 PM
 *       Compiler:  icc/icpc
 *
 *         Author:  Guanglei Wang (glwang), wgleiuni@gmail.com
 *   Organization:  Arizona State University
 *
 * =====================================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include "DefClass.h"
#include "Global.h"

extern double rd(){
    return ((double)rand()/RAND_MAX);
}

int main (int argc, char * argv[]){
    srand(time(NULL));
    if (argc==1) numdt=1000;
    else {
        numdt=atoi(argv[1]);
    }
    Device one;
    one.go();
    return 0;
}
