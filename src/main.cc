// Copyright 2014 Yu Jing<yujing5b5d@gmail.com>
//#include "ThreadPool.h"

#include "blesimrank/blesimrank.hh"


int main(int argc, char* argv[]) {
    Blesimrank br("synthetic",11,10,5,100,10000,20,0.01);
    br.preproc();
    br.save();
    return 0;
}

