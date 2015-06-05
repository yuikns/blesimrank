// Copyright 2014 Yu Jing<yujing5b5d@gmail.com>
//#include "ThreadPool.h"

#include <cstdlib>

#include "blesimrank/blesimrank.hh"
#include "argcv/timer/timer.hh"

using argcv::timer::timer;

int main(int argc, char* argv[]) {
    if (argc < 9) {
        fprintf(stderr, "usage:\n\t%s dataset T P Q R1 R2 K theta\n",argv[0]);
    } else {
        timer t;
        t.label("start");
        Blesimrank br(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]),
                      atoi(argv[7]), atof(argv[8]));
        t.label("init");
        printf("data loaded , time cost : %f ms \n",t.between("start","init"));fflush(NULL);
        br.preproc();
        t.label("pp");
        printf("data pre-proc , time cost : %f ms \n",t.between("pp","start"));fflush(NULL);
        br.save();
        t.label("done");
        printf("query and save all , time cost : %f ms \n",t.between("done","pp"));fflush(NULL);
    }
    return 0;
}
