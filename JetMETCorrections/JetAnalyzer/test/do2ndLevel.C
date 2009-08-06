#include "TreeAnalyzer.h"
#include <iostream>


void doOneLoop(int i);

void do2ndLevel(int i_begin=1, int i_end=50) {

  
  for(int iLoop=i_begin; iLoop<i_end; ++iLoop) {

    std::cout << "Starting Loop n. " << iLoop << std::endl;
    doOneLoop(iLoop);

  }

} 


void doOneLoop(int i) {

  TreeAnalyzer* t = new TreeAnalyzer(i);
  t->Loop();
  delete t;

}
