
#include "glip.h"


int main(){
  glip::sImage32 *I,*AC;
  glip::sAdjRel *A;
  char filename[512];
  int T;

  printf("Enter the filename: ");  
  scanf("%s", filename);
  printf("Volume threshold T: ");
  scanf("%d",&T);
  I = glip::Image32::Read(filename); //(char *)"./dat/cheese.pgm");
  A = glip::AdjRel::Circular(1.0);
  AC = glip::Image32::VolumeClosing(A, I, T);

  glip::Image32::Write(AC, (char *)"VC.pgm");

  glip::Image32::Destroy(&I);
  glip::Image32::Destroy(&AC);
  glip::AdjRel::Destroy(&A);
  return 0;
}

