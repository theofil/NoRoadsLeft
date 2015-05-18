#include <cstdio>
#include <cstdlib>
#include "print_tree.h"

//  g++ -o idhep idhep.c

int main (int argc, char **argv) {

  int id=0;
  if (argc<2) {
    puts("Usage: idhep <idhep>\n");
    exit(0);
  }

  id=atoi(argv[1]);
  printf("%d = %s\n",id,RealName(id));

  exit(0);

}
