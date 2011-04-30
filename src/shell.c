#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
int main()
{
  char *text = calloc(1024,sizeof(char));
  text = getusershell();
  printf("%s\n",text);
  return 0;
}
    
