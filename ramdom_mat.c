#include <stdio.h>
#include <stdlib.h>

int main() {
  int i, j, dim;

  scanf("%d", &dim);
  printf("%d\n", dim);
  
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      if(i==j)
	printf("%d ", 0);
      else
	printf("%d ", rand()%10);
    }
    printf("\n");
  }

  return 0;
}
