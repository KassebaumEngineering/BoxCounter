#include <stdio.h>
#include <math.h>

/* #define DEBUG */

int getline( char *line ){
    int i;
    char c;

    i = 0;
    while ( (c = getchar()) != EOF ){
        if( c == '\n' ){
            line[i] = 0;
            return c;
        } else {
            line[i++] = c;
            line[i] = 0;
        }
    }
    return c;
};

main()
{
    float dx, dxy;
    char line[100];
    int i,j;

    printf("#Du = 2* D1( e1.q8.Adata ) - D1( e2.q8.Adata )\n----\n");

    i=1; j=1;
    printf("delay %d:\n",j++);
    while ( getline(line) != EOF ){
        if( line[0] != '#' ){
            sscanf(line, "%f %f", &dx, &dxy);
#ifdef DEBUG
printf("%s\n",line);
printf("dx = %f, dxy = %f\n",dx, dxy);
#endif
            printf("%d: %f    %f    %f\n",i++, dx, dxy, (2*dx - dxy) );
            if(i==9) {
                i=1;
                printf("\ndelay %d:\n",j++);
            }
        }
	}
}