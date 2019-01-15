//
// Created by Filippos Kasioulis on 13/01/2019.
//

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <time.h>


struct Point {
    float x;
    float y;
    float z;

};


void PrintPoints(struct Point *arr,int number ) {
    for (int i = 0; i < number; i++) {
        printf("Point %d: x=%f y=%f z=%f\n", i,arr[i].x, arr[i].y, arr[i].z);
    }
}

void generatePoints( struct Point *arr ,int number){
    srand(time(NULL));// randomize seed
    for (int i=0;i<number;i++){
        arr[i].x=((float)rand()/(float)RAND_MAX);
        arr[i].y=((float)rand()/(float)RAND_MAX);
        arr[i].z=((float)rand()/(float)RAND_MAX);

    }

}


int main(int argc,char** argv){

int numberOfPoints=atoi(argv[1]);
numberOfPoints=pow(2,numberOfPoints);
int dimOfGrid=atoi(argv[2]);
dimOfGrid=pow(2,dimOfGrid);




struct Point *points=malloc(numberOfPoints*sizeof(struct Point));
generatePoints(points,numberOfPoints);
PrintPoints(points,numberOfPoints);

return 0;


}