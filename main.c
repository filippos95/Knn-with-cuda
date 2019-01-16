//
// Created by Filippos Kasioulis on 13/01/2019.
//

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


 struct Point {
    float x;
    float y;
    float z;

};


typedef struct{
    float max_dist;
    float min_dist;
    struct Point k;


}kNeighbour;

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

float distanceOfPoints( struct Point p1,struct Point p2){
    float dist=0;
    dist=pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2);
    dist=sqrt(dist);
    return dist;


}


void putPointsBlocks(int *pointsToBlock,struct Point *allPoints,int *perBlockPoints,int numberOfPoints,int gridDim,int blockLength) {
         int blockId=0;
         for (int i=0; i<numberOfPoints;i++){
           blockId=floor(allPoints[i].x/blockLength)+floor(allPoints[i].y/blockLength)*gridDim+floor(allPoints[i].z/blockLength)*gridDim*gridDim);
           pointsToBlock[i]=blockId;
           perblockpoints[blockId]++;

         }
}

int main(int argc,char** argv){

int numberOfcPoints=atoi(argv[1]);
numberOfcPoints=pow(2,numberOfcPoints);
int dimOfGrid=atoi(argv[2]);
dimOfGrid=pow(2,dimOfGrid);
int numberOfBlocks=pow(dimOfGrid,3);
int numberOfqpoints=numberOfcPoints;
struct Point *cpoints=malloc(numberOfnPoints*sizeof(struct Point));
struct Point *qpoints=malloc(numberOfqpoints*sizeof(struct Point));
int *pointsctoblock=malloc(numberOfcPointsPoints*sizeof(int));
int *pointsqtoblock=malloc(numberOfqpoints*sizeof(int));
int *perblockcpoints=malloc(numberOfBlocks*sizeof(int));
int *perblockqpoints=malloc(numberOfBlocks*sizeof(int));

generatePoints(cpoints,numberOfcPoints);
generatePoints(qpoints,numberOfqpoints);
//PrintPoints(npoints,numberOfcPoints);
//PrintPoints(qpoints,numberOfqpoints);
float block_length= ((float)1)/((float)dimOfGrid);

putPointsBlocks(pointsctoblock,cpoints,perblockcpoints,numberOfcPoints,dimOfGrid,block_length);
putPointsBlocks(pointsqtoblock,qpoints,perblockqpoints,numberOfqpoints,dimOfGrid,block_length);




return 0;

}