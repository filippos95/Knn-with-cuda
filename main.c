//
// Created by Filippos Kasioulis on 13/01/2019.
//

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>



struct Point {
    float x;
    float y;
    float z;

};



typedef  struct {
    int  num_of_blocks;
    int  *neighbor_blocks;

}neighborBlock;

void PrintPoints(struct Point *arr,struct Point *arr2,int number ) {
    for (int i = 0; i < number; i++) {

        printf("Nearest to Point %d: x=%f y=%f z=%f--->  x=%f y=%f z=%f\n", i,arr[i].x, arr[i].y, arr[i].z,arr2[i].x, arr2[i].y, arr2[i].z);

    }
}

void generatePoints( struct Point* arr ,int number){
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


void putPointsBlocks(int *pointsToBlock,struct Point* allPoints,int* perBlockPoints,int numberOfPoints,int gridDim,float blockLength,struct Point* pointsToBlockDim) {
    int blockId=0;
    for (int i=0; i<numberOfPoints;i++){

        blockId=floor(allPoints[i].x/blockLength)+floor(allPoints[i].y/blockLength)*gridDim+floor(allPoints[i].z/blockLength)*gridDim*gridDim;
        pointsToBlockDim[i].x=floor(allPoints[i].x/blockLength);
        pointsToBlockDim[i].y=floor(allPoints[i].y/blockLength);
        pointsToBlockDim[i].z=floor(allPoints[i].z/blockLength);
        pointsToBlock[i]=blockId;
        perBlockPoints[blockId]++;

    }
}



void getNeighbourBlocks(int gridDim,struct Point blockDim,neighborBlock* neighbor) {
    int *temp = malloc( 27 * sizeof(int));//27 the max number of neighbours
    int size=0;

    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                if (!(i == 0 & j == 0 & k == 0)) {
                    int tempx = blockDim.x;
                    int tempy = blockDim.y;
                    int tempz = blockDim.z;
                    if (tempx < 0 | tempy < 0 | tempz < 0 | tempx >= gridDim | tempz >= gridDim | tempy >= gridDim) {
                        temp[size]=tempx+tempy*gridDim+tempy*gridDim*gridDim;
                        size++;
                    }
                }
            }
        }
    }
    (*neighbor).num_of_blocks=size;
    (*neighbor).neighbor_blocks = malloc(size*sizeof(int));
    memcpy(neighbor->neighbor_blocks , temp , size*sizeof(int));
    free(temp);
}

struct  Point searchMinDistanceblock(neighborBlock* neig,struct Point* allPoints,int numberOfPoints,int* pointsToblock,struct Point pointtosearch  ){
    float minDist=999;
    float tempDist;
    struct Point minCand;
    for (int i=0;i<numberOfPoints;i++){
        for (int j=0;j<neig->num_of_blocks;j++){
            if (neig->neighbor_blocks[j]==pointsToblock[i]){
                tempDist=distanceOfPoints(allPoints[i],pointtosearch);
                if (tempDist<minDist){
                    minDist=tempDist;
                    minCand=allPoints[i];
                }
            }

        }
    }

    return minCand;

}




int main(int argc,char** argv) {


    int numberOfcPoints = 6;
    numberOfcPoints = pow(2, numberOfcPoints);
    int dimOfGrid = 1;
    dimOfGrid = pow(2, dimOfGrid);
    int numberOfBlocks = pow(dimOfGrid, 3);
    int numberOfqpoints = numberOfcPoints;
    struct Point *cpoints = malloc(numberOfcPoints * sizeof(struct Point));
    struct Point *qpoints = malloc(numberOfqpoints * sizeof(struct Point));
    struct Point *pointsqtoblockDim = malloc(numberOfqpoints * sizeof(struct Point));
    struct Point *pointsctoblockDim = malloc(numberOfcPoints * sizeof(struct Point));

    struct Point *knn = malloc(numberOfqpoints * sizeof(struct Point));

    int *pointsctoblock = malloc(numberOfcPoints * sizeof(int));
    int *pointsqtoblock = malloc(numberOfqpoints * sizeof(int));
    int *perblockcpoints = malloc(numberOfBlocks * sizeof(int));
    int *perblockqpoints = malloc(numberOfBlocks * sizeof(int));

    generatePoints(cpoints, numberOfcPoints);
    generatePoints(qpoints, numberOfqpoints);

    float block_length = ((float) 1) / ((float) dimOfGrid);

//call function for fragmentation
    putPointsBlocks(pointsctoblock, cpoints, perblockcpoints, numberOfcPoints, dimOfGrid, block_length,
                    pointsqtoblockDim);
    putPointsBlocks(pointsqtoblock, qpoints, perblockqpoints, numberOfqpoints, dimOfGrid, block_length,
                    pointsctoblockDim);
//find the primary candidates of each queries
    for (int q = 0; q < numberOfqpoints; q++) {
        int blockId;
        struct Point blockDim;
        pointsqtoblock[q] = blockId;
        pointsqtoblockDim[q] = blockDim;
        struct Point *primaryCandidates = malloc(perblockcpoints[blockId] * sizeof(struct Point));
        int k = 0;

        for (int c = 0; c < numberOfcPoints; c++) {
            if (pointsctoblock[c] == blockId) {
                primaryCandidates[k] = cpoints[c];
                k++;
            }
        }
        float tempDistance;
        struct Point minCandidateCord;
        float minDistCand = 999;
        for (int i = 0; i < perblockcpoints[blockId]; i++) {
            tempDistance = distanceOfPoints(qpoints[q], primaryCandidates[i]);
            if (tempDistance < minDistCand) {
                minDistCand = tempDistance;
                minCandidateCord = primaryCandidates[i];
            }
        }
        float min_from_bounds = 999;
        float tempDistancex = fmod(qpoints[q].x, (block_length * blockDim.x));
        if (tempDistancex > block_length - tempDistancex)
            tempDistance = block_length - tempDistancex;
        min_from_bounds = tempDistancex;
        float tempDistancey = fmod(qpoints[q].y, (block_length * blockDim.y));
        if (tempDistancey > block_length - tempDistancey)
            tempDistancey = block_length - tempDistancey;
        if (tempDistancey < min_from_bounds)
            min_from_bounds = tempDistancey;

        float tempDistancez = fmod(qpoints[q].z, (block_length * blockDim.z));
        if (tempDistancez > block_length - tempDistancez)
            tempDistancez = block_length - tempDistancez;
        if (tempDistancez < min_from_bounds)
            min_from_bounds = tempDistancez;

        if (min_from_bounds < minDistCand) {
            neighborBlock neighbour;
            getNeighbourBlocks(dimOfGrid, blockDim, &neighbour);
            float minDistNeig;
            struct Point minCordNeig;
            minCordNeig = searchMinDistanceblock(&neighbour, cpoints, numberOfcPoints, pointsqtoblock, qpoints[q]);


            if (distanceOfPoints(minCordNeig, qpoints[q]) < minDistCand)
                minCandidateCord = minCordNeig;
        }
        knn[q] = minCandidateCord;


    }


    printf("------------K-nn----------\n");
    PrintPoints(qpoints, knn, numberOfqpoints);
    free(cpoints);
    free(qpoints);
    free(perblockcpoints);
    free(perblockqpoints);
    free(pointsqtoblockDim);
    free(pointsctoblockDim);
    free(knn);
    free(pointsqtoblock);
    free(pointsctoblock);
    return 0;

}