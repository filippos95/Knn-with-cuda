
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

struct Point2{
    int x;
    int y;
    int z;
};



typedef  struct {
    int  num_of_blocks;
    int  *neighbor_blocks;

}neighborBlock;

void PrintPoints(struct Point *arr,int number){
    for (int i = 0; i < number; i++) {

        printf("Point %d: x=%f y=%f z=%f \n", i,arr[i].x, arr[i].y, arr[i].z);
    }
}

void PrintKnn(struct Point *arr,struct Point *arr2,int number ) {
    for (int i = 0; i < number; i++) {

        printf("Nearest to Point %d: x=%f y=%f z=%f--->  x=%f y=%f z=%f\n", i,arr[i].x, arr[i].y, arr[i].z,arr2[i].x, arr2[i].y, arr2[i].z);

    }
}

void generatePoints( struct Point* arr ,int number,int seed){
    srand(seed);// randomize seed
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


void putPointsBlocks(int *pointsToBlock,struct Point* allPoints,int* perBlockPoints,int numberOfPoints,int gridDim,float blockLength,struct Point2* pointsToBlockDim) {
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



void getNeighbourBlocks(int gridDim,struct Point2 blockDim,neighborBlock* neighbor) {
    int *temp = malloc( 27 * sizeof(int));//27 the max number of neighbours
    int size=0;

    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                if (!(i == 0 & j == 0 & k == 0)) {
                    int tempx = blockDim.x+i;
                    int tempy = blockDim.y+j;
                    int tempz = blockDim.z+k;
                    if (!(tempx < 0 | tempy < 0 | tempz < 0 | tempx >= gridDim | tempz >= gridDim | tempy >= gridDim)) {
                        temp[size]=tempx+tempy*gridDim+tempz*gridDim*gridDim;
                        size++;
                    }
                }
            }
        }
    }
    (*neighbor).num_of_blocks=size;
    //printf("%d %d %d \n ",blockDim.x,blockDim.y,blockDim.z);
    (*neighbor).neighbor_blocks = malloc(size*sizeof(int));
    memcpy(neighbor->neighbor_blocks , temp , size*sizeof(int));
    free(temp);
}

struct  Point searchMinDistanceblock(neighborBlock* neig,struct Point* allPoints,int numberOfPoints,int* pointsToblock,struct Point pointtosearch  ){
    float minDist=999;
    float tempDist;
    struct Point minCand;
    //printf("number of blocks %d\n",neig->num_of_blocks);
    for (int j=0;j<neig->num_of_blocks;j++){
        for (int i=0;i<numberOfPoints;i++){
            if (neig->neighbor_blocks[j]==pointsToblock[i]){
                //printf("blocks %d",neig->neighbor_blocks[j]);
                tempDist=distanceOfPoints(allPoints[i],pointtosearch);
                if (tempDist<minDist){
                    minDist=tempDist;
                    minCand.x=allPoints[i].x;
                    minCand.y=allPoints[i].y;
                    minCand.z=allPoints[i].z;
                }
            }

        }
    }
    //printf("mincand from neighbors %f %f %f\n",minCand.x,minCand.y,minCand.z);
    free(neig->neighbor_blocks);
    return minCand;

}

void validation(struct Point* knn,struct Point* qpoints,struct  Point* cpoints,int numberPoints){
    float success=0;
    float tempDist;
    float minDist;
    struct Point minCordCan;
    for (int i=0;i<numberPoints;i++){
        tempDist=0;
        minDist=99;
        for(int j=0;j<numberPoints;j++) {
            tempDist = distanceOfPoints(qpoints[i], cpoints[j]);
            if (tempDist < minDist) {
                minDist = tempDist;
                minCordCan.x = cpoints[j].x;
                minCordCan.y=cpoints[j].y;
                minCordCan.z=cpoints[j].z;

            }
        }

        if ((minCordCan.x==knn[i].x)&(minCordCan.y==knn[i].y)&(minCordCan.z==knn[i].z))
        {
            printf("ok %d %f %f %f \n ",i,minCordCan.x,minCordCan.y,minCordCan.z);
            success++;
        }
        else
            printf("not ok %d %f %f %f \n ",i,minCordCan.x,minCordCan.y,minCordCan.z);
    }

    printf("Success ratee %.1f\n",(success/numberPoints)*100);

}



int main(int argc,char** argv) {


    int numberOfcPoints = 15;
    numberOfcPoints = pow(2, numberOfcPoints);
    int dimOfGrid = 4;
    dimOfGrid = pow(2, dimOfGrid);
    int numberOfBlocks = pow(dimOfGrid, 3);
    int numberOfqpoints = numberOfcPoints;
    struct Point *cpoints = malloc(numberOfcPoints * sizeof(struct Point));
    struct Point *qpoints = malloc(numberOfqpoints * sizeof(struct Point));
    struct Point2 *pointsqtoblockDim = malloc(numberOfqpoints * sizeof(struct Point2));
    struct Point2 *pointsctoblockDim = malloc(numberOfcPoints * sizeof(struct Point2));

    struct Point *knn = malloc(numberOfqpoints * sizeof(struct Point));

    int *pointsctoblock = malloc(numberOfcPoints * sizeof(int));
    int *pointsqtoblock = malloc(numberOfqpoints * sizeof(int));
    int *perblockcpoints = malloc(numberOfBlocks * sizeof(int));
    int *perblockqpoints = malloc(numberOfBlocks * sizeof(int));

    generatePoints(cpoints, numberOfcPoints,1);
    generatePoints(qpoints, numberOfqpoints,2);
    printf("-------------C  POINTS-----------\n");
    PrintPoints(cpoints,numberOfcPoints);
    printf("-------------Q POINTS------------\n");
    PrintPoints(qpoints,numberOfqpoints);

    float block_length = ((float) 1) / ((float) dimOfGrid);

//call function for fragmentation
    putPointsBlocks(pointsctoblock, cpoints, perblockcpoints, numberOfcPoints, dimOfGrid, block_length,
                    pointsctoblockDim);
    putPointsBlocks(pointsqtoblock, qpoints, perblockqpoints, numberOfqpoints, dimOfGrid, block_length,
                    pointsqtoblockDim);
//find the primary candidates of each queries
    int count=0;
    for (int q = 0; q < numberOfqpoints; q++) {
        int blockId;
        struct Point2 blockDim;
        blockId=pointsqtoblock[q];
        blockDim.x=pointsqtoblockDim[q].x;
        blockDim.y=pointsqtoblockDim[q].y;
        blockDim.z=pointsqtoblockDim[q].z;
        //printf(" %d %d %d  \n",blockDim.x,blockDim.y,blockDim.z);

        //printf(" %d %d %d  \n",blockDim.x,blockDim.y,blockDim.z);
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
        for (int i = 0; i < k; i++) {
            tempDistance = distanceOfPoints(qpoints[q], primaryCandidates[i]);
            if (tempDistance < minDistCand) {
                minDistCand = tempDistance;
                minCandidateCord.x = primaryCandidates[i].x;
                minCandidateCord.y=primaryCandidates[i].y;
                minCandidateCord.z=primaryCandidates[i].z;
            }
        }
        printf("Point %d minCandidateCord  %f %f %f\n" ,q,minCandidateCord.x,minCandidateCord.y,minCandidateCord.z);
        float min_from_bounds = 999;
        float tempDistancex = qpoints[q].x-(block_length * blockDim.x);
        if ((tempDistancex > block_length - tempDistancex)&&(blockDim.x<dimOfGrid-1))
            tempDistancex = block_length - tempDistancex;
        min_from_bounds = tempDistancex;
        float tempDistancey = qpoints[q].y-(block_length * blockDim.y);
        if ((tempDistancey > block_length - tempDistancey)&&(blockDim.y<dimOfGrid-1))
            tempDistancey = block_length - tempDistancey;
        if (tempDistancey < min_from_bounds)
            min_from_bounds = tempDistancey;


        float tempDistancez = qpoints[q].z-(block_length * blockDim.z);
        if ((tempDistancez > block_length - tempDistancez)&&(blockDim.z<dimOfGrid-1))
            tempDistancez = block_length - tempDistancez;
        if (tempDistancez < min_from_bounds)
            min_from_bounds = tempDistancez;


        if (min_from_bounds < minDistCand) {
            printf("%d\n",q);
            neighborBlock neighbour;
            getNeighbourBlocks(dimOfGrid, blockDim, &neighbour);
            //printf("number of neighbors %d\n", neighbour.num_of_blocks);
            struct Point minCordNeig;
            minCordNeig = searchMinDistanceblock(&neighbour, cpoints, numberOfcPoints, pointsctoblock, qpoints[q]);


            if (distanceOfPoints(minCordNeig, qpoints[q]) < minDistCand)
                minCandidateCord = minCordNeig;

        }
        knn[q] = minCandidateCord;
        free(primaryCandidates);

    }


    printf("--------------------------------------K-nn---------------------------------\n");
    PrintKnn(qpoints, knn, numberOfqpoints);

    validation(knn,qpoints,cpoints,numberOfqpoints);
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