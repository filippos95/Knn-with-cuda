



//
// Created by Filippos Kasioulis on 25/01/2019.
//

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <sys/time.h>


typedef  struct  {
    float x;
    float y;
    float z;

}Point;




__global__ void knn_search(Point* allcPoints,Point* allqPoints,int* perBlockcPoints,int* perBlockqPoints,int* startingPointc,int* startingPointq,Point* knn,float* knn_dist){


    __shared__ Point shrMem[2028];
    int blockId=blockIdx.x+blockIdx.y*gridDim.x+blockIdx.z*gridDim.x*gridDim.y;

    Point qpoint;
    Point cpoint;
    int c,q;
    int number_qpoints=perBlockqPoints[blockId];
    int number_cpoints=perBlockcPoints[blockId];
    int i;

for(i=threadIdx.x;i<number_cpoints;i+=blockDim.x) {
    c = i + startingPointc[blockId];
    shrMem[i] = allcPoints[c];
}
    __syncthreads();

 for(i=threadIdx.x;i<number_qpoints;i+=blocDim.x){

        q = i + startingPointq[blockId];
        qpoint = allqPoints[q];
        float tempDist = 0;
        float minCanDist = 999;
        Point minCanPoint;
        for (int c = 0; c < number_cpoints; c++) {

            tempDist = pow((shrMem[c].x - qpoint.x), 2) + pow((shrMem[c].y - qpoint.y), 2) +
                       pow((shrMem[c].z - qpoint.z), 2);
            tempDist = sqrt(tempDist);


            if (tempDist < minCanDist) {
                minCanDist = tempDist;
                minCanPoint = shrMem[i];

            }
        }

        knn[q] = minCanPoint;
        knn_dist[q] = minCanDist;
    }

    //float minBounds = 999;
    //float tempDistancex = qpoints[q].x-(block_length * blockId.x);
    //if ((tempDistancex > block_length - tempDistancex)&&(blockId.x<gridDim.x-1))
    //  tempDistancex = block_length - tempDistancex;
    //min_from_bounds = tempDistancex;
    //float tempDistancey = qpoints[q].y-(block_length * blockId.y);
    //if ((tempDistancey > block_length - tempDistancey)&&(blockId.y<gridDim.y-1))
    //   tempDistancey = block_length - tempDistancey;
    //if (tempDistancey < min_from_bounds)
    //min_from_bounds = tempDistancey;

    // float tempDistancez = qpoints[q].z-(block_length * blockDim.z);
    //if ((tempDistancez > block_length - tempDistancez)&&(blockId.z<gridDim.z-1))
    //tempDistancez = block_length - tempDistancez;
    // if (tempDistancez < min_from_bounds)
    // min_from_bounds = tempDistancez;


    int neighbors[27];
    int size = 0;
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                if (!(i == 0 & j == 0 & k == 0)) {
                    int tempx = blockIdx.x + i;
                    int tempy = blockIdx.y + j;
                    int tempz = blockIdx.z + k;
                    if (!(tempx < 0 | tempy < 0 | tempz < 0 | tempx >= gridDim.x | tempz >= gridDim.y |
                          tempy >= gridDim.z)) {
                        neighbors[size] = tempx + tempy * gridDim.x + tempz * gridDim.x * gridDim.y;
                        size++;
                    }
                }
            }
        }
    }




    float minNeigDist=999;
    Point minNeig;

    for (int k=0;k<size;k++){

        int neighborId=neighbors[k];
        number_cpoints=perBlockcPoints[neighborId];
        number_qpoints=perBlockqPoints[neighborId];

       for(i=threadIdx.x;i<number_cpoints;i=i+blockDim.x){
            c=i+startingPointc[neighborId];
            shrMem[i]=allcPoints[c];
        }
        __syncthreads();
        for(i=threadIdx.x;i<number_qpoints;i+blockDim.x){
            q=i+startingPointq[neighborId];
            qpoint=allqPoints[q];
            float tempDist;
            for(int k=0;k<number_cpoints;k++){

                tempDist = pow((shrMem[k].x - qpoint.x), 2) + pow((shrMem[k].y - qpoint.y), 2) +
                           pow((shrMem[k].z - qpoint.z), 2);
                tempDist = sqrt(tempDist);

                if (tempDist < minNeigDist) {
                    minNeigDist = tempDist;
                    minNeig =shrMem[k];
                }
            }
            if (minNeigDist < knn_dist[q]) {
                knn_dist[q] = minNeigDist;
                knn[q] = minNeig;
        }
     }
    }
}


void PrintPoints(Point *arr,int number){
    for (int i = 0; i < number; i++) {

        printf("Point %d: x=%f y=%f z=%f \n", i,arr[i].x, arr[i].y, arr[i].z);
    }
}

void PrintKnn( Point *arr, Point *arr2,int number ) {
    for (int i = 0; i < number; i++) {

        printf("Nearest to Point %d: x=%f y=%f z=%f--->  x=%f y=%f z=%f\n", i,arr[i].x, arr[i].y, arr[i].z,arr2[i].x, arr2[i].y, arr2[i].z);

    }
}

void generatePoints(Point* arr ,int number,int seed){
    srand(seed);// randomize seed
    for (int i=0;i<number;i++){
        arr[i].x=((float)rand()/(float)RAND_MAX);
        arr[i].y=((float)rand()/(float)RAND_MAX);
        arr[i].z=((float)rand()/(float)RAND_MAX);

    }

}

float distanceOfPoints(  Point p1, Point p2){
    float dist=0;
    dist=pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2);
    dist=sqrt(dist);
    return dist;
}


void putPointsBlocks(int *pointsToBlock, Point* allPoints,int* perBlockPoints,int numberOfPoints,int gridDim,float blockLength) {
    int blockId=0;
    for (int i=0; i<numberOfPoints;i++){

        blockId=floor(allPoints[i].x/blockLength)+floor(allPoints[i].y/blockLength)*gridDim+floor(allPoints[i].z/blockLength)*gridDim*gridDim;
        pointsToBlock[i]=blockId;
        perBlockPoints[blockId]++;

    }
}



void arrangePointsbyblock(int* pointsToBlock, Point* allPoints, Point* newPoints,int numberOfPoints,int numberOfBlocks,int* perBlockPoints,int* startingPoint){
    int k=0;
    int sum=0;
    for(int j=0;j<numberOfBlocks;j++){
        for (int i=0;i<numberOfPoints;i++) {
            if (pointsToBlock[i] == j) {
                newPoints[k] = allPoints[i];
                k++;

            }
        }
        startingPoint[j]=sum;
        sum=sum+perBlockPoints[j];
    }
}



int main(int argc,char** argv){

    int numberOfcPoints = 22;
    numberOfcPoints = pow(2, numberOfcPoints);
    int dimOfGrid = 6;
    dimOfGrid = pow(2, dimOfGrid);
    int numberOfBlocks = pow(dimOfGrid, 3);
    int numberOfqpoints = numberOfcPoints;
    Point* cpoints = (Point*)malloc(numberOfcPoints * sizeof( Point));
    Point* qpoints = (Point*)malloc(numberOfqpoints * sizeof( Point));

    Point* knn = (Point*)malloc(numberOfqpoints * sizeof( Point));
    Point* arrangecpoints=(Point*)malloc(numberOfcPoints*sizeof( Point));
    Point* arrangeqpoints=(Point*)malloc(numberOfqpoints*sizeof(Point));
    int* pointsctoblock = (int*)malloc(numberOfcPoints * sizeof(int));
    int* pointsqtoblock = (int*)malloc(numberOfqpoints * sizeof(int));
    int* perblockcpoints = (int*)malloc(numberOfBlocks * sizeof(int));
    int* perblockqpoints = (int*)malloc(numberOfBlocks * sizeof(int));
    int* startingpoint_c=(int*)malloc(numberOfBlocks*sizeof(int));
    int* startingpoint_q=(int*)malloc(numberOfBlocks*sizeof(int));
    float* knn_Dist=(float*)malloc(numberOfqpoints*sizeof(float));

    generatePoints(cpoints,numberOfcPoints,1);
    generatePoints(qpoints,numberOfqpoints,2);
    //printf("-------------C  POINTS-----------\n");
    //PrintPoints(cpoints,numberOfcPoints);
    //printf("-------------Q POINTS------------\n");
    //PrintPoints(qpoints,numberOfqpoints);
    struct timeval start_t,end_t;
    double par_time;
    gettimeofday(&start_t,NULL);
    float block_length = ((float) 1) / ((float) dimOfGrid);

    //call function for fragmentation
    putPointsBlocks(pointsctoblock, cpoints, perblockcpoints, numberOfcPoints, dimOfGrid, block_length);
    putPointsBlocks(pointsqtoblock, qpoints, perblockqpoints, numberOfqpoints, dimOfGrid, block_length);
    arrangePointsbyblock(pointsctoblock,cpoints,arrangecpoints,numberOfcPoints,numberOfBlocks,perblockcpoints,startingpoint_c);
    arrangePointsbyblock(pointsqtoblock,qpoints,arrangeqpoints,numberOfqpoints,numberOfBlocks,perblockqpoints,startingpoint_q);
    free(pointsctoblock);
    free(pointsqtoblock);
    cudaError_t err;
    Point* arrangecpointsDev;
    err=cudaMalloc(&arrangecpointsDev,numberOfcPoints*sizeof(Point));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy(arrangecpointsDev,arrangecpoints,numberOfcPoints*sizeof(Point),cudaMemcpyHostToDevice);
    Point* arrangeqpointsDev;
    err=cudaMalloc(&arrangeqpointsDev,numberOfqpoints*sizeof(Point));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy(arrangeqpointsDev,arrangeqpoints,numberOfqpoints*sizeof(Point),cudaMemcpyHostToDevice);

    int* startingpointDev_c;
    err=cudaMalloc(&startingpointDev_c,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy(startingpointDev_c,startingpoint_c,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);
    int* startingpointDev_q;

    err=cudaMalloc(&startingpointDev_q,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy(startingpointDev_q,startingpoint_q,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);

    int* perblockcpointsDev;
    err=cudaMalloc(&perblockcpointsDev,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cudaMemcpy(perblockcpointsDev,perblockcpoints,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);
    int *perblockqpointsDev;
    err=cudaMalloc(&perblockqpointsDev,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cudaMemcpy(perblockqpointsDev,perblockqpoints,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);
    Point* knn_Dev=NULL;
    err=cudaMalloc(&knn_Dev,numberOfqpoints*sizeof(Point));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    float* knnDist_Dev=NULL;
    err=cudaMalloc(&knnDist_Dev,numberOfqpoints*sizeof(float));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }


    knn_search<<<dim3(dimOfGrid,dimOfGrid,dimOfGrid),1024>>>(arrangecpointsDev,arrangeqpointsDev,perblockcpointsDev,perblockqpointsDev,startingpointDev_c,startingpointDev_q,knn_Dev,knnDist_Dev);
    gettimeofday(&end_t,NULL);
    par_time = (double)((end_t.tv_usec - start_t.tv_usec)/1.0e6
                        + end_t.tv_sec - start_t.tv_sec);

    cudaMemcpy(knn_Dist,knnDist_Dev,numberOfqpoints*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(knn,knn_Dev,numberOfqpoints*sizeof(Point),cudaMemcpyDeviceToHost);
    printf("Time to )
    //PrintKnn(arrangeqpoints,knn,numberOfqpoints);


    free(cpoints);
    free(qpoints);
    free(arrangecpoints);
    free(arrangeqpoints);
    free(perblockcpoints);
    free(perblockqpoints);
    free(startingpoint_c);
    free(startingpoint_q);
    free(knn);
    free(knn_Dist);

    cudaDeviceReset();

    return 0;
}






