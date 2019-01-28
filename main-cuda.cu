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


typedef  struct  {
    float x;
    float y;
    float z;

}Point;




__global__ void knn(float* allcPoints,float* allqPoints,int * perBlockcPoints,int* perBlockqPoints,int* startingPointc,int* startingPointq,int* knn,int* knn){

    __shared__ Point shrMem[1024];
    int blockId=blockIdx.x+blockIdx.y*gridDim.x+blockIdx.z*gridDim.x*gridDim.y;

    Point qpoint;
    Point cpoint;
    int c,q;
    int number_qpoints=perBlockqPoints[blockId];
    int number_cpoints=perBlcokcPoints[blockId];
    int i=threadId.x;


    if(i<number_cpoints){
         c=i+startingPoint[blockId];
         shrMem[threadId.x]=allcPoints[c];

    }
    __syncthreads();

    if(i<number_qpoints) {

        q = i + startingPointq[blockId];
        qpoint = allqpoints[q];
        float tempDist = 0;
        float minCanDist=999;
        struct Point minCanPoint;
        for (int c = 0; c < number_cpoints; c++) {

            tempDist = pow((shrMem[c].x - qpoint.x), 2) + pow((shrMem[c].y - qpoint.y), 2) +
                       pow((shrMem[c].z - qpoint.z), 2);
            tempDist = sqrt(tempDist);


            if (tempDist < minCanDis) {
                minCanDis = tempDist;
                minCanPoint = shrMem[c];

            }
        }

        knn[q]=minCanPoint;
        knn_dist[q]=minCanDist;
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
        int size=0;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {
                    if (!(i == 0 & j == 0 & k == 0)) {
                        int tempx = blockIdx.x+i;
                        int tempy = blockIdx.y+j;
                        int tempz = blockIdx.z+k;
                        if (!(tempx < 0 | tempy < 0 | tempz < 0 | tempx >= gridDim | tempz >= gridDim | tempy >= gridDim)) {
                            neighbors[size]=tempx+tempy*gridDim+tempz*gridDim*gridDim;
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

     if(i<number_cpoints]){
         c=i+startingPointc[neighborId];
         shrMem[i]=allcPoints[c];
     }
     __syncthreads();
     if(i<number_qpoints){
      q=i+startingPointq[neighborId];
      qpoint=allqPoints[q];
      for(int c=0;c<number_qpoints;c++){

          tempDist = pow((shrMem[c].x - qpoint.x), 2) + pow((shrMem[c].y - qpoint.y), 2) +
                     pow((shrMem[c].z - qpoint.z), 2);
          tempDist = sqrt(tempDist);

          if (tempDist < minNeigDist) {
              minNeigDis = tempDist;
              minNeig =shrMem[c];
          }
      }
     }
    }
    if(i<number_qpoints) {
        if (minNeigDist < knn[q]) {
            knn[q] = minNeigDist;
            knn[q] = minNeig;
        }
    }
 }


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
        pointsToBlock[i]=blockId;
        perBlockPoints[blockId]++;

    }
}



void arrangePointsbyblock(int *pointsToBlock,struct Point* allPoints,struct Point* newPoints,int numberOfPoints,int numberOfBlocks,int* perBlockPoints,int* startingPoint){
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
    struct Point *cpoints = malloc(numberOfcPoints * sizeof(struct Point));
    struct Point *qpoints = malloc(numberOfqpoints * sizeof(struct Point));

    //struct Point *knn = malloc(numberOfqpoints * sizeof( Point));
    struct Point *arrangecpoints=malloc(numberOfcPoints*sizeof( Point));
    struct Point *arrangeqpoints==malloc(numberOfqPoints*sizeof(Point));
    int *pointsctoblock = malloc(numberOfcPoints * sizeof(int));
    int *pointsqtoblock = malloc(numberOfqpoints * sizeof(int));
    int *perblockcpoints = malloc(numberOfBlocks * sizeof(int));
    int *perblockqpoints = malloc(numberOfBlocks * sizeof(int));
    int *startingpoint_c=malloc(numberOfBlocks*sizeof(int));
    int *startingpoint_q=malloc(numberOfBlocks*sizeof(int));

    generatePoints(cpoints, numberOfcPoints,1);
    generatePoints(qpoints, numberOfqpoints,2);
    //printf("-------------C  POINTS-----------\n");
    //PrintPoints(cpoints,numberOfcPoints);
    //printf("-------------Q POINTS------------\n");
    //PrintPoints(qpoints,numberOfqpoints);
    struct timeval start_t,end_t;
    gettimeofday(&start_t,NULL);
    float block_length = ((float) 1) / ((float) dimOfGrid);

    //call function for fragmentation
    putPointsBlocks(pointsctoblock, cpoints, perblockcpoints, numberOfcPoints, dimOfGrid, block_length);
    putPointsBlocks(pointsqtoblock, qpoints, perblockqpoints, numberOfqpoints, dimOfGrid, block_length,pointsqtoblockDim);
    arrangePointsbyblock(pointsctoblock,cpoints,arrangecpoints,numberOfcPoints,numberOfBlocks,perblockcpoints,startingpoint_c);
    arrangePointsbyblock(pointsqtoblock,qpoints,arrangeqpoints,numberOfqpoints,numberOfBlocks,perblockqpoints,startingpoint_q);
    free(pointsctoblock);
    free(pointsqtoblock);
    cudaError_t err;
    err=cudaMalloc(&arrangecpointsDev,numberOfcPoints*sizeof(Point));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy=(arrangecpointsDev,arrangecpoints,numberOfcPoints*sizeof(Point),cudaMemcpyHostToDevice);

    err=cudaMalloc(&arrangeqpointsDev,numberOfcPoints*sizeof(Point));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy=(arrangeqpointsDev,arrangeqpoints,numberOfqPoints*sizeof(Point),cudaMemcpyHostToDevice);


    err=cudaMalloc(&startingpointDev_c,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy(startingpointDev_c,startingpoint_c,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);

    err=cudaMalloc(&startingpointDev_q,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)                             // `cudaSuccess` is provided by CUDA.
    {
        printf("Error: %s\n", cudaGetErrorString(err)); // `cudaGetErrorString` is provided by CUDA.
    }
    cudaMemcpy(startingpointDev_q,startingpoint_q,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);


    err=cudaMalloc(&perblockcpointsDev,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cudaMemcpy(perblockcpointsDev,perblockcpoints,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);

    err=cudaMalloc(&perblockqpointsDev,numberOfBlocks*sizeof(int));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cudaMemcpy(perblockqpointsDev,perblockqpoints,numberOfBlocks*sizeof(int),cudaMemcpyHostToDevice);

    err=cudaMalloc(&knn_Dev,numberOfqpoints*sizeof(Point));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    gridDim=dim3(dimOfGrid,dimOfGrid,dimOfGrid);

    err=cudaMalloc(&knnDist_Dev,numberOfqpoints*sizeof(float));
    if (err != cudaSuccess)
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    gridDim=dim3(dimOfGrid,dimOfGrid,dimOfGrid);

    knn<<<griDim,1024>>>(arrangecpointsDev,arrangeqpointsDev,perblockcpointsDev,perblockqpointsDev,perblockqpointsDev,startingpointDev_c,startingpointDev_q,knn_Dev,knnDist_dev)
    getimeofday(&end_t,NULL);
    par_time = (double)((end_t.tv_usec - start_t.tv_usec)/1.0e6
                        + endw_t.tv_sec - start_t.tv_sec);

    cudaMemcpy(knn_Dist,knnDist_dev,numberOfqpoints*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(knn,knnDev,numberOfqpoints*sizeof(Point),cudaMemcpyDeviceToHost);

    PrintKnn(qpoints,knn,numberOfqpoints);


    free(cpoints);
    free(qpoints);
    free(arrangecpoints);
    free(arrangeqpoints);
    free(perblockcpoints);
    free(perblockqpoints);
    free(startingpoint_c);
    free(startingpoint_q);
    
    cudaDeviceReset();
    return 0;
}










