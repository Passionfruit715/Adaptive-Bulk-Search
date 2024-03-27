#include <float.h>
#include <stdio.h>

#define true 1
#define false 0

__global__ void localSearch(int **targetBuffer_d, int **solutionBuffer_d, float *solutionBufferValue_d, int **permutationPool_d,
                            int **permutationPoolInverse_d, float ***W_d, int numBits, int numPermutations, int numSegments)
{
    // Energy difference delta_i declared in shared memory (see definition of delta_i in paper)
    // Array size is static, either set an upper bound, or use extern keyword to dynamically specify size (where there should be modifications on calling the kernel function)
    // We use static array size for now, can change later

    // Record the current solution, best solution, and their corresponding values in shared memory, so that all the threads can see them
    __shared__ float delta[512];
    __shared__ int currentSolution[512] = {0}; // Initialize the current solution to be a vector of zeros
    __shared__ int bestSolution[512] = {0};
    __shared__ float min_delta = FLT_MAX;
    __shared__ float currentValue = 0;
    __shared__ float bestValue = FLT_MAX;

    // newly added for straight search and cyclic-min
    __shared__ int k = 0;
    __shared__ int flag = numBits;
    __shared__ int isInitial = true;
    __shared__ int numSegBits = (numBits - 1) / numSegments + 1;
    // Record target solution we want to approach
    // And do bit permutation, the permutation used is blockIdx.x mod numPermutations
    __shared__ int targetSolution[512];
    int pIndex = blockIdx.x % numPermutations;

    // Permutation
    if (threadIdx.x < numBits)
        targetSolution[permutationPool_d[pIndex][threadIdx.x]] = targetBuffer_d[blockIdx.x][threadIdx.x];
    __syncthreads();
    // REMARK: the weight matrix used here is W_d[pIndex] (it is a pointer to pointer of float)
    // One question needs discussion: should we copy the permutated weight matrix on shared memory of this block? The matrix is large

    // Implement straight search (Exactly the same as Algorithm 3 in paper)
    while (flag != 0) // flag is set to be numBits, when every elem of array does not change, we get T from 0.
    {
        flag = numBits; // reset flag
        min_delta = FLT_MAX;

        __syncthreads();

        if (currentSolution[threadIdx.x] != targetSolution[threadIdx.x]) // only calculate when the bit is different from target
        {
            if (isInitial == true) // Initialization from vector 0; using formula in the literature to compute
            {
                delta[threadIdx.x] = W_d[pIndex][threadIdx.x][threadIdx.x];
                for (int j = 0; j < numBits; j++)
                    delta[threadIdx.x] += 2 * W_d[pIndex][threadIdx.x][j] * currentSolution[j];
                delta[threadIdx.x] *= -2 * currentSolution[threadIdx.x] + 1; // \phi(x)
            }
            __syncthreads();
            isInitial = false;

            if (delta[threadIdx.x] < min_delta)
            {
                min_delta = delta[threadIdx.x];
                k = threadIdx.x;
            }
        }

        __syncthreads();

        // update each delta i
        if (thread.Idx == k) // we only need to update once in each block, rather than for all threads
        {
            currentValue += min_delta;
            if (currentValue < bestValue)
            {
                bestValue = currentValue;
                for (int i = 0; i < numBits; i++) // Since in one thread we cannot use parallel computing, we have to use loop
                    bestSolution[i] = currentSolution[i];
            }
            delta[threadIdx.x] *= -1;
        }
        else
        {
            delta[threadIdx.x] += 2 * W_d[pIndex][threadIdx.x][k] * (-2 * currentSolution[threadIdx.x] + 1) * (-2 * currentSolution[k] + 1);
        }

        // flip with minimum delta k
        if (threadIdx.x == 0)
            currentSolution[k] = 1 - currentSolution[k];

        __syncthreads();

        // judge whether to get out of while-loop
        if (currentSolution[threadIdx.x] == targetSolution[threadIdx.x])
            flag--;

        __syncthreads();
    }

    // Now that we know energy value of the target solution,
    // Implement cyclic-Min
    // Question: each segment is dealt with sequentially, where no parallelism exists between segments, wasting other threads
    for (int i = 0; i < numBits; i += numSegBits)
    {
        min_delta = FLT_MAX;
        __syncthreads();

        if (threadIdx.x >= i && threadIdx.x < i + numSegBits)
        {
            delta[threadIdx.x] = W_d[pIndex][threadIdx.x][threadIdx.x];
            for (int j = 0; j < numBits; j++)
                delta[threadIdx.x] += 2 * W_d[pIndex][threadIdx.x][j] * currentSolution[j];
            delta[threadIdx.x] *= -2 * currentSolution[threadIdx.x] + 1;

            if (delta[threadIdx.x] < min_delta)
            {
                min_delta = delta[threadIdx.x];
                k = threadIdx.x;
            }

            __syncthreads();

            if (thread.Idx == k)
                currentValue += min_delta;

            if (currentValue < bestValue)
            {
                bestValue = currentValue;
                if (threadIdx.x < numBits)
                    bestSolution[threadIdx.x] = currentSolution[threadIdx.x];
            }

            __syncthreads();

            if (thread.Idx == k)
                currentSolution[k] = 1 - currentSolution[k];
        }
    }

    // Reverse back the bestSolution, and write it on solution buffer
    if (threadIdx.x < numBits)
    {
        solutionBuffer_d[blockIdx.x][threadIdx.x] = bestSolution[permutationPoolInverse_d[pIndex][threadIdx.x]];
    }
    solutionBufferValue_d[blockIdx.x] = bestValue;
}