#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define POP_SIZE 8         // number of chromosome in population
#define CHROM_LENGTH 8     // size of chromosome (number of genes)
#define MAX_ITER 300       // maximum generation (iteration)
#define MUTATION_RATE 0.05 // rate of mutation
#define CROSS_PROB 0.75    // probability of crossover
#define THRESHOLD 0.75     // stopping propoortion
#define MAX_OCCUR 272      // max key for each chromosome

const int W[8][8] = {
    {4, -3, 2, 7, 5, -1, 0, 2},
    {-3, 8, -5, 1, 0, -2, 3, -4},
    {2, -5, 7, -6, 1, 3, -4, 5},
    {7, 1, -6, 9, 3, -2, 2, -3},
    {5, 0, 1, 3, -10, -4, 3, -2},
    {-1, -2, 3, -2, -4, 8, 1, -1},
    {0, 3, -4, 2, 3, 1, 7, -6},
    {2, -4, 5, -3, -2, -1, -6, 9}
};

int occur[MAX_OCCUR];

typedef struct
{
    // using type of int to store genes
    int genes[CHROM_LENGTH];
    double fitness;
} Chromosome;

Chromosome population[POP_SIZE];
Chromosome new_population[POP_SIZE];

void initialize_population()
{
    for (int i = 0; i < POP_SIZE; i++)
    {
        for (int j = 0; j < CHROM_LENGTH; j++)
            population[i].genes[j] = rand() % 2; // 0 or 1
        population[i].fitness = 0.0;
    }
}

// x^TWx
double evaluate(Chromosome chrom)
{
    double value = 0.0;
    for (int i = 0; i < CHROM_LENGTH; i++)
        for (int j = 0; j < CHROM_LENGTH; j++)
            value += chrom.genes[i] * chrom.genes[j] * W[i][j];
    return value;
}

int *count()
{
    memset(occur, 0, sizeof(occur));

    int *flag = (int *)malloc(2 * sizeof(int));

    for (int i = 0; i < POP_SIZE; i++)
    {
        int key = 0;
        for (int j = 0; j < CHROM_LENGTH; j++)
        {
            key *= 2;
            key += population[i].genes[j];
        }

        occur[key] += 1;

        if ((double)occur[key] / POP_SIZE >= THRESHOLD)
        {
            flag[0] = 1;
            flag[1] = i;
            return flag;
        }
    }

    int max_prop = 0;

    for (int i = 0; i < POP_SIZE; i++)
        if (occur[i] > max_prop)
            flag[1] = i;
    flag[0] = 0;

    return flag;
}

// Proportionate roulette wheel selection
void selection()
{
    double total_fitness = 0.0;
    double cum_prob[POP_SIZE] = {0.0};
    double max_fitness = 0.0;

    for (int i = 0; i < POP_SIZE; i++)
    {
        population[i].fitness = evaluate(population[i]);

        // fine the maximum fitness
        if (population[i].fitness > max_fitness)
            max_fitness = population[i].fitness;
    }

    for (int i = 0; i < POP_SIZE; i++)
    {

        // Since we want to minimize our energy, solution with low energy should have high fitness
        // We inverse the fitness by subtracting current fitness from the maximum one, which also ensures all fitness greater or equal to 0
        // +1 to ensure strctly greater than 0 (those with low fitness also have chance to be chosen)
        population[i].fitness = max_fitness - population[i].fitness + 1;
        total_fitness += population[i].fitness;
    }

    for (int i = 0; i < POP_SIZE; i++)
    {
        // cumulative probability
        cum_prob[i] = (i == 0) ? (population[i].fitness / total_fitness) : (cum_prob[i - 1] + population[i].fitness / total_fitness);
    }

    for (int i = 0; i < POP_SIZE; i++)
    {
        double rand_num = (double)rand() / RAND_MAX;
        for (int j = 0; j < POP_SIZE; j++)
        {

            // low fitness with low probability to be chosen
            if (rand_num <= cum_prob[j])
            {
                new_population[i] = population[j];
                break;
            }
        }
    }

    for (int i = 0; i < POP_SIZE; i++)
        population[i] = new_population[i];
}

void crossover()
{
    for (int i = 0; i < POP_SIZE; i += 2)
    {
        double rand_num = (double)rand() / RAND_MAX;
        if (rand_num < CROSS_PROB)
        {
            int crossover_point = rand() % CHROM_LENGTH;

            // even chromosomes with odd ones to crossover
            for (int j = crossover_point; j < CHROM_LENGTH; j++)
            {
                int temp = population[i].genes[j];
                population[i].genes[j] = population[i + 1].genes[j];
                population[i + 1].genes[j] = temp;
            }
        }
    }
}

void mutation()
{
    for (int i = 0; i < POP_SIZE; i++)
    {
        for (int j = 0; j < CHROM_LENGTH; j++)
        {
            double rand_num = (double)rand() / RAND_MAX;
            if (rand_num < MUTATION_RATE)
                population[i].genes[j] = 1 - population[i].genes[j]; // 1 to 0 or 0 to 1
        }
    }
}

int main()
{
    int *flag;
    printf("Hello World");

    srand(time(NULL));

    initialize_population();

    for (int i = 0; i < MAX_ITER; i++)
    {
        flag = count();

        if (flag[0] == 1)
            break;

        selection();
        crossover();
        mutation();
    }

    for (int i = 0; i < CHROM_LENGTH; i++)
        printf("%d ", population[flag[1]].genes[i]);
    printf("\n");

    free(flag);

    return 0;
}