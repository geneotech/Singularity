#include "prime_evolution.h"
#include "cgp.h"
#include <iostream>

void _reportPrimeChromo(struct chromosome* bestChromo) {
	printChromosome(bestChromo, 0);
	saveChromosomeDot(bestChromo, 0, "bestChromo.dot");
	saveChromosomeLatex(bestChromo, 0, "bestChromo.tex");
	std::cout << "Active nodes:" << getActiveNodes(bestChromo) << std::endl;
}

double primeFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {

	int i, j;
	double squareError = 0;

	if (getNumChromosomeInputs(chromo) != getNumDataSetInputs(data)) {
		printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
		printf("Terminating.\n");
		exit(0);
	}

	if (getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)) {
		printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
		printf("Terminating.\n");
		exit(0);
	}

	for (i = 0; i<getNumDataSetSamples(data); i++) {

		executeChromosome(chromo, getDataSetSampleInputs(data, i));

		for (j = 0; j<getNumChromosomeOutputs(chromo); j++) {

			squareError += pow(getDataSetSampleOutput(data, i, j) - getChromosomeOutput(chromo, j), 2);
		}
	}

	return squareError / (getNumDataSetSamples(data) * getNumDataSetOutputs(data)) + getActiveNodes(chromo)/10000.0;
}

void run_prime_cgp() {
	struct parameters *params = NULL;
	struct dataSet *trainingData = NULL;
	struct chromosome *chromo = NULL;

	reportChromo = _reportPrimeChromo;

	int numInputs = 1;
	int numNodes = 25;
	int numOutputs = 1;
	int nodeArity = 2;

	int numGens = 1000000000;
	double targetFitness = 0.1;
	int updateFrequency = 5000;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	addNodeFunction(params, "add,sub,mul,div,sin,cos,sqrt,sq,cube,tan,pi,ln,e");

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);

	printParameters(params);

	trainingData = initialiseDataSetFromFile("F:/Projects/Singularity/cgp/dataSets/primeSequence.data");
	setCustomFitnessFunction(params, primeFitness, "MSE");

	chromo = runCGP(params, trainingData, numGens);
	_reportPrimeChromo(chromo);

	freeDataSet(trainingData);
	freeChromosome(chromo);
	freeParameters(params);

	int abcc;
	std::cin >> abcc;
}