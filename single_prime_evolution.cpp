#include "prime_evolution.h"
#include "cgp.h"
#include <iostream>

void _reportSinglePrimeChromo(struct chromosome* bestChromo) {
	printChromosome(bestChromo, 0);
	saveChromosomeDot(bestChromo, 0, "bestChromo.dot");
	saveChromosomeLatex(bestChromo, 0, "bestChromo.tex");
	std::cout << "Active nodes:" << getActiveNodes(bestChromo) << std::endl;
}

double singlePrimeFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {
	unsigned long long unprime = 81;
	double nprime = unprime;

	executeChromosome(chromo, &nprime);

	double p1 = getChromosomeOutput(chromo, 0);
	double d11 = std::abs(p1 - 419);
	if(d11 < 1.0)
		return d11 + getActiveNodes(chromo);
	else
		return d11 + getActiveNodes(chromo)/120.0;
}

void run_single_prime_cgp() {
	struct parameters *params = NULL;
	struct dataSet *trainingData = NULL;
	struct chromosome *chromo = NULL;

	reportChromo = _reportSinglePrimeChromo;

	int numInputs = 1;
	int numNodes = 25;
	int numOutputs = 1;
	int nodeArity = 2;

	int numGens = 1000000000;
	double targetFitness = 0.0;
	int updateFrequency = 5000;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	addNodeFunction(params, "add,sub,mul,div,sin,cos,sqrt,sq,cube,tan,pi,ln,e");

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);

	printParameters(params);

	setCustomFitnessFunction(params, singlePrimeFitness, "MSE");

	chromo = runCGP(params, trainingData, numGens);
	_reportSinglePrimeChromo(chromo);

	freeDataSet(trainingData);
	freeChromosome(chromo);
	freeParameters(params);

	int abcc;
	std::cin >> abcc;
}