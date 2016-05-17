#include <tuple>
#include <iostream>

#include "cgp.h"
#include <algorithm>
using namespace std;
#include "prime_evolution.h"
#include "primes.h"
#include "modulity.h"
#include "tests/tests.h"

double multipliedError(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {
	unsigned long long usemiprime = 7109 * 3301;
	double semiprime = usemiprime;
	
	executeChromosome(chromo, &semiprime);
	
	double p1 = getChromosomeOutput(chromo, 0);
	double p2 = getChromosomeOutput(chromo, 1);
	double fact1 = 7109;
	double fact2 = 3301;
	double d11 = std::abs(p1 - fact1);
	double d22 = std::abs(p2 - fact2);
	double d12 = std::abs(p1 - fact2);
	double d21 = std::abs(p2 - fact1);

	double dd[4] = { d11, d22, d12, d21 };

	return d11 + d22 + getActiveNodes(chromo) / 120.0; *std::min_element(dd, dd + 4);

	unsigned long long ip1 = std::abs(p1);
	unsigned long long ip2 = std::abs(p2);
	if (ip1 == 1) ip1 = 0;
	if (ip2 == 1) ip2 = 0;

	unsigned long long prod = ip1 * ip2;
	
	double fit = prod > usemiprime ? prod - usemiprime : usemiprime - prod;
	
	if (fit == 0.0) {
		cout << "Factors: " << ip1 << " " << ip2 << endl;
		return 0.0;
	}

	return fit * getActiveNodes(chromo);
}

void _reportChromo(struct chromosome* bestChromo) {
	printChromosome(bestChromo, 0);
	saveChromosomeDot(bestChromo, 0, "bestChromo.dot");
	saveChromosomeLatex(bestChromo, 0, "bestChromo.tex");
	std::cout << "( " << multipliedError(nullptr, bestChromo, nullptr) << " ) dblfactors:" << std::fixed;

	double o0 = getChromosomeOutput(bestChromo, 0);
	double o1 = getChromosomeOutput(bestChromo, 1);
	std::cout << o0 << " " << o1 << std::endl;

	{
		unsigned long long usemiprime = 2791859169072907;
		double semiprime = usemiprime;
		executeChromosome(bestChromo, &semiprime);

		double p1 = getChromosomeOutput(bestChromo, 0);
		double p2 = getChromosomeOutput(bestChromo, 1);

		std::cout << "Other case: " << getChromosomeOutput(bestChromo, 0) << " " << getChromosomeOutput(bestChromo, 1) << std::endl;
	}

	{
		unsigned long long usemiprime = 112589788668691;
		double semiprime = usemiprime;
		executeChromosome(bestChromo, &semiprime);

		double p1 = getChromosomeOutput(bestChromo, 0);
		double p2 = getChromosomeOutput(bestChromo, 1);

		std::cout << "Other case: " << getChromosomeOutput(bestChromo, 0) << " " << getChromosomeOutput(bestChromo, 1) << std::endl;
	}

	std::cout << "Active nodes:" << getActiveNodes(bestChromo) << std::endl;
}

int main() {
	is_prime_lookup_init();

	//run_modulity_cgp();
	//run_single_prime_cgp();
	analyticmodulity();
	return 0;

	struct parameters *params = NULL;
	struct chromosome *chromo = NULL;
	
	reportChromo = _reportChromo;

	int numInputs = 1;
	int numNodes = 30;
	int numOutputs = 2;
	int nodeArity = 2;

	int numGens = 100000000;
	double targetFitness = 0.0;
	int updateFrequency = 20000;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
	setCustomFitnessFunction(params, multipliedError, "SEPRI");

	addNodeFunction(params, "add,sub,mul,div,sin,sqrt,sq,cube,tan,pi,ln,e");

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);

	printParameters(params);

	chromo = runCGP(params, nullptr, numGens);

	printChromosome(chromo, 0);

	saveChromosomeDot(chromo, 0, "chromo.dot");
	saveChromosomeLatex(chromo, 0, "chromo.tex");

	freeChromosome(chromo);
	freeParameters(params);

	int abcc;
	cin >> abcc;
}