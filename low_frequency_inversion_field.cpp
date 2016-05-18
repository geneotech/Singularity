#include "prime_evolution.h"
#include "cgp.h"
#include <iostream>

#include "img.h"
#include "modulity.h"

#include "factorizator.h"
#include <fstream>
#include <algorithm>

static int bestNum = 0;

static bool recurrent = false;

void _reportLFIFChromo(struct chromosome* bestChromo) {
	printChromosome(bestChromo, 0);

	std::ostringstream ss;
	ss << "lfif/lfif" << bestNum++;

	auto dot = ss.str() + ".dot";
	auto tex = ss.str() + ".tex";
	auto png = ss.str() + ".png";

	saveChromosomeDot(bestChromo, 0, dot.c_str());
	if (!recurrent)
		saveChromosomeLatex(bestChromo, 0, tex.c_str());

	std::cout << "Active nodes:" << getActiveNodes(bestChromo) << std::endl;
}

#include "randval.h"

struct factored_pair {
	int factor_1;
	int factor_2;
	int next_composite_factor_1;
	int next_composite_factor_2;
};

std::vector<factored_pair> pairs_uniquely_twofactored;
std::vector<int> composites;

double LFIFFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {
	double sumError = 0.0;

	int trialcnt = pairs_uniquely_twofactored.size();
	int firsttrial = trialcnt / 5 * 4;

	static thread_local double in[2];

	for (int i = firsttrial; i< trialcnt; i++) {
		//int x = randval(1, 1920);
		//int y = randval(1, x);

		auto& p = pairs_uniquely_twofactored[i];
		int x = p.factor_1;
		int y = p.factor_2;
		
		{
			in[0] = x;
			in[1] = y;

			executeChromosome(chromo, in);
			sumError += std::abs(p.next_composite_factor_1 - getChromosomeOutput(chromo, 0));
			sumError += std::abs(p.next_composite_factor_2 - getChromosomeOutput(chromo, 1));
		}
	}

	return sumError;
}

#define FACTORED_COUNT 100000

full_factors_info factorizations[FACTORED_COUNT];
using namespace std;



static double _M(const int numInputs, const double *inputs, const double *connectionWeights) {
	auto sm = std::min(inputs[0], inputs[1]);

	if (sm <= 0.0)
		return 0.0;
	
	return modulity(std::max(inputs[0], inputs[1]), sm);
}

static double _gcd_l(const int numInputs, const double *inputs, const double *connectionWeights) {
	return gcd_length(inputs[0], inputs[1]);
}

static double _gcd(const int numInputs, const double *inputs, const double *connectionWeights) {
	return gcd(inputs[0], inputs[1]);
}

void run_low_frequency_inversion_field() {
	for (int i = 0; i < FACTORED_COUNT; ++i) {
		factorizations[i] = get_full_factor_info(i);
	}

	pairs_uniquely_twofactored.reserve(FACTORED_COUNT);
	composites.reserve(FACTORED_COUNT);

	for (int i = 2; i < FACTORED_COUNT; ++i) {
		if (factorizations[i].prime_factors_sequence.size() > 1) {
			composites.push_back(i);
		}
	}

	for (int i = 0; i < composites.size() - 1; ++i) {
		auto& c1 = factorizations[composites[i]];
		auto& c2 = factorizations[composites[i + 1]];

		if (
			c2.number - c1.number == 1 
			&&
			c1.has_unique_twofactorization()
			&&
			c2.has_unique_twofactorization()
			)
		{
			pairs_uniquely_twofactored.push_back({
				c1.smaller_twofactor(),
				c1.bigger_twofactor(),
				c2.smaller_twofactor(),
				c2.bigger_twofactor(),
			});
		}
	}

	{
		std::ofstream gen_file("adjacent_twofactors.txt");

		for (auto& p : pairs_uniquely_twofactored) {
			gen_file << "(" << p.factor_1 << ", " << p.factor_2 << ") -> (" << p.next_composite_factor_1 << ", " << p.next_composite_factor_2 << ")\n";
		}
	}

	struct parameters *params = NULL;
	struct chromosome *chromo = NULL;

	reportChromo = _reportLFIFChromo;

	int numInputs = 2;
	int numNodes = 100;
	int numOutputs = 2;
	int nodeArity = 2;

	int numGens = 1000000000;
	double targetFitness = 0.1;
	int updateFrequency = 5000;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	setNumThreads(params, 6);

	//addNodeFunction(params, "add,sub,mul,div,sin,cos,sqrt,sq,ln,1,pi,e");
	addNodeFunction(params, "add,sub,mul,div,sin,sqrt,ln,1,pi");
	addCustomNodeFunction(params, _M, "M", 2);
	addCustomNodeFunction(params, _gcd, "gcd", 2);
	addCustomNodeFunction(params, _gcd_l , "gcd_l", 2);

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);

	setMutationRate(params, 0.1);
	if (recurrent)
		setRecurrentConnectionProbability(params, 0.1);

	printParameters(params);

	setCustomFitnessFunction(params, LFIFFitness, "LFIF");

	chromo = runCGP(params, NULL, numGens);
	_reportLFIFChromo(chromo);

	freeChromosome(chromo);
	freeParameters(params);

	int abcc;
	std::cin >> abcc;
}