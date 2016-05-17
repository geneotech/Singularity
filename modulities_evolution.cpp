#include "prime_evolution.h"
#include "cgp.h"
#include <iostream>

#include "img.h"
#include "modulity.h"

int bestNum = 0;

bool recurrent = false;

void _reportModulityChromo(struct chromosome* bestChromo) {
	printChromosome(bestChromo, 0);

	int offset = 1;
	int width = 1920;
	int height = 1080;
	int xscale = 1;

	setwh(width, height);
	double in[2] = { 0, 0 };

	for (int j = offset, n = 0; j <= width + offset; ++j, n = (j)) {
		for (int i = 2; i <= std::min(n, height); ++i) {
			in[0] = n;
			in[1] = i;

			executeChromosome(bestChromo, in);

			auto val = getChromosomeOutput(bestChromo, 0) + 1;
			
			//std::abs((gcd_length(n, 1) - 1) - modulity(n, i) + 1);

			if (val > 0)
				setpix(j - offset, i, 255 / val);
		}
	}

	std::ostringstream ss;
	ss << "chromosomepics/symphony_of_moduli" << bestNum++;


	auto dot = ss.str() + ".dot";
	auto tex = ss.str() + ".tex";
	auto png = ss.str() + ".png";

	saveChromosomeDot(bestChromo, 0, dot.c_str());
	if (!recurrent)
		saveChromosomeLatex(bestChromo, 0, tex.c_str());

	lodepng::encode(png.c_str(), img, w, h);

	std::cout << "Active nodes:" << getActiveNodes(bestChromo) << std::endl;
}

#include "randval.h"

int bigoff_x = 9000;
int bigoff_y = 4500;
int smoff_x = 29000;
int smoff_y = 1500;

int modulities[2000][2000];
int bigmodulities[2000][2000];

std::pair<int, int> trials[80000];

double modulityFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {
	double sumError = 0.0;

	int trialcnt = 79000;

	static thread_local double in[2];

	for (int i = 0; i< trialcnt; i++) {
		//int x = randval(1, 1920);
		//int y = randval(1, x);
		int x = trials[i].first;
		int y = trials[i].second;

		//{
		//	in[0] = x + smoff_x;
		//	in[1] = y + smoff_y;
		//
		//	executeChromosome(chromo, in);
		//	sumError += std::abs(modulities[x][y] - getChromosomeOutput(chromo, 0));
		//}
		{
			in[0] = x + bigoff_x;
			in[1] = y + bigoff_y;

			executeChromosome(chromo, in);
			sumError += std::abs(bigmodulities[x][y] - getChromosomeOutput(chromo, 0));
		}
	}

	return sumError;
}

void run_modulity_cgp() {
	for (int i = 1; i < 2000; ++i)
		for (int j = 1; j < 2000; ++j)
			modulities[i][j] = modulity(i + smoff_x, j + smoff_y);

	for (int i = 1; i < 2000; ++i)
		for (int j = 1; j < 2000; ++j)
			bigmodulities[i][j] = modulity(i+ bigoff_x, j + bigoff_y);

	//for (auto& t : trials) {
	//	t.first = randval(1, 1920);
	//	t.second = randval(1, t.first);
	//}

	int ti = 0;
	for (int i = 1; i < 400; ++i)
		for (int j = 1; j <= i; ++j) {
			trials[ti].first = i;
			trials[ti].second = j;
			++ti;
		}

	std::cout << "Trials: " << ti << std::endl;

	struct parameters *params = NULL;
	struct chromosome *chromo = NULL;

	reportChromo = _reportModulityChromo;

	int numInputs = 2;
	int numNodes = 200;
	int numOutputs = 1;
	int nodeArity = 2;

	int numGens = 1000000000;
	double targetFitness = 0.1;
	int updateFrequency = 1000;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	setNumThreads(params, 6);

	addNodeFunction(params, "add,sub,mul,div,sin,cos,sqrt,sq,cube,tan,ln,1,pi,e");

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);

	setMutationRate(params, 0.5);
	if(recurrent)
	setRecurrentConnectionProbability(params, 0.1);

	printParameters(params);

	setCustomFitnessFunction(params, modulityFitness, "MSE");

	chromo = runCGP(params, NULL, numGens);
	_reportModulityChromo(chromo);

	freeChromosome(chromo);
	freeParameters(params);

	int abcc;
	std::cin >> abcc;
}