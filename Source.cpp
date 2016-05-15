#include <tuple>
#include <iostream>

using namespace std;

int main() {
	std::tuple<int, double, double> tt;

	std::get<0>(tt) = 23;
	std::get<1>(tt) = 23.0;
	std::get<2>(tt) = 24.0;

	std::cout << std::get<int>(tt);

	int abc;
	std::cin >> abc;
}