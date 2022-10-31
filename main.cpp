#include<vector>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
#include<map>
#include<exception>

int nakop = 1;

double null(double number) {
	if (number == -0) {
		number = 0;
		return number;
	}
	return number;
}

void printer(std::vector<std::vector<double>>& matr, const std::vector<std::string>& basis,
	const std::vector<std::string>& free) {
	std::cout << "\n\n\t ";
	for (size_t i = 0; i < free.size(); ++i) {
		std::cout << free[i] << "	 ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < basis.size(); ++i) {
		std::cout << basis[i] << "	";
		for (size_t j = 0; j < matr[i].size(); ++j) {
			if (matr[i][j] == -0) {
				matr[i][j] = 0;
			}
			std::cout << std::setw(5) << std::fixed << std::setprecision(2) << round(matr[i][j] * 100) / 100 << std::right << "\t";
		}
		std::cout << "\n";
	}
}

int find_column(std::vector<double>& F) {
	int index_r_column = -1;
	for (size_t j = 1; j < F.size(); ++j) {
		if (F[j] > 0) {
			index_r_column = j;
			break;
		}
	}
	return index_r_column;
}

int find_row(std::vector<std::vector<double>>& sympl, int index_r_column) {
	int index_r_row = -1;
	double min = INT64_MAX;
	for (size_t i = 0; i < sympl.size() - 1; ++i) {
		if (sympl[i][0] >= 0 && sympl[i][index_r_column] < 0) {
			continue;
		}
		double znach = sympl[i][0] / sympl[i][index_r_column];
		if (znach >= 0 && znach < min) {
			min = znach;
			index_r_row = i;
		}
	}

	return index_r_row;
}

void transformation(std::vector<std::vector<double>>& sympl, int index_r_row, int index_r_column) {
	double razr = sympl[index_r_row][index_r_column];
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (i != index_r_row) {
			for (size_t j = 0; j < sympl[i].size(); ++j) {
				if (j != index_r_column) {
					sympl[i][j] = sympl[i][j] - (sympl[index_r_row][j] * sympl[i][index_r_column] / razr);
				}
			}
		}
	}
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (i != index_r_row) {
			sympl[i][index_r_column] = -1 * sympl[i][index_r_column] / razr;
		}
	}
	for (size_t j = 0; j < sympl[index_r_row].size(); ++j) {
		if (j != index_r_column) {
			sympl[index_r_row][j] = sympl[index_r_row][j] / razr;
		}
	}
	sympl[index_r_row][index_r_column] = 1 / razr;
}

void searher(std::vector<std::vector<double>>& matr, std::vector<std::string>& free, std::vector<std::string>& basis) {
	int index_r_column;
	for (size_t i = 0; i < matr.size() - 1; ++i) {
		if (matr[i][0] < 0) {
			bool otr = false;
			for (size_t j = 1; j < matr[i].size(); ++j) {
				if (matr[i][j] < 0) {
					index_r_column = j;
					size_t index_r_raw = find_row(matr, index_r_column);
					transformation(matr, index_r_raw, index_r_column);
					std::swap(free[index_r_column], basis[index_r_raw]);
					otr = true;
					break;
				}
			}
			if (otr) {
				printer(matr, basis, free);
				i = static_cast<size_t>(-1);
			}
			else {
				throw std::logic_error("No sollution to this task :(");
			}
		}
	}
}

void method(std::vector<std::vector<double>>& matr, std::vector<std::string>& free,
	std::vector<std::string>& basis, std::string rezult) {

	if (rezult == "min") {
		for (auto& el : matr[matr.size() - 1]) {
			el = -el;
		}
	}

	std::map<std::string, double> m_result;

	for (size_t i = 1; i < free.size(); ++i) {
		m_result[free[i]] = 0;
	}

	printer(matr, basis, free);
	searher(matr, free, basis);

	int index_r_column = find_column(matr[matr.size() - 1]);

	while (index_r_column != -1) {
		int index_r_row = find_row(matr, index_r_column);
		if (index_r_row == -1) {
			std::cout << R"(Infinity number of sollution \(-_-)/)";
			break;
		}
		std::swap(basis[index_r_row], free[index_r_column]);
		transformation(matr, index_r_row, index_r_column);
		printer(matr, basis, free);
		index_r_column = find_column(matr[matr.size() - 1]);
	}
	std::cout << "Answer: " << std::endl;
	for (size_t i = 0; i < basis.size() - 1; ++i) {
		m_result.at(free[i]) = matr[i][0];
	}
	for (const auto& el : m_result) {
		std::cout << el.first << " = " << el.second << std::endl;
	}
	if (rezult == "max") {
		matr[matr.size() - 1][0] *= -1;
	}
	std::cout << "F = " << matr[matr.size() - 1][0] << std::endl;
}

class Simplex_tabels {
	std::vector<std::string> basis;
	std::vector<std::string> free;
	std::vector<std::vector<double>> matr;
	std::string rezult;
	std::map<std::string, double> m_result;
	std::map<std::string, double> the_int_result;
	double promezutochniy_result_func;
	int maximum = INT16_MIN;

public:

	Simplex_tabels(std::vector<std::string>& _basis, std::vector<std::string>& _free,
		std::vector<std::vector<double>>& _matr, std::string _rezult) {
		basis = _basis;
		free = _free;
		matr = _matr;
		rezult = _rezult;
	}

	void Simpl_method() {

		if (rezult == "min") {
			for (auto& el : matr[matr.size() - 1]) {
				el = -el;
			}
		}

		for (size_t i = 1; i < free.size(); ++i) {
			m_result[free[i]] = 0;
		}

		printer(matr, basis, free);
		searher(matr, free, basis);

		int index_r_column = find_column(matr[matr.size() - 1]);

		while (index_r_column != -1) {

			int index_r_row = find_row(matr, index_r_column);

			if (index_r_row == -1) {
				throw std::logic_error("Infinity number of sollution \(-_-)/");
			}

			std::swap(basis[index_r_row], free[index_r_column]);

			transformation(matr, index_r_row, index_r_column);

			printer(matr, basis, free);

			index_r_column = find_column(matr[matr.size() - 1]);
		}
		std::cout << "Answer:\n";
		for (size_t i = 0; i < basis.size() - 1; ++i) {
			try {
				m_result.at(basis[i]) = matr[i][0];
			}
			catch (std::exception&) {}
		}
		for (const auto& el : m_result) {
			std::cout << el.first << " = " << el.second << std::endl;
		}
		promezutochniy_result_func = matr[matr.size() - 1][0];
		if (rezult == "max") {
			promezutochniy_result_func *= -1;
		}
		std::cout << "F = " << promezutochniy_result_func << std::endl;
	}

	void add_new_less_row(char index_of_var, int new_limit) {
		std::vector<double> new_row(free.size());
		new_row[0] = new_limit;
		int index = static_cast<int>(index_of_var) - 48;
		new_row[index] = 1;
		basis.insert(basis.end() - 1, "x" + std::to_string(basis.size() + free.size() - 1));
		matr.insert(matr.end() - 1, new_row);
	}

	void add_new_greater_row(char index_of_var, int new_limit) {
		std::vector<double> new_row(free.size());
		new_row[0] = -new_limit;
		int index = static_cast<int>(index_of_var) - 48;
		new_row[index] = -1;
		basis.insert(basis.end() - 1, "x" + std::to_string(basis.size() + free.size() - 1));
		matr.insert(matr.end() - 1, new_row);
	}

	Simplex_tabels& operator=(const Simplex_tabels& ob) {
		matr = ob.matr;
		free = ob.free;
		basis = ob.basis;
		rezult = ob.rezult;
		return *this;
	}

	friend void method_BUB(Simplex_tabels& ob);
};

void method_BUB(Simplex_tabels& ob) {
	Simplex_tabels source = ob;
	try {
		ob.Simpl_method();
	}
	catch (std::exception& error) {
		std::cout << std::endl << error.what() << std::endl;
		return;
	}

	double E = 1e-6;
	std::pair<std::string, double> branching_var = { std::string(), 0 };
	for (auto& el : ob.m_result) {
		if (fabs(el.second - std::round(el.second)) > E) {
			branching_var = el;
			break;
		}
	}

	if (!branching_var.first.empty()) {
		ob = source;
		std::cout << "branching_variable: " << branching_var.first << std::endl;
		std::cout << branching_var.first << " >= " << ceil(branching_var.second) << std::endl;

		ob.add_new_greater_row(branching_var.first[1], ceil(branching_var.second));
		method_BUB(ob);

		std::cout << std::endl << branching_var.first << " <= " << floor(branching_var.second) << std::endl;

		ob = source;
		ob.add_new_less_row(branching_var.first[1], floor(branching_var.second));
		method_BUB(ob);

		ob = source;
	}
	else {
		if (ob.maximum < ob.promezutochniy_result_func) {
			ob.maximum = static_cast<int>(ob.promezutochniy_result_func);
			ob.the_int_result = ob.m_result;
		}
	}
	++nakop;
	if (nakop == 4) {
		std::cout << std::endl;
		for (auto& t : ob.the_int_result) {
			std::cout << t.first << "=" << t.second << "\t";
		}
		std::cout << "F=" << ob.maximum << std::endl << std::endl;
	}
}

int counter(std::vector<int>& variables, std::vector<double>& ex) {
	int res = 0;
	for (size_t i = 0; i < variables.size(); ++i) {
		res += static_cast<int>(variables[i] * ex[i + 1]);
	}
	return res;
}
bool logic(std::vector<int>& variables, std::vector<std::vector<double>>& matr) {
	for (size_t i = 0; i < variables.size(); ++i) {
		int res = counter(variables, matr[i]);
		if (res > matr[i][0]) {
			return false;
		}
	}
	return true;
}

void bruteforce(std::vector<std::vector<double>> matr) {
	std::cout << "Bruteforce: \n";
	auto F = matr[matr.size() - 1];
	std::vector<int> variables = { 0, 0 ,0 };
	matr.pop_back();
	while (logic(variables, matr)) {
		while (logic(variables, matr)) {
			while (logic(variables, matr)) {
				std::cout << "[";
				for (size_t i = 0; i < variables.size(); ++i) {
					std::cout << variables[i] << " ";
				}
				std::cout << "] = ";
				std::cout << counter(variables, F) << std::endl;
				++variables[2];
			}
			variables[2] = 0;
			++variables[1];
		}
		variables[1] = 0;
		++variables[0];
	}
}

int main() {
	std::vector<std::vector<double>> matr = {
		{3,1,1,1},
		{5,1,4,0},
		{7,0,0.5,3},
		{0,3,3,7}
	};
	std::vector<std::string> free = { "sv","x1","x2","x3" };
	std::vector<std::string> basis = { "x4","x5","x6","F" };
	std::string rezult = "max";

	Simplex_tabels ob(basis, free, matr, rezult);

	method_BUB(ob);
	bruteforce(matr);
	return 0;
}
