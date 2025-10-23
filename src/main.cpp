
#include "discretizations/PdeDiscretization.hpp"
#include "domains/Real.hpp"

int main() {
    auto a = PdeDiscretization<Real>(4, 4);

    auto strrep = a.to_json_string();
    std::cout << strrep << std::endl;

    return 0;
}
