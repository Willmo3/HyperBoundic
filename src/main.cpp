
#include <fstream>

#include "meshes/RectangularMesh.hpp"
#include "domains/Real.hpp"
#include "Wixed/WixedForm.hpp"

int main() {
    std::ofstream file("out.json");

    auto a = RectangularMesh<WixedForm>(4, 4);
    file << a.to_json_string();
    file.close();

    return 0;
}
