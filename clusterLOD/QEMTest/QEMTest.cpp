#include "../HalfEdgeMesh.hpp"
#include <iostream>
#include <cassert>

#include <glm/gtx/string_cast.hpp>

void QEMTestcase() {

}

int main() {
    std::cout << "Running QEMTest..." << std::endl;
    try {
        QEMTestcase();
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Exception caught: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "[ERROR] Unknown exception caught." << std::endl;
        return 1;
    }
    return 0;
}
