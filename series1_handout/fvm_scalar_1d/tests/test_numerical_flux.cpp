#include <gtest/gtest.h>

#include <ancse/numerical_flux.hpp>

template <class NumericalFlux>
void check_consistency(const NumericalFlux &nf) {
    auto model = make_dummy_model();

    auto to_check = std::vector<double>{1.0, 2.0, 3.3251, -1.332};
    for (auto u : to_check) {
        ASSERT_DOUBLE_EQ(model.flux(u), nf(u, u));
    }
}

TEST(TestCentralFlux, consistency) {
    auto model = make_dummy_model();
    auto central_flux = CentralFlux(model);

    check_consistency(central_flux);
}


