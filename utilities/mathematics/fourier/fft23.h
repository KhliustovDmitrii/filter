#ifndef FFT23_H
#define FFT23_H

#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstring>

namespace filter::math
{
class Fourier {
public:
    explicit Fourier(size_t n_request)
        : spec_len(1),
          order2(0),
          order3(0)
    {
        // Factor out powers of 2, then powers of 3
        size_t nn = n_request;
        while (nn > 0 && (nn % 2 == 0)) {
            nn /= 2;
            ++order2;
            spec_len *= 2;
        }
        while (nn > 0 && (nn % 3 == 0)) {
            nn /= 3;
            ++order3;
            spec_len *= 3;
        }
        // spec_len = 2^order2 * 3^order3  (largest factor-of-{2,3} portion of n_request)

        // Allocate (inv_index needs 2*n for the ping-pong buffer)
        inv_index.assign(2 * spec_len, 0);
        wkn.assign(spec_len, {0.0, 0.0});
        xn.assign(spec_len, {0.0, 0.0});
        fn.assign(spec_len, {0.0, 0.0});

        // Build digit-reversed permutation
        calculate_inv_ind();

        // Precompute twiddle factors
        const double angle_step = -2.0 * M_PI / static_cast<double>(spec_len);
        for (size_t i = 0; i < spec_len; ++i) {
            double angle = angle_step * static_cast<double>(i);
            wkn[i] = std::complex<double>(std::cos(angle), std::sin(angle));
        }
    }

    size_t size() const { return spec_len; }

    void set_xn(const std::vector<std::complex<double>>& xn_) {
        if (xn_.size() != spec_len)
            throw std::logic_error("FFT order mismatch");
        xn = xn_;
    }

    void set_xn(std::vector<std::complex<double>>&& xn_) {
        if (xn_.size() != spec_len)
            throw std::logic_error("FFT order mismatch");
        xn = std::move(xn_);
    }

    const std::vector<std::complex<double>>& get_fn() const { return fn; }

    void compute(bool inv);

    std::vector<std::complex<double>> xn;         // time-domain input

private:
    size_t spec_len;              // transform size  = 2^order2 * 3^order3
    size_t order2;                // number of radix-2 factors
    size_t order3;                // number of radix-3 factors

    std::vector<int>                  inv_index;  // 2*spec_len (ping-pong buffer)
    std::vector<std::complex<double>> wkn;        // twiddle factors
    std::vector<std::complex<double>> fn;         // frequency-domain result

    void calculate_inv_ind();
};
} // filter::math
#endif
