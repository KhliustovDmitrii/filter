#include "fft23.h"

namespace filter::math
{
void Fourier::calculate_inv_ind()
{
    const int n = static_cast<int>(spec_len);
    const int deo = (order2 == 0) ? 1 : 0;

    // Identity permutation in the first half
    for (int i = 0; i < n; ++i)
        inv_index[i] = i;

    int num   = 1;
    int order = 0;
    int k = 0, l = 0;          // will track which half is source / dest

    // Radix-3 digit-reversal passes
    for (; order < static_cast<int>(order3) - deo; ++order) {
        l = (order & 1);
        k = !l;
        if (k) k *= n;
        else   l *= n;
        int m = n / num - 1;
        for (int i = 0; i < num; ++i) {
            int offs = i * (m + 1);
            for (int j = 0; j < m; ++j) {
                inv_index[j + offs + k] = inv_index[(j * 3) % m + offs + l];
            }
            inv_index[m + offs + k] = inv_index[m + offs + l];
        }
        num *= 3;
    }

    // Radix-2 digit-reversal passes
    for (; order < static_cast<int>(order2 + order3) - 1; ++order) {
        l = (order & 1);
        k = !l;
        if (k) k *= n;
        else   l *= n;
        int m = n / num - 1;
        for (int i = 0; i < num; ++i) {
            int offs = i * (m + 1);
            for (int j = 0; j < m; ++j) {
                inv_index[j + offs + k] = inv_index[(j * 2) % m + offs + l];
            }
            inv_index[m + offs + k] = inv_index[m + offs + l];
        }
        num *= 2;
    }

    // Copy result to the first half if the final pass left it in the second half
    if (k) {
        std::memcpy(inv_index.data(),
                    inv_index.data() + k,
                    static_cast<size_t>(k) * sizeof(int));
    }
}

void Fourier::compute(bool inv)
{
    static const double SQRT3_2 = std::sqrt(3.0) / 2.0;
    const int N = static_cast<int>(spec_len);

    // Helper: twiddle-factor index
    auto twiddle = [N](int k, int n_cur) -> int {
        return k * (N / n_cur);
    };

    // Reorder input by digit-reversed permutation
    for (int m = 0; m < N; ++m)
        fn[m] = xn[inv_index[m]];

    // Radix-2 butterfly stages
    int n = 1;
    for (size_t stage = 0; stage < order2; ++stage) {
        int nd2 = n;
        n += n;                         // n *= 2
        for (int k = 0; k < nd2; ++k) {
            int ind = twiddle(k, n);
            for (int m = k; m < N; m += n) {
                int mpnd2 = m + nd2;
                std::complex<double> temp = wkn[ind] * fn[mpnd2];
                fn[mpnd2] = fn[m] - temp;
                fn[m]     = fn[m] + temp;
            }
        }
    }

    // Radix-3 butterfly stages
    for (size_t stage = 0; stage < order3; ++stage) {
        int nd2 = n;
        n += 2 * n;                     // n *= 3
        for (int k = 0; k < nd2; ++k) {
            int ind  = twiddle(k, n);
            int ind1 = twiddle(2 * k, n);
            for (int m = k; m < N; m += n) {
                int mpnd2 = m + 2 * nd2;   // third element
                int mpnd3 = m + nd2;        // second element
                std::complex<double> temp  = wkn[ind]  * fn[mpnd3];
                std::complex<double> temp1 = wkn[ind1] * fn[mpnd2];
                std::complex<double> temp2 = temp + temp1;          // sum
                temp = temp - temp1;                                // diff

                // i * diff * sqrt(3)/2
                temp1 = std::complex<double>(
                    -std::imag(temp) * SQRT3_2,
                     std::real(temp) * SQRT3_2
                );

                // half sum
                temp = std::complex<double>(
                    0.5 * std::real(temp2),
                    0.5 * std::imag(temp2)
                );

                fn[mpnd3] = fn[m] - temp - temp1;
                fn[mpnd2] = fn[m] - temp + temp1;
                fn[m]     = fn[m] + temp2;
            }
        }
    }

    // Post-processing (real-spectrum extraction)
    if (!inv) {
        for (int m = 1; m < N / 2; ++m) {
            fn[m] += std::conj(fn[N - m]);
            fn[N - m] = {0.0, 0.0};
            fn[m] = std::conj(fn[m]) / static_cast<double>(N);
        }
        fn[0]     = {0.0, 0.0};
        fn[N / 2] = {0.0, 0.0};
    }
}
} // filter::math