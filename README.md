# Filter

**A zero-dependency C++ framework for non-linear dynamic estimation and inverse problems.**

Define a forward model. Plug it in. Recover its hidden state from noisy measurements using numerically stable Extended or Unscented Kalman filtering — without pulling in a single external library.

---

## Motivation

A surprising range of problems in science and engineering reduce to the same formulation: given a (non-linear) forward model

$$y = F(x, \theta),$$

noisy measurements $\hat{y} \approx F(x_0, \theta)$, and a set of known parameters $\theta$, recover a good estimate of the hidden state $x_0$.

The same mathematical scaffolding that fits a polynomial to scattered points also inverts airborne electromagnetic (AEM) survey data over hundreds of square kilometres. Re-implementing estimation algorithms for each new problem is wasteful — and, in practice, it's where the numerical pitfalls live. **Filter** factors the problem cleanly into two pieces: a small, stable, well-tested estimation core, and a user-defined forward model that plugs into it through a thin abstract interface.

## Why this framework

- **Model-agnostic by design.** Every forward model derives from a single abstract base class. The filters know nothing about your model beyond the interface — swap from a polynomial to a three-dimensional geophysical conductivity model without touching a line of estimator code.
- **Numerically careful EKF.** The Extended Kalman Filter is implemented with two stability-oriented refinements that matter on real data:
  - **Sequential (component-by-component) measurement processing**, which sidesteps the inversion of the full innovation covariance matrix.
  - **Cholesky factorization of the covariance matrix $P$**, propagated in factored form rather than as $P$ itself. This preserves positive-definiteness and dramatically improves conditioning on long runs and ill-posed problems, where the textbook formulation tends to blow up.
- **Unscented Kalman Filter** for strongly non-linear models where EKF's Jacobian linearization becomes a liability.
- **Zero external dependencies.** Pure C++. No Eigen, no BLAS, no Boost, no Python glue. Clone, build with CMake, deploy anywhere a C++ toolchain runs.
- **Validated on real-world geophysics.** Ships with a full airborne-electromagnetics inversion pipeline that has produced results in three peer-reviewed publications (below). It has been exercised on noisy field data at scale, not just on synthetic benchmarks.

## Repository layout

```
filter/
├── types/       Abstract interfaces — Model, Filter, and the basic types
│                that algorithms and models agree on.
├── core/        Estimator implementations: EKF (with sequential update and
│                Cholesky-factored covariance), UKF, and supporting routines.
├── utilities/   Self-contained math and I/O helpers — FFT, spline routines,
│                text parsing for data ingestion.
├── examples/    Concrete models and end-to-end data-processing pipelines:
│   ├── polynomial/   a minimal toy model — the shortest path from a fresh
│   │                 clone to a working estimate.
│   └── aem/          airborne electromagnetics inversion — the real-world
│                     application that the published work is built on.
└── CMakeLists.txt
```

The separation between `types/` and `core/` is deliberate: the abstract interfaces are small and stable, while the algorithm side is free to grow. New estimators can be added without perturbing any existing model code.

## Quick start

```bash
git clone https://github.com/KhliustovDmitrii/filter.git
cd filter
cmake -B build
cmake --build build
```

Run the examples:

```bash
./build/examples/polynomial   # polynomial-fit demonstration
./build/examples/aem          # AEM inversion demonstration
```

No dependencies to install. A standard modern C++ compiler and CMake are all that is required.

## Defining your own model

A new estimation problem is, by design, a new class. The abstract `Model` in `types/` specifies the forward map and whatever derivatives your chosen estimator needs; user code supplies the body:

```cpp
#include "types/model.h"

class My_Model : public Model
{
public:

    // forward function y = F(x, theta)
    virtual void response(std::vector<double> &resp_arr) const override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) const override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) const override;
};
```

Instantiate a filter from `core/`, feed it the model and the measurement stream, and read out the estimate. The two examples in `examples/` — the polynomial fit and the AEM inversion — are the recommended references for what a full pipeline looks like end-to-end.

## Published work

The framework is the computational backbone of the author's research in dynamic estimation. It has been validated on non-trivial geophysical inverse problems in the following papers:

- Karshakov & Khliustov (2023). *Inversion of frequency-domain airborne electromagnetic data.* [PDF](https://www.geotechnologies.ru/publications/2023_AEM_kar_khliu.pdf)
- Khliustov et al. (2025). *Induced polarization effects in frequency-domain AEM data.* [PDF](https://www.geotechnologies.ru/publications/2025_11_Induced%20polarization%20effects%20in%20frequency%20domain%20AEM%20data.pdf)
- Khliustov et al. (2025). *The impact of AEM system configuration on inversion results.* [PDF](https://www.geotechnologies.ru/publications/2025_9_The%20impact%20of%20AEM%20system%20configuration%20on%20inversion%20results.pdf)

## Roadmap

- **SVD-based data denoising.** A pre-processing stage for noisy field data, improving the robustness of downstream inversion.
- **Gamma-spectrometry model.** Extending the framework's reach to a second geophysical modality, exercising the abstract model interface on a fundamentally different forward physics.

## License

Apache License 2.0 — see [`LICENSE`](LICENSE).

## Citation

If this framework or its algorithms support your work, please cite the publications listed above and this repository.

## Contact

Maintained by [Dmitrii Khliustov](https://github.com/KhliustovDmitrii). Issues and pull requests are welcome.
