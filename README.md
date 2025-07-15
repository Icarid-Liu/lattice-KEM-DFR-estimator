# Lattice KEM DFR Estimator

A lightweight tool for estimating the Decryption Failure Rate (DFR) in lattice-based Key Encapsulation Mechanisms (KEMs). This tool provides efficient analysis of cryptographic failure rates across various lattice-based schemes.


### Structure

`DFR_utils.sage`: Core utility functions for DFR estimation and analysis

`KEM_DFR_benchmark.ipynb`: Interactive Jupyter notebook providing comparative benchmarks across different lattice-based KEMs (DAWN, NEV, and Kyber)

### Prerequisites

- SageMath (for running `.sage` files)
- Jupyter Notebook (for the benchmark notebook)

### Usage

1. Load the utility functions:
   ```python
   load('DFR_utils.sage')
   ```

2. An example to estimate DFR of DAWN-alpha under NIST-I:
   ```python
    # DAWN-alpha NIST-I
    n = 512
    q = 769
    kg = 160
    kf = 64
    ks = 96
    ke = 160
    
    Ds = build_ternary_distribution(n, ks, ks)
    Df = build_ternary_distribution(n, kf, kf)
    De_dc = build_uniform_distribution(-3, 3)
    Dm = build_ternary_distribution(n, n//4, 0)
    De = build_ternary_distribution(n, ke, ke)
    De_e_dc = add_distribution(De, De_dc)
    Dg = build_ternary_distribution(n, kg, kg)
    
    Dfe = calculate_convolution_pdf(n, Df, De_e_dc)
    Dgs = calculate_convolution_pdf(n, Dg, Ds)
    Dfm = calculate_convolution_pdf(n, Df, Dm)
    Dgs_fe = add_distribution(Dgs, Dfe)
    Dgs_fe2 = add_distribution(Dgs_fe, Dgs_fe)
    Dgs_fe2_fm = add_distribution(Dgs_fe2, Dfm)
    
    _, P1 = calculate_DRF(Dgs_fe2_fm, q//2, n, 1)
    
    Dz = Dgs_fe2_fm
    Dz2 = add_distribution(Dz, Dz)
    _, P2 = calculate_DRF(Dz2, q, 3 * n, 0)
    
    P = P1 + P2
    print('log(P1, 2) =', RR(log(P1, 2)))
    print('log(P2, 2) =', RR(log(P2, 2)))
    print('log(DFR, 2) =', RR(log(P, 2)))
   ```

## Acknowledgments

This project builds upon foundational work in lattice-based cryptography analysis. Some functions are adapted from the [leaky-LWE-Estimator](https://github.com/lducas/leaky-LWE-Estimator) repository.
