# Asymptotics of the Coefficients of $\prod_{k\ge1}(1+kx^k)$ (OEIS A022629)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24528285.svg)](https://doi.org/10.5281/zenodo.24528285)

**Preprint (Zenodo DOI):** [https://doi.org/10.5281/zenodo.24528285](https://doi.org/10.5281/zenodo.24528285)

Let $a(n)=[x^n]\prod_{k\ge1}(1+kx^k)$ ([A022629](https://oeis.org/A022629)), the sum over all partitions of $n$ into distinct parts of the product of the parts. Vaclav Kotesovec conjectured in 2018 that $\log a(n)\sim\sqrt{n/2}\,(\log(2n)-2)$. The manuscript proves

$$\log K!\le\log a(n)\le\log K!+O\Bigl(\frac{\sqrt n}{\log n}\Bigr),\qquad \log a(n)=\sqrt{\tfrac n2}\,\bigl(\log(2n)-2\bigr)+O\Bigl(\frac{\sqrt n}{\log n}\Bigr),$$

where $K$ is the largest integer with $K(K+1)/2\le n$, and the full asymptotic formula

$$a(n)=\frac{e^{\Phi_n(s_n)}}{\sqrt{2\pi V(s_n)}}\Bigl(1+O\bigl(n^{-1/4}\sqrt{\log n}\bigr)\Bigr),\qquad \Phi_n(s)=sn+\sum_{k\ge1}\log\bigl(1+ke^{-sk}\bigr),$$

where $s_n$ is the minimiser of $\Phi_n$ and $V(s)=\sum_k k^2p_k(1-p_k)$ with $p_k=ke^{-sk}/(1+ke^{-sk})$. The proof is a local limit theorem for a sum of independent Bernoulli variables.

## Contents

- `paper.tex`, `paper.pdf`: the manuscript.
- `verify_a022629.py`: computes $a(n)$ exactly, compares it with the OEIS data and the b-file (terms 0 to 10000), and evaluates both theorems.

## Reproducing

```bash
# standard library only; downloads b022629.txt from the OEIS once
python3 verify_a022629.py      # prints the table of the paper and ALL CHECKS PASSED (about 2 s)
```

## Status

Preprint, not yet peer-reviewed. Parts of the derivations and code were developed with AI assistance; every computational claim is backed by the script above.

## License

Manuscript text: CC BY 4.0. Code: MIT License.
