
This is a repo of codes used for IPDPS 2023 best paper nominee listed
in https://www.ipdps.org/ipdps2023/2023-advance-program.html under
"PLENARY SESSION: Best Papers" on Tuesday, May 16, 2023.
The paper details are as follows: "PAQR: Pivoting Avoiding QR factorization" by
Wissam Sid-Lakhdar, Sebastien Cayrols, Daniel Bielich, Ahmad Abdelfattah,
Piotr Luszczek, Mark Gates, and Stanimire Tomov (University of Tennessee at Knoxville);
Hans Johansen and David Williams-Young (Lawrence Berkeley National Laboratory);
Timothy Davis (Texas A&M University);
and Jack Dongarra and Hartwig Anzt (University of Tennessee at Knoxville)
DOI: 10.1109/IPDPS54959.2023.00040
https://doi.ieeecomputersociety.org/10.1109/IPDPS54959.2023.00040

# Originally github.com/drslmw/paqr_SC22
Originally it was an artifact description for article evaluation for SC22

Four components are present to reproduce the results in the paper.

1) paqr_matlab.tar.gz
2) docker wissamsid/paqr_lapack
3) docker abdelfattah83/paqr_batch_cuda:latest
4) paqr_scalapack.tar.gz

## 1) MATLAB

In order to reproduce MATLAB results.

```
tar xz paqr_matlab.tar.gz
cd matlab/src
matlab
>>> test(100);
```

The "100" should be replaced by "1000" to get results of Table1.

NOTE: The version of MATLAB used is `R2020a Update 2 (9.8.0.1380330) 64-bit (glnxa64)`.

## 2) LAPACK

```
docker pull wissamsid/paqr_lapack
docker run -it --entrypoint=bash --privileged --userns=host wissamsid/paqr_lapack
sh test_paqr_lapack.sh 100 100
```

In order to reproduce the results of Table2, replace "100 100" by "10000 10000".

NOTE: The version of LAPACK used is `3.10.0`.

## 3) MAGMA

In order to reproduce the results of Table3, run the following:

```
docker pull abdelfattah83/paqr_batch_cuda:latest
docker run -it --entrypoint=bash --privileged --userns=host abdelfattah83/paqr_batch_cuda:latest
sh test_paqr_batch.sh
```

NOTE: The MAGMA version used is `2.5.0`.

## 4) ScaLAPACK

The following "module load" are examples only.  This only highlights that cmake, MKL and Open-MPI should be present in the environment.

```
tar xzf paqr_scalapack.tar.gz
module load cmake intel-oneapi-mkl openmpi
cd scalapack
bash test.sh
```

The five results lines correspond to: warmup QR 1, warmup QR 2, QR, PAQR, QRCP.

NOTE: The ScaLAPACK version used is `2.2`.

## Appendix) Software stack for LAPACK and ScaLAPACK

The ScaLAPACK code was compiled with:
- GCC/10.2.0
- Intel-oneapi-MKL/2022.0.2
- Open-MPI/4.1.2
- cmake/3.22.2
