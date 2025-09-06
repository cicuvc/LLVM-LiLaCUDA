# LLVM CUDA compiler with linear layout builtins

This repository contians the source code of a modified version of LLVM/clang with linear layout solvers/miscellaneous/builtins for CUDA template libraries. 

Supported clang builtins:

1. Swizzle solve

    ```c++
    __builtin_swizzle_solve(uint32_t n, uint32_t r, <patterns>...) -> uint32_t
    ```
    For shared memory with size $2^{n}$, assuming threads use 16 bytes aligned vectorized instruction, the swizzle matrix on $F_2$ maps flatten tensor address of $n$ bits to $n$ bits shared memory address, while the first 3 rows (sub-matrix $S$) determine if conflicts exist. Access patterns can be represented by $n\times 3$ matrices $P_1,...,P_m$，requiring that for each $P_i$， $S P_i$ is invertible. This builtin solves matrix $S$ for given patterns while trying to minimize the number of non-zero diagonal lines in $S$ to simplify swizzle calculation. The builtin function produce the r-th row of solution. If no solution for given patterns, it returns 0.

    Supported backends:

    - Rule based: Non-general solver. Produce constructed solution for simple cases where all $P_i$ are formed by 3 consecutive one-hot column vector. 

    - Baseline solver: General solver with $\mathcal{O}(2^{n})$ worst-case time, guarantee solution found if exists, may not be optimal.

    - (In progress) Search-based solver: Memorized bi-direction search with state-compression, $\mathcal{O}(2^{4m}\cdot n)$ time, acceptable when $m\le5$, solution is guaranteed to be optimal.

    - (In progress) SA solver: Simulated Annealing solver, guarantee solution found if exists, may not be optimal but better than baseline solver.
    
2. Intra-warp shuffle decomposition

   decompose layout transition matrix into composition of 2 intra-thread register shuffles and 1 inter-thread shuffles.
   
   Not implemented yet. Interface design is under revision .

3. Miscellaneous

   - builtin GEMV on $F_2$ domain
     
     ```c++
     __builtin_f2_gemv(uint32_t vec, uint32_t <matrix rows>...) -> uint32_t
     ```

Enable builtins with compiler option `-fcuda-linear-layout-extensions`

-----

# The LLVM Compiler Infrastructure

[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/llvm/llvm-project/badge)](https://securityscorecards.dev/viewer/?uri=github.com/llvm/llvm-project)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/8273/badge)](https://www.bestpractices.dev/projects/8273)
[![libc++](https://github.com/llvm/llvm-project/actions/workflows/libcxx-build-and-test.yaml/badge.svg?branch=main&event=schedule)](https://github.com/llvm/llvm-project/actions/workflows/libcxx-build-and-test.yaml?query=event%3Aschedule)

Welcome to the LLVM project!

This repository contains the source code for LLVM, a toolkit for the
construction of highly optimized compilers, optimizers, and run-time
environments.

The LLVM project has multiple components. The core of the project is
itself called "LLVM". This contains all of the tools, libraries, and header
files needed to process intermediate representations and convert them into
object files. Tools include an assembler, disassembler, bitcode analyzer, and
bitcode optimizer.

C-like languages use the [Clang](https://clang.llvm.org/) frontend. This
component compiles C, C++, Objective-C, and Objective-C++ code into LLVM bitcode
-- and from there into object files, using LLVM.

Other components include:
the [libc++ C++ standard library](https://libcxx.llvm.org),
the [LLD linker](https://lld.llvm.org), and more.

## Getting the Source Code and Building LLVM

Consult the
[Getting Started with LLVM](https://llvm.org/docs/GettingStarted.html#getting-the-source-code-and-building-llvm)
page for information on building and running LLVM.

For information on how to contribute to the LLVM project, please take a look at
the [Contributing to LLVM](https://llvm.org/docs/Contributing.html) guide.

## Getting in touch

Join the [LLVM Discourse forums](https://discourse.llvm.org/), [Discord
chat](https://discord.gg/xS7Z362),
[LLVM Office Hours](https://llvm.org/docs/GettingInvolved.html#office-hours) or
[Regular sync-ups](https://llvm.org/docs/GettingInvolved.html#online-sync-ups).

The LLVM project has adopted a [code of conduct](https://llvm.org/docs/CodeOfConduct.html) for
participants to all modes of communication within the project.
