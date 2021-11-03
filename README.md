Name: Akhil Shajan (Self Project)


OpenMP parallelization of QUICK program

QUICK is an open source, GPU enabled, ab initio and density functional theory program developed by Götz lab at University of California San Diego and Merz lab at Michigan State University. It has different features like Hartree-Fock energy calculations, Density functional theory calculations (LDA, GGA and Hybrid-GGA functionals available), Gradient and geometry optimization calculations and Mulliken charge analysis. It also supports QM/MM calculations with Amber21.
Fortran API is used in QUICK as QM energy and force engine and also MPI parallelization for CPU platforms is also implemented. Massively parallel GPU implementation via CUDA for Nvidia GPUs and Multi-GPU support via MPI + CUDA, also across multiple compute nodes is achieved.
Even though MPI+CUDA gives a very good speed, it is only achieved across multiple node and GPUs. But in cases where user with single node with multiple cores and multiple GPUs and the cores would fight for access to the memory at the same time. By sharing memory, OpenMP threads can reduce the memory footprint of the application, that means faster working code.
OpenMP is an Application Program Interface (API), jointly defined by a group of major computer hardware and software vendors. OpenMP provides a portable, scalable model for developers of shared memory parallel applications. The API supports C/C++ and Fortran on a wide variety of architectures. OpenMP is designed for multi-processor/core, shared memory machines. OpenMP programs accomplish parallelism exclusively through the use of threads. A thread of execution is the smallest unit of processing that can be scheduled by an operating system. The idea of a subroutine that can be scheduled to run autonomously might help explain what a thread is. Threads exist within the resources of a single process. Without the process, they cease to exist. Typically, the number of threads match the number of machine processors/cores. However, the actual use of threads is up to the application.



![image](https://user-images.githubusercontent.com/87541836/140169421-25453492-8efa-46d6-a0be-616bc8d4f1a4.png)
