# Krylov Preconditioners Benchmark

Порівняння ітеративних методів Крилова (CG, BiCGSTAB, GMRES тощо) та різних передобумовлювачів (Jacobi, SSOR, ILU, AMG) для розв'язання систем лінійних рівнянь.

## Особливості
- Власні обгортки `Vector` та `Matrix` над Eigen
- Абстрактний інтерфейс для солверів і прекондиціонерів
- Підтримка PETSc та Trilinos (планується)
- Детальне логування часу, ітерацій та похибок

## Зібрано методи
- **Солвери**: CG, BiCG, BiCGSTAB, GMRES(m)
- **Прекондиціонери**: NoPreconditioner, Jacobi, SSOR, ILU (в процесі)

## Технології
- C++20
- Eigen3
- PETSc / Trilinos (в процесі)

## Збірка

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
Автор
Оля Богданенко
