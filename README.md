# Krylov Solvers & Preconditioners Benchmark

**Дослідження** з порівняльного аналізу ітераційних методів Крилова та методів прекондиціонування для розв'язання розріджених систем лінійних алгебраїчних рівнянь (СЛАР).

Проєкт призначений для бенчмаркінгу та оцінки ефективності (час, пам'ять, ітерації, норми похибок) власних реалізацій алгоритмів у порівнянні з промисловими HPC-бібліотеками **PETSc** та **Trilinos**.

---

## Особливості та Можливості

*   **Мультиплатформенна інтеграція:** Робота з матрицями через власні обгортки `Vector` та `Matrix` над **Eigen3**, з експортом у структури PETSc та Tpetra (Trilinos).
*   **Крос-фреймворкові адаптери (Generic Framework Adapters):** Реалізовано спеціальні оболонки (`PCShell` для PETSc та `Tpetra::Operator` для Trilinos), які дозволяють підключати **кастомні (і не їхні рідні) прекондиціонери** безпосередньо до солверів сторонніх промислових бібліотек.
*   **Профілювання:** Моніторинг пікового та поточного споживання RAM (`MemoryTracker`) для оцінки overhead кожного прекондиціонера та метода.
*   **Аналіз збіжності:** Логування повної історії зміни норми нев'язки на кожній ітерації для побудови графіків.

---

## Підтримувані методи / прекондиціонери

### 1. Солвери (Krylov Subspace Methods)
*   **Власні реалізації (C++20):** CG, BiCG, BiCGSTAB, GMRES(m)
*   **Промислові ядра:** PETSc KSP & Trilinos Belos

### 2. Прекондиціонери (Preconditioners)

| Метод | Власна реалізація (Eigen) | PETSc (PC) | Trilinos (Ifpack2) |
| :--- | :---: | :---: | :---: |
| **None** | ✅ | ✅ | ✅ |
| **Jacobi** | ✅ | ✅ | ✅ |
| **Gauss-Seidel** | ✅ | ✅ | ✅ |
| **SOR / SSOR** | ✅ | ✅ | ✅ |
| **ILU(0) / ILU(p)** | ❌ | ✅ | ✅ (RILUK) |
| **ILU(p, $\tau$)** | ❌ | ✅ | ✅ (ILUT) |
| **ASM / Schwarz** | ❌ | ✅ | ✅ |
| **Chebyshev** | ❌ | ❌ | ✅ |

---

## Технологічний стек

*   **Мова:** C++20 (стандарт активовано)
*   **Керування пакетами:** Spack (для зручного керування бібліотеками)
*   **Математичні бекенди:**
    *   Eigen3 (RowMajor, інтеграція MatrixMarket)
    *   PETSc
    *   Trilinos
    *   MPI & Kokkos *(в процесі інтеграції)*

---

## Планується

*   Розподіл на ядра
*   Деякі прекондиціонери з PC PETSc

---

## Приклад виводу
```text
rows: 81920 cols: 81920
compute matrices patch ../matrices/shallow_water1.mtx
Method    Precond         Status      Iter    Time solve    Time prec     Memory solve   Memory prec    Norm          Rel_error     
BiCG      P ILU(0)        Failed      13      2.87738e+00   3.09330e-02   1.25788e+07    8.83507e+06    1.71521e-03   2.27285e-16   
BiCG      T ILU(0)        Failed      13      3.35912e+00   5.98088e-01   1.29106e+07    1.56754e+07    1.71521e-03   2.27466e-16   
BiCG      P ILU(p)        Success     7       1.48047e+00   5.44520e-02   9.94099e+06    2.92864e+06    3.38839e-04   1.93523e-16   
BiCG      T ILU(p)        Success     7       1.71297e+00   8.28376e-01   1.31072e+06    1.47948e+07    3.38839e-04   1.93241e-16   
BiCG      P ILU(p,tau)    Success     7       1.48838e+00   4.62961e-02   0.00000e+00    0.00000e+00    3.38839e-04   1.93523e-16   
BiCG      T ILU(p, tau)   Failed      8       2.11110e+00   6.17541e-01   0.00000e+00    1.04079e+07    8.61038e-03   2.21130e-16   
BiCG      P ASM           Failed      13      2.91430e+00   4.43272e-02   1.31072e+05    9.17504e+05    1.71521e-03   2.27285e-16   
BiCG      P ICC           Failed      13      2.87808e+00   1.89920e-02   6.55360e+04    2.66240e+05    1.71521e-03   2.27089e-16   
BiCG      T ILU(k)        Failed      8       2.10539e+00   6.04738e-01   0.00000e+00    8.84736e+05    8.56969e-03   2.21277e-16   
BiCG      T Chebyshev     Failed      45      1.13542e+01   4.05294e-01   1.31072e+05    1.01581e+06    1.66323e-03   3.49073e-16   
BiCG      T Schwarz(ilut) Failed      13      3.46366e+00   5.63107e-01   4.58752e+05    3.68640e+05    1.71521e-03   2.26798e-16
```

---

<img width="2951" height="1754" alt="convergence_BiCG" src="https://github.com/user-attachments/assets/1a1af87e-74f8-49b9-947f-729b44489aa9" />

---

## Швидкий старт & Збірка

Для збірки проєкту на Linux/використовується менеджер пакетів **Spack** та **CMake**.

### 1. Встановлення залежностей (через Spack)
```bash
spack install petsc
spack install trilinos +ifpack2 +tpetra
spack install eigen
```
## Структура проєкту
```bash
├── CMakeLists.txt              # Конфігурація збірки та пошук пакетів через Spack
├── main.cpp                    # Контролер бенчмарку та конфігурація тестів
├── include/                    # Заголовні файли (.h)
	├── Preconditioners/        # Заголовні файли прекондиціонерів
	├── Solvers/                # Заголовні файли солверів
	└── Усі інші заголовні файли
├── src/                        # Реалізація логіки (.cpp)
    ├── Preconditioners/        # Власні, PETSc та Trilinos прекондиціонери
    ├── Solvers/                # Власні, PETSc та Trilinos вирішувачі
    ├── Matrix.cpp / Vector.cpp # Математичні обгортки
    ├── Logger.cpp              # Логування результатів / збереження у .json
    └── MemoryTracker.cpp       # Утиліта трекінгу RAM
└── matrices/                   # Матриці для дослідження у форматі .mtx
```

## Клонування та компіляція
```bash
git clone https://github.com/olhabohdanenko/cpp-krylov-solvers.git
cd cpp-krylov-solvers
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## Запуск
```bash
./SLAE_S
```
