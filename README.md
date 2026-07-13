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

Method         Precond         Status      Iter    Time solve    Time prec     Memory solve b Memory prec b  Norm          Rel_error     
C CG           P ILU(0)        Success     11      9.93030e-01   3.20274e-02   4.71859e+06    8.81869e+06    5.82331e-01   7.16815e-15   
C CG           T ILU(0)        Success     11      1.06377e+00   3.20865e-01   3.80109e+06    5.95558e+06    5.82331e-01   7.16795e-15   
C CG           P ILU(p)        Success     6       5.36609e-01   5.32260e-02   0.00000e+00    0.00000e+00    8.11318e-02   1.02244e-15   
C CG           T ILU(p)        Success     6       5.92350e-01   4.56329e-01   0.00000e+00    8.98662e+06    8.11318e-02   1.02232e-15   
C CG           P ILU(p,tau)    Success     6       5.51519e-01   5.36543e-02   0.00000e+00    0.00000e+00    8.11318e-02   1.02244e-15   
C CG           T ILU(p, tau)   Success     8       7.77224e-01   2.68130e-01   0.00000e+00    1.17883e+07    1.55118e-01   1.72429e-15   
C CG           P ASM           Success     11      1.02277e+00   4.48941e-02   1.31072e+05    9.17504e+05    5.82331e-01   7.16815e-15   
C CG           P ICC           Success     11      1.01264e+00   2.42484e-02   6.55360e+04    2.62144e+05    5.82331e-01   7.16695e-15   
C CG           T ILU(k)        Success     8       7.83839e-01   2.51233e-01   0.00000e+00    8.47872e+05    1.56296e-01   1.73942e-15   
C CG           T Chebyshev     Success     38      3.75656e+00   6.86388e-02   1.31072e+05    9.95328e+05    6.18374e-01   4.65842e-15   
C CG           T Schwarz(ilut) Success     11      1.13014e+00   2.18153e-01   3.93216e+05    4.05504e+05    5.82331e-01   7.16828e-15

Method         Precond         Status      Iter    Time solve    Time prec     Memory solve b Memory prec b  Norm          Rel_error     
P cg           P ILU(0)        Success     11      7.41694e-02   2.60598e-02   1.63840e+05    0.00000e+00    5.70227e-11   7.16732e-15   
P cg           T ILU(0)        Success     11      1.80179e-01   2.39386e-01   0.00000e+00    0.00000e+00    5.70227e-11   7.16632e-15   
P cg           P ILU(p)        Success     6       3.53529e-02   5.36298e-02   0.00000e+00    0.00000e+00    8.27325e-12   1.02278e-15   
P cg           T ILU(p)        Success     6       1.14309e-01   4.43004e-01   0.00000e+00    0.00000e+00    8.27325e-12   1.02301e-15   
P cg           P ILU(p,tau)    Success     6       3.55439e-02   5.28246e-02   0.00000e+00    0.00000e+00    8.27325e-12   1.02278e-15   
P cg           T ILU(p, tau)   Success     8       1.66225e-01   2.44753e-01   0.00000e+00    7.08608e+06    1.41124e-11   1.72577e-15   
P cg           P ASM           Success     11      8.81157e-02   4.87087e-02   0.00000e+00    0.00000e+00    5.70227e-11   7.16732e-15   
P cg           P ICC           Success     11      6.79996e-02   2.27596e-02   0.00000e+00    0.00000e+00    5.70227e-11   7.16692e-15   
P cg           T ILU(k)        Success     8       1.54880e-01   2.48470e-01   0.00000e+00    3.19488e+05    1.42230e-11   1.73968e-15   
P cg           T Chebyshev     Success     37      5.79564e-01   4.95400e-02   0.00000e+00    0.00000e+00    1.16854e-10   1.38484e-14   
P cg           T Schwarz(ilut) Success     11      2.46143e-01   2.16301e-01   0.00000e+00    4.21888e+05    5.70227e-11   7.16736e-15

Method         Precond         Status      Iter    Time solve    Time prec     Memory solve b Memory prec b  Norm          Rel_error     
T Block CG     P ILU(0)        Success     11      3.24532e-01   2.68790e-02   1.92102e+06    0.00000e+00    6.11074e-15   7.16713e-15   
T Block CG     T ILU(0)        Success     11      7.73366e-02   2.47639e-01   0.00000e+00    0.00000e+00    6.11074e-15   7.16622e-15   
T Block CG     P ILU(p)        Success     6       2.00165e-01   5.22212e-02   0.00000e+00    0.00000e+00    8.51363e-16   1.02340e-15   
T Block CG     T ILU(p)        Success     6       5.93273e-02   4.42613e-01   0.00000e+00    0.00000e+00    8.51363e-16   1.02411e-15   
T Block CG     P ILU(p,tau)    Success     6       2.13023e-01   5.58053e-02   0.00000e+00    0.00000e+00    8.51363e-16   1.02340e-15   
T Block CG     T ILU(p, tau)   Success     8       8.10579e-02   2.49835e-01   0.00000e+00    2.00294e+06    1.62774e-15   1.72562e-15   
T Block CG     P ASM           Success     11      3.52192e-01   4.80966e-02   0.00000e+00    0.00000e+00    6.11074e-15   7.16713e-15   
T Block CG     P ICC           Success     11      3.30685e-01   2.31016e-02   0.00000e+00    0.00000e+00    6.11074e-15   7.16676e-15   
T Block CG     T ILU(k)        Success     8       8.18560e-02   2.44390e-01   0.00000e+00    1.63021e+06    1.64011e-15   1.73909e-15   
T Block CG     T Chebyshev     Success     38      2.45969e-01   4.89320e-02   0.00000e+00    0.00000e+00    9.48525e-15   6.79351e-15   
T Block CG     T Schwarz(ilut) Success     11      1.48360e-01   2.04300e-01   0.00000e+00    3.31776e+05    6.11074e-15   7.16701e-15

Method         Precond         Status      Iter    Time solve    Time prec     Memory solve b Memory prec b  Norm          Rel_error     
T tfqmr        P ILU(0)        Success     6       3.49892e-01   2.59833e-02   1.31072e+05    0.00000e+00    7.92545e-15   2.58172e-15   
T tfqmr        T ILU(0)        Success     6       8.94229e-02   2.39531e-01   0.00000e+00    0.00000e+00    7.92545e-15   2.58034e-15   
T tfqmr        P ILU(p)        Success     3       2.05952e-01   5.12049e-02   0.00000e+00    0.00000e+00    8.89890e-15   3.40674e-15   
T tfqmr        T ILU(p)        Success     3       6.41265e-02   4.44986e-01   0.00000e+00    0.00000e+00    8.89890e-15   3.40694e-15   
T tfqmr        P ILU(p,tau)    Success     3       2.05034e-01   5.11419e-02   0.00000e+00    0.00000e+00    8.89890e-15   3.40674e-15   
T tfqmr        T ILU(p, tau)   Success     4       8.59944e-02   2.34634e-01   0.00000e+00    3.19488e+05    4.94956e-16   3.41361e-16   
T tfqmr        P ASM           Success     6       3.97113e-01   4.63414e-02   0.00000e+00    0.00000e+00    7.92545e-15   2.58172e-15   
T tfqmr        P ICC           Success     6       3.62767e-01   2.43373e-02   0.00000e+00    0.00000e+00    7.92545e-15   2.58028e-15   
T tfqmr        T ILU(k)        Success     4       8.66926e-02   2.50576e-01   0.00000e+00    5.57056e+05    4.92162e-16   4.16666e-16   
T tfqmr        T Chebyshev     Success     21      2.97275e-01   5.07859e-02   0.00000e+00    0.00000e+00    4.10102e-15   5.44251e-15   
T tfqmr        T Schwarz(ilut) Success     6       1.65954e-01   2.15996e-01   0.00000e+00    3.97312e+05    7.92545e-15   2.58280e-15
```

---

<img width="2951" height="1754" alt="convergence_C CG" src="https://github.com/user-attachments/assets/5f576ba0-18d9-4afe-bda4-5b6ed66a5526" />

<img width="2950" height="1754" alt="convergence_P cg" src="https://github.com/user-attachments/assets/96f857cf-6db9-4f83-8ac7-57ad37bf1190" />

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
