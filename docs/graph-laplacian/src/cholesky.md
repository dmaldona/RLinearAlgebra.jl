# Cholesky Decompositions

Below, we consider decomposing a symmetric, positive semi-definite (SPSD) matrix
into the product that uses a lower triangular matrix and a diagonal matrix. We
also consider the corresponding block variant of this decomposition.

We will refer to this decomposition as the **Cholesky decomposition**. However,
we note that this is nonstandard. Specifically, this decomposition is a
generalization of the Cholesky decomposition, and is often referred to as the
**LDL decomposition**.

As mentioned, we will consider the scalar Cholesky decomposition and the block
Cholesky decomposition of spsd matrices. These decompositions are characterized
as follows. For the **Scalar Cholesky Decomposition** of a SPSD matrix $A$,
there exists a lower triangular matrix $L$ whose diagonal entries are ones and
a diagonal matrix $D$ whose entries are non-negative, such that

$$A = LDL',$$

where $L'$ denotes the transpose of $L$.

For the **Block Cholesky Decomposition** of a SPSD matrix $A$, there exists a
lower triangular matrix $L$ whose block diagonal entries are identity matrices
(possibly of different dimension), and a block diagonal matrix $D$ with SPSD
blocks, such that

$$A = LDL'.$$

There are several ways to construct the Cholesky decomposition and we refer to

- Golub, Gene H., and Charles F. Van Loan. Matrix computations. Vol. 3. JHU press, 2013.

for several derivations. Here, we will make use of the outer product (or Schur
Complement) derivation of the Cholesky decomposition.

!!! note

    This derivation has an equal number of operations in comparison to the more
    common derivation of the Cholesky decomposition, but often requires almost
    twice as many accesses to the entries of $A$. However, this formulation will
    be useful for our applications to Graph Laplacians.

## Preliminary Facts

To proceed in our derivation of the Cholesky Decompositions, we will need two
facts that can be readily verified, as we will describe below.

**Fact 1:** let

$$A = \begin{bmatrix}
C_{11} & C_{12} \\
C_{21} & C_{22}
\end{bmatrix},$$

be some partitioning of the matrix $A$. Note, $C_{21} = C_{12}'$. If $C_{11}$ is
positive definite, then the Schur Complement of $A$ with respect to $C_{11}$,

$$C_{22} - C_{21}C_{11}^{-1}C_{12},$$

is SPSD. The proof of this fact follows from first writing a quadratic equation
in terms of the partitioned components of $A$ and the corresponding partition of
an arbitrary vector $v$. Since $A$ is positive semi-definite, this quadratic
equation is non-negative for any choice of $v$. In particular, if we minimize
over the component of $v$ that corresponds to $C_{11}$ and substitute back into
the quadratic for this choice, then we will conclude that the Schur complement
is positive semi-definite. Symmetry follows readily.

**Fact 2:** If $L_1, L_2$ are lower triangular matrices of the same dimension,
then $L_1L_2$ is also a lower triangular matrix. This fact can be established by
directly calculating the upper triangular components of $L_1L_2$.

## Scalar Cholesky Decomposition  

We begin by partitioning the matrix $A$ into

$$A = \begin{bmatrix}
B_{11} & B_{12} \\
B_{21} & B_{22}
\end{bmatrix},$$

where $B_{11}$ is a positive scalar. Note, if the upper left corner of $A$ is
zero (it cannot be negative since $A$ is SPSD), we can either pivot the entries
of $A$ or just start with a diagonal entry of $A$ that is positive. If there are
not any diagonal entries of $A$ that are positive, then $A$ must be zero and we
do not need to decompose it.

We will now construct the Cholesky decomposition of $A$ incrementally. Define

$$L_1 = \begin{bmatrix}
1 & \\
B_{21}B_{11}^{-1} & I
\end{bmatrix},$$

where the empty entry is a row of zeros and $I$ is the identity matrix that is
one dimension smaller than $A$. Also, define

$$D_1 = \begin{bmatrix}
B_{11} & \\
   & B_{22} - B_{21}B_{11}^{-1} B_{12}
\end{bmatrix},$$

where the empty entries are zeros. Then, we can verify by direct calculation,

$$A = L_1 D_1 L_1'.$$

Notice, $L_1$ is a lower triangular matrix and $B_{22} - B_{21}B_{11}^{-1} B_{12}$ is SPSD by **Fact 1** above.

We can now proceed by repeating our argument to define a lower triangular $L_2$,
and block diagonal matrix $D_2$ whose first block is $1 \times 1$ such that $A_1$
can be written as $L_2 D_2 L_2'$. Putting the two steps together, we have

$$A = L_1 \begin{bmatrix}
B_{11} & \\
    & L_2 D_2 L_2'
\end{bmatrix}L_1' =
L_1 \begin{bmatrix}
1 & \\
 & L_2
\end{bmatrix}
\begin{bmatrix}
B_{11} & \\
    & D_2
\end{bmatrix}
\begin{bmatrix}
1 & \\
    & L_2'
\end{bmatrix}
L_1'.$$

By **Fact 2**, the product of $L_1$ with

$$\begin{bmatrix}
1 & \\
& L_2
\end{bmatrix}$$

is lower triangular as both terms in the product are lower triangular. We now
repeat this reasoning on the second block of $D_2$, and so on until we are left
with a diagonal matrix. This diagonal matrix is the $D$ in the Cholesky
Decomposition. The product of the lower triangular matrices generate the $L$ in
the Cholesky Decomposition.

## Example of Scalar Cholesky Decomposition

Consider decomposing the matrix,

$$A = \begin{bmatrix}
2 & -1 & -1 & 0 & 0 \\
-1 & 4 & -1 & -1 & -1 \\
-1 & -1 & 3 & -1 & 0 \\
0 & -1 & -1 & 3 & -1 \\
0 & -1 & 0 & -1 & 2
\end{bmatrix}.$$

Starting from the top left corner and following our procedure,

$$L_1 = \begin{bmatrix}
1 & \\
-1/2 & 1 \\
-1/2 & & 1 \\
 0 & & & 1 \\
 0 & & & & 1
\end{bmatrix},$$

and

$$\begin{aligned}
D_1 &= \begin{bmatrix}
2 &  &  &  &  \\
 & 4 & -1 & -1 & -1 \\
 & -1 & 3 & -1 & 0 \\
 & -1 & -1 & 3 & -1 \\
 & -1 & 0 & -1 & 2
\end{bmatrix} - \begin{bmatrix}
 &  &  &  &  \\
 & -1/2 & -1/2 & 0 & 0 \\
 & -1/2 & -1/2 & 0 & 0 \\
 & 0 & 0 & 0 & 0 \\
 & 0 & 0 & 0 & 0
\end{bmatrix}, \\
&= \begin{bmatrix}
2 &  &  &  &  \\
 & 7/2 & -3/2 & -1 & -1 \\
 & -3/2 & 5/2 & -1 & 0 \\
 & -1 & -1 & 3 & -1 \\
 & -1 & 0 & -1 & 2
\end{bmatrix}
\end{aligned}$$

where missing entries are zero.

Repeating the procedure for the block component of $D_1$, we calculate

$$L_2 = \begin{bmatrix}
1 & \\
-3/7 & 1 \\
-2/7 & & 1 \\
-2/7 & & & 1,
\end{bmatrix}$$

and

$$\begin{aligned}
D_2 &= \begin{bmatrix}
7/2 & & & \\
&  5/2 & -1 & 0 \\
&  -1 & 3 & -1 \\
&  0 & -1 & 2
\end{bmatrix} -
\begin{bmatrix}
0 & & & \\
  & 9/14 & 3/7 & 3/7 \\
  & 3/7 & 2/7 & 2/7 \\
  & 3/7 & 2/7 & 2/7
\end{bmatrix} \\
&=\begin{bmatrix}
7/2 & & & \\
&  13/7 & -10/7 & -3/7 \\
&  -10/7 & 19/7 & -9/7 \\
&  -3/7 & -9/7 & 12/7
\end{bmatrix}
\end{aligned}$$

Repeating the procedure for the block component of $D_2$, we calculate

$$L_3 = \begin{bmatrix}
1 & \\
-10/13 & 1 \\
-3/13 & & 1
\end{bmatrix},$$

and

$$\begin{aligned}
D_3 &= \begin{bmatrix}
13/7 & & \\
 & 19/7 & -9/7 \\
 & -9/7 & 12/7
\end{bmatrix} - \begin{bmatrix}
0 & & \\
  & 100/91 & 30/91 \\
  & 30/91 & 9/91
\end{bmatrix} \\
&= \begin{bmatrix}
13/7 & &  \\
&147/91 & -147/91 \\
&-147/91 & 147/91
\end{bmatrix}.
\end{aligned}$$

For the last block diagonal of $D_3$, we conclude

$$L_4 = \begin{bmatrix}
1 & \\
-1 & 1
\end{bmatrix},$$

and

$$\begin{aligned}
D_4 &= \begin{bmatrix}
147/91 &  \\
 & 147/91 - 147/91
\end{bmatrix} = \begin{bmatrix}
147/91 & \\
 & 0
\end{bmatrix}.
\end{aligned}$$

Finally, let $I_k$ denote the $k \times k$ identity matrix. We can compute

$$\begin{aligned}
L &= L_1 \begin{bmatrix}
1 & \\
  & L_2
\end{bmatrix} \begin{bmatrix}
I_2 & \\
  & L_3
\end{bmatrix} \begin{bmatrix}
I_3 & \\
  & L_4
\end{bmatrix} \\
&= \begin{bmatrix}
1 & \\
-1/2 & 1 \\
-1/2 & - 3/7 & 1 \\
0 & -2/7 & -10/13 & 1 \\
0 & -2/7 & -3/13 & -1 & 1
\end{bmatrix}.
\end{aligned}$$

Moreover, from $D_1,\ldots,D_4$, we can recover
$$D = \begin{bmatrix}
2 & \\
  & 7/2 \\
  & & 13/7 \\
  & & & 147/91 \\
  & & & & 0
\end{bmatrix}.$$

We can verify by direct calculation that $A = LDL'$.

## Block Cholesky Decomposition

We will see that the block Cholesky Decomposition can be derived almost
identically to the scalar Cholesky Decomposition. Again, we begin by
partitioning the matrix $A$ into

$$A = \begin{bmatrix}
B_{11} & B_{12} \\
B_{21} & B_{22}
\end{bmatrix},$$

where $B_{11}$ is positive definite matrix of some dimension $d$ (smaller than
the rank of $A$). Note, it is difficult to know beforehand whether $B_{11}$ is
positive definite, so we need to know something about the blocks of our matrix
before proceeding.

We will now construct the Cholesky decomposition of $A$ incrementally, just
as we did before. Let $n$ denote the number of rows in $A$. Define

$$L_1 = \begin{bmatrix}
I_d & \\
B_{21}B_{11}^{-1} & I_{n-d}
\end{bmatrix},$$

where the empty entry is a row of zeros and $I_k$ is the $k \times k$ identity
matrix. Also, define

$$D_1 = \begin{bmatrix}
B_{11} & \\
   & B_{22} - B_{21}B_{11}^{-1} B_{12}
\end{bmatrix},$$

where the empty entries are zeros. Then, we can verify by direct calculation,

$$A = L_1 D_1 L_1'.$$

Notice, $L_1$ is a lower triangular matrix with block identity matrices on its
diagonal. Moreover, $B_{22} - B_{21}B_{11}^{-1} B_{12}$ is SPSD by **Fact 1**
above.

We can now proceed by partitioning $B_{22} - B_{21}B_{11}^{-1} B_{12}$ just as
we did for $A$, and computing $L_2$ and $D_2$ so that $L_2$ is a lower
triangular matrix of dimension $n-d \times n-d$ with identity matrices on its
block diagonal, and so that $D_2$ is block diagonal. We can then compute $L$
and $D$ just as we did in the case of the scalar Cholesky Decomposition.

## Example of Block Cholesky Decomposition

Consider decomposing the matrix,

$$A = \begin{bmatrix}
2 & -1 & -1 & 0 & 0 \\
-1 & 4 & -1 & -1 & -1 \\
-1 & -1 & 3 & -1 & 0 \\
0 & -1 & -1 & 3 & -1 \\
0 & -1 & 0 & -1 & 2
\end{bmatrix}.$$

We will work with $2 \times 2$ blocks until we have exhausted all such blocks.
From our derivation,

$$L_1 = \begin{bmatrix}
1 & \\
  & 1 \\
-5/7 & -3/7 & 1 \\
-1/7 & -2/7 &  & 1 \\
-1/7 & -2/7 &  & & 1
\end{bmatrix},$$

and

$$\begin{aligned}
D_1 &= \begin{bmatrix}
2 & -1 & & & \\
-1 & 4 & & & \\
 & & 3 & -1 & 0 \\
 & & -1 & 3 & -1 \\
 & & 0 & -1 & 2
\end{bmatrix} - \begin{bmatrix}
& & & & & \\
& & & & & \\
& & 8/7 & 3/7 & 3/7 \\
& & 3/7 & 2/7 & 2/7 \\
& & 3/7 & 2/7 & 2/7
\end{bmatrix} \\
&= \begin{bmatrix}
2 & -1 & & & \\
-1 & 4 & & & \\
 & & 13/7 & -10/7 & -3/7 \\
 & & -10/7 & 19/7 & -9/7 \\
 & & -3/7 & -9/7 & 12/7
\end{bmatrix}
\end{aligned}$$

We now repeat the same exercise on the upper left $2 \times 2$ block matrix in
the $3 \times 3$ block of $D_1$. From our derivation,

$$L_2 = \begin{bmatrix}
-1 & -1
\end{bmatrix},
$$

and

$$\begin{aligned}
D_2 &= \begin{bmatrix}
13/7 & -10/7 &  \\
-10/7 & 19/7 &  \\
 & & 0
\end{bmatrix}
\end{aligned}$$

Now, computing
$$\begin{aligned}
L = L_1 \begin{bmatrix}
1 \\
0 & 1 \\
0 & 0 & 1 \\
0 & 0 & 0 & 1 \\  
0 & 0 & -1 & -1 & 1
\end{bmatrix} = \begin{bmatrix}
1 & \\
  & 1 \\
-5/7 & -3/7 & 1 \\
-1/7 & -2/7 &  & 1 \\
-1/7 & -2/7 &  -1 & -1 & 1
\end{bmatrix}
\end{aligned}$$

and letting

$$D = \begin{bmatrix}
2 & -1 & & & \\
-1 & 4 & & & \\
& & 13/7 & -10/7 &  \\
& & -10/7 & 19/7 &  \\
& & & & 0
\end{bmatrix},$$

we have $A = LDL'$.
