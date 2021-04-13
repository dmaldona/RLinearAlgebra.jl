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

Notice, $L_1$ is a lower triangular matrix and $A_1 := B_{22} - B_{21}B_{11}^{-1} B_{12}$ is SPSD by **Fact 1** above.

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
