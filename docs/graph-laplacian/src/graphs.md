# Graphs and Graph Laplacians

Below, we discuss the key terminology and results when discussing graphs and
their Laplacians. We will also cite important results and demonstrate them on
simple examples.

## Directed Graphs

A **directed graph** is a pair $G = (V,E)$,
- where $V = \lbrace v_1,\ldots,v_m \rbrace$ is a collection of *nodes*,
- and $E \subset V \times V$ is a set of ordered pairs called *edges*.

Consider the following example of a graph $G = (V,E)$.

!!! danger

    INSERT EXAMPLE

Here,
- the five nodes in the graph are $V = \lbrace v_1,v_2,v_3,v_4,v_5 \rbrace$,
- and the edges are $E = \lbrace e_1=(v_1,v_2), e_2=(v_1,v_3), e_3=(v_3,v_2), e_4=(v_4,v_2), e_5=(v_2,v_5), e_6=(v_5,v_4), e_7=(v_3,v_4) \rbrace$.

### Properties of Edges
For a given edge in a directed graph, we are often interested in identifying the **source** node
and the **target** node. For example, edge $e_1$'s source node is $v_1$ and its
target node is $v_2$. For mathematical purposes, we define functions associated
with the source and target as follows.
- The source function, $s: (u,v) \mapsto u$. For instance, $s(e_3) = v_3$.
- The target function, $t: (u,v) \mapsto v$. For instance, $t(e_3) = v_2$.

We can summarize the edges of a directed graph in Algebraic form using a matrix:
the **incidence matrix** of a directed graph $G$ is a matrix $B \in \lbrace -1,
0, 1 \rbrace^{|V| \times |E|}$ where

$$B_{ij} = \begin{cases}
1 & s(e_j) = v_i \\
-1 & t(e_j) = v_i \\
0 & \textrm{otherwise}.
\end{cases}$$

In our example, the incidence matrix is

$$B = \begin{bmatrix}
1 & 1 & 0 & 0 & 0 & 0 & 0 \\
-1 & 0 & -1 & -1 & 1 & 0 & 0 \\
0 & -1 & 1 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 1 & 0 & -1 & -1 \\
0 & 0 & 0 & 0 & -1 & 1 & 0
\end{bmatrix}$$

### Properties of Nodes
For a given node in a directed graph, we are often interested in identifying
its degree. Then, the **degree of a node $v$** is the number of edges that
either originate from $v$ or terminate at $v$.

To formalize this notion in mathematical notation, first define the indicator
function for an event $\mathcal{A}$ as

$$\chi(\mathcal{A}) =
\begin{cases}
1 & \mathcal{A}~\textrm{is true}. \\
0 & \mathcal{A}~\textrm{is false}.
\end{cases}$$

The **degree of node $v$** is

$$d(v) = \sum_{ e \in E} \chi( s(e) = v ) + \chi(t(e) = v).$$

In our example, $d(v_2) = 4$.

!!! warning

    If we had edges that formed a loop (i.e., $s(e)=t(e) = v$) then the above
    definition of degree would count such an edge twice. Some authors allow for
    this while others preclude it.

The **degree matrix** of a graph is a diagonal matrix of dimension $|V| \times |V|$, where the entries on the diagonal are the degrees of the nodes. The degree matrix is denoted by $D(G)$.

From our example, ordered by the index of the nodes,

$$D(G) = \mathrm{diag}( 2, 4, 3, 3, 2).$$

## Undirected Graphs

!!! note

    Undirected graphs are often referred to simply as "graphs". We will adopt this
    nomenclature.

Similar to a directed graph, a **graph** is a pair $G=(V,E)$
- where $V = \lbrace v_1,\ldots,v_m \rbrace$ is a collection of *nodes*;
- and $E \subset V \times V$ is a set of *unordered* pairs called *edges*, where
the elements of the unordered pair are distinct.

Compare the definition of an edge for a graph to that of the directed graph.
One key difference is that the pair of nodes that form an edge are unordered
for a graph, whereas the ordering of the nodes in the edge is important in a
directed graph. A second key difference is that the pair of nodes that form
an edge for a graph cannot be the same, whereas we have not precluded this in
the case of directed graphs.

Consider the following example of a graph.

!!! danger

    Insert Example

Here,
- the five nodes of the graph are $V = \lbrace v_1,v_2,v_3,v_4,v_5 \rbrace$;
- and the edges are $E = \lbrace a=\lbrace v_1,v_2 \rbrace, b=\lbrace v_1,v_3\rbrace, c=\lbrace v_2,v_3\rbrace, d = \lbrace v_2, v_4\rbrace, e=\lbrace v_2,v_5\rbrace, f=\lbrace v_4, v_5 \rbrace, g=\lbrace v_3, v_4 \rbrace \rbrace$.

Compare this example to the one for directed graphs. We see that all of the
nodes are connected by the same edges, but the edges in the graph are not
directed as they are for the directed graph. As a result, when we list the
edges of the graph, the ordering of the nodes is not important. On the other
hand, when we list the edges of the directed graph, the ordering of the nodes is
important. We represent this distinction by using $\lbrace, \rbrace$ for graph,
and $(,)$ for a directed graph.

### Properties of Edges

For a graph, we can no longer define a source and a target for an edge as we did
for a directed graph. However, we can still define the **incidence matrix** as
a matrix $B \in \lbrace 0, 1 \rbrace^{|V| \times |E|}$, where

$$B_{ij} = \chi( v_i \in e_j ).$$

The incidence matrix for our example is

$$B = \begin{bmatrix}
1 & 1 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 1 & 1 & 1 & 0 & 0 \\
0 & 1 & 1 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 1 & 0 & 1 & 1 \\
0 & 0 & 0 & 0 & 1 & 1 & 0
\end{bmatrix}$$

Notice, the incidence matrix for the graph is just the entry-wise absolute value
of the incidence matrix for the directed graph. This is exactly what we should
expect as we have not changed how the nodes are connected. Later, we will see
that this underlies an integral property of Graph Laplacians.

For graphs, we also have the notion of a **path**. Roughly, a path from a node
$u$ to a node $v$ is a sequence of edges that where the initial edge contains
$u$ and the terminal edge in the sequence contains $v$, and there are no
repeated nodes. Formally, a **path** from $u$ to $v$ is a sequence of edges,
$(e_1,\ldots,e_k) \subset E$, where
- node $u$ is an element of $e_1$;
- node $v$ is an element of $e_k$;
- the set $e_{i-1} \cap e_i$ contains only one element (i.e., no reversals);
- and $e_{i} \cap e_j = \emptyset$ if $|i - j| > 1$ (i.e., no repeated nodes).

!!! note

    We can modify the definition to allow for closed paths from $u$ to $u$ by
    switching the third and fourth requirements by identifying $0$ with $k$.
    Therefore, $e_{0} = e_k$ and so we require $e_1 \cap e_k$ to have only
    one element, $u$. Similarly, $e_1 \cap e_k \neq \emptyset$ does not violate
    the fourth condition since $e_1 \cap e_k = e_1 \cap e_0$ and $|1 - 0| = 1$.

In our example, all of the paths from $v_1$ to $v_3$ are $(b)$, $(a,c)$,
$(a,d,g)$, and $(a,e,f,g)$. Note, $(a,e,f,d,c)$ is not a path as $v_2 \in a$
and $v_2 \in d$, which violates the fourth condition of a path. From this
example, we also see that paths can have different lengths. The **length** of a
path is the number of elements in its sequence of edges.

Using this notion of a path, we can define several other properties. We say
a graph is **connected** if, for any distinct $u,v \in V$, there is a path from
$u$ to $v$. Our example is a connected graph as there is a path from every
node to every other node.

A graph is a **tree** if it is connected but there are no closed paths (also
called *cycles*). Our example is not a tree as there are many cycles: for
example, $(a,c,b)$. However, our graph contains a **spanning tree**, which is a
subgraph of $G=(V,E)$, denoted by $T=(V,E')$, where the vertices are kept the
same but $E' \subset E$ such that $T$ is a tree. In our example, a spanning
tree can be generated in multiple ways. As examples,
- a spanning tree can be made with $E' = \lbrace a, c, d, e \rbrace$
- a spanning tree can be made with $E' = \lbrace b, g, d, f \rbrace$.

### Properties of Nodes

We will treat the nodes of the graph in the same way that we treated the nodes
of a directed graph. For any node, we can define its **degree** as the number of
edges that are connected to the node. In mathematical notation, the **degree of
node $v$** is

$$d(v) = | \lbrace e \in E : v \in e \rbrace |.$$

We can similarly construct **the degree matrix** to summarize the degrees of all
nodes in the graph. For our example, the degree matrix, ordered by the indices
of the nodes, is

$$D(G) = \mathrm{diag}( 2, 4, 3, 3, 2),$$

which is identical, unsurprisingly, to the degree matrix of the directed graph.


## Adjacency Matrix & Weighted Graphs

Up until now, we have presented a few ways of describing a (directed) graph that
have focused primarily on the relationships between edges and nodes. For
example, we introduced the incidence matrix, which describes which edges
correspond to which nodes (and, in the case of a directed graph, the direction
of these relationships). Similarly, we used the degree matrix to summarize the
number of edges that are connected to a node.

Our goal now is to introduce a matrix to describe the relationship between
nodes; specifically, this matrix will describe which nodes are neighbors (i.e.,
they are connected by a path of length $1$.). Define the **adjacency matrix** of
a (directed) graph $G=(V,E)$ by the matrix $A$ where

$$A_{ij} = \begin{cases}
1 & \exists e \in E: ~ v_i,v_j \not\in e, ~ i \neq j, \\
0 & \textrm{otherwise}.
\end{cases}$$

For the directed and undirected graphs shown above,

$$A = \begin{bmatrix}
0 & 1 & 1 & 0 & 0 \\
1 & 0 & 1 & 1 & 1 \\
1 & 1 & 0 & 1 & 0 \\
0 & 1 & 1 & 0 & 1 \\
0 & 1 & 0 & 1 & 0
\end{bmatrix}.$$

We see that the adjacency matrix fully describes a graph $G$, as we know how
many nodes there are and to which nodes they are connected. In our example,
the value $A_{12} = 1$, which implies that $v_1$ is connected by an edge to
node $v_2$. Unfortunately, for a directed graph $G$, $A$ loses the direction
associated with the edges of a directed graph.

The adjacency matrix also gives us a starting point for generalizing the graphs
that we have been discussing to the case where the edges have some weight, where
the weights represent the relative strength of the connection between two nodes
in the graph that share an edge.

To be specific, a **weighted graph** is a pair $G=(V,W)$,
- where $V = \lbrace v_1,\ldots,v_m \rbrace$ is a collection of *nodes*;
- and $W \in \mathbb{R}_{\geq 0}^{|V| \times |V|}$, where $W_{ii} = 0$ for  all $i=1,\ldots,m$ and $W$ is symmetric, is called the *weight matrix*.

Notice that $W$ generalizes the role of the adjacency matrix. Moreover, notice
that in a weighted graph, if $W_{ij} = 0$, then there is not an edge
between $v_i$ and $v_j$. On the otherhand, if $W_{ij} > 0$, then there is an
edge between $v_i$ and $v_j$.

Using this definition of a weighted graph, we can also define the **degree of a
node $v_i$** by

$$d(v_i) = \sum_{j=1}^m W_{ij},$$

which no longer needs to be an integer for a weighted graph as it was for a
directed or undirected graph. We can also define the **degree matrix** of a
weighted graph to be the matrix whose diagonal entries are the degrees of each
node.

## Graph Laplacian

It turns our that the difference of the degree matrix and the adjacency matrix
(or weight matrix) often shows up in applications as the coefficient matrix of a
linear system that needs to be solved, or as tool to describe the connectivity
of a graph, say, for clustering nodes.

Formally, this difference, called the **Graph Laplacian** of a (directed or
undirected) graphed $G=(V,E)$, is

$$L = D - A,$$

- where $D$ is the degree matrix of graph $G$;
- and $A$ is the adjacency matrix of graph $G$.

Analogously, the **Graph Laplacian** of a weighted graph, $G=(V,W)$, is

$$L = D - W,$$

where $D$ is the degree matrix of the weighted graph.

In the case of a directed graph, the Graph Laplacian is also related to the
incidence matrix, $B$, by the relationship

$$L = B B'.$$

Interestingly, it turns out that if the orientations of a directed graph are
arbitrarily changed (for example, flipping the direction of $e_2$ and $e_5$),
the *Graph Laplacian remains unchanged.* In other words, the Graph Laplacian
reflects the underlying connectivity of the graph and is invariant to the
orientation. We can use this idea to specify the corresponding relationship
between the Graph Laplacian of a (weighted or undirected) graph and a variant
of said graph's incidence matrix.

For a graph $G = (V,E)$, let $G_\sigma$ denote the graph with the same set of
nodes, but whose edges are now arbitrarily oriented according to some rule
$\sigma$ that takes in an edge $e \in E$ and specifies an oriented edge $e'$,
where the elements of $e$ and $e'$ are the same. Then, we can specify an
incidence matrix $B(G_\sigma)$ corresponding to $G_\sigma$. Then,

$$L(G) = L(G_\sigma) = B(G_\sigma) B(G_\sigma)'.$$

Similarly, for a weighted graph $G = (V,W)$, we can

1. Extract the underlying (unweighted, undirected) graph $\bar G$
2. Using the underlying graph, specify an orientation function and create an underlying directed graph, $\bar G_\sigma$.
3. Define the incidence matrix for any node $v_i$ in $\bar G_{\sigma}$ and edge $\bar e_j$ in $\bar G_{\sigma}$ by

$$B_{ij}(\bar G_{\sigma}) = \begin{cases}
\sqrt{W_{ij}} & s( \bar e_j ) = v_i, \\
-\sqrt(W_{ij}) & t( \bar e_j) = v_i, \\
0 & \textrm{otherwise.}
\end{cases}$$

Then,
$$L(G) = L(\bar G_\sigma) = B( \bar G_\sigma) B( \bar G_\sigma)'.$$

!!! danger

    Insert Example 
