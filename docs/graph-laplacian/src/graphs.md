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
- A spanning tree can be made with $E' = \lbrace b, g, d, f \rbrace$.

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


## Graph Laplacian
