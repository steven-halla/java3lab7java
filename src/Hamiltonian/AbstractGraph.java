package Hamiltonian;
import java.util.*;
/**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*<pre>
 * Class        AbstractGraph.java
 * Description  Crucial abstract class in the triad of defining data structures
 *              in Java with interfaces, abstract classes and concrete classes.
 *              It implements the Graph interface which in turn contains all 
 *              the common operations of graphs. It overrides all the methods
 *              in the interface Graph and it contains an inner class Tree, very
 *              differnt from a Tree interface defined in Chapter 25 of the text.
 *              It is a parent class for the UnweightedGraph and WeightedGraph
 *              concrete classes.
 * Platform     jdk 1.8.0_241; NetBeans IDE 11.3; PC Windows 10
 * Course       CS 143
 * Hourse       1 hours and 9 minutes
 * Date         4/5/2021
 * History Log  7/18/2018, 5/7/2020
 * @author	<i>Niko Culevski</i>
 * @version 	%1% %2%
 * @param       <V> generic type
 *</pre>
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
public abstract class AbstractGraph<V> implements Graph<V> 
{
    protected List<V> vertices = new ArrayList<>(); // Store vertices
    protected List<List<Edge>> neighbors = new ArrayList<>(); // Adjacency lists

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Constructor  AbstractGraph()-default constructor
     * Description  Not used
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    protected AbstractGraph() 
    {
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Constructor  AbstractGraph()-overloaded constructor
     * Description  Construct a graph from vertices and edges stored in arrays.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       vertices V{}
     * @param       edges int[][]
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    protected AbstractGraph(V[] vertices, int[][] edges) 
    {
        for (V vertex : vertices)
        {
            addVertex(vertex);
        }
        createAdjacencyLists(edges, vertices.length);
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Constructor  AbstractGraph()-overloaded constructor
     * Description  Construct a graph from vertices and edges stored in List.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       vertices List
     * @param       edges List
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    protected AbstractGraph(List<V> vertices, List<Edge> edges) 
    {
        for (int i = 0; i < vertices.size(); i++)
            addVertex(vertices.get(i));

        createAdjacencyLists(edges, vertices.size());
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Constructor  AbstractGraph()-overloaded constructor
     * Description  Construct a graph for integer vertices 0, 1, 2 and edge list.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       edges List
     * @param       numberOfVertices int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    protected AbstractGraph(List<Edge> edges, int numberOfVertices) 
    {
        for (int i = 0; i < numberOfVertices; i++) 
            addVertex((V)(new Integer(i))); // vertices is {0, 1, ...}

        createAdjacencyLists(edges, numberOfVertices);
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Constructor  AbstractGraph()-overloaded constructor
     * Description  Construct a graph from integer vertices 0, 1,... and edge array.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       edges int[][]
     * @param       numberOfVertices int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    protected AbstractGraph(int[][] edges, int numberOfVertices) 
    {
        for (int i = 0; i < numberOfVertices; i++) 
            addVertex((V)(new Integer(i))); // vertices is {0, 1, ...}

        createAdjacencyLists(edges, numberOfVertices);
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       createAdjacencyLists
     * Description  Create adjacency lists for each vertex.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       edges int[][]
     * @param       numberOfVertices int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/    
    private void createAdjacencyLists(int[][] edges, int numberOfVertices) 
    {
        for (int i = 0; i < edges.length; i++) 
        {
            addEdge(edges[i][0], edges[i][1]);
        }
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       createAdjacencyLists
     * Description  Create adjacency lists for each vertex.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       edges List<Edges>
     * @param       numberOfVertices int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/   
    private void createAdjacencyLists(List<Edge> edges, int numberOfVertices) 
    {
        for (Edge edge: edges) 
        {
            addEdge(edge.u, edge.v);
        }
    }

    @Override
    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getSize
     * Description  Overridden method to return number of vertices in the graph.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      size int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    public int getSize() 
    {
      return vertices.size();
    }

    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getVertices
     * Description  Overridden method to return number of vertices as a List 
     *              in the graph.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      list List<V>
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    public List<V> getVertices() 
    {
      return vertices;
    }

    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getVertex
     * Description  Return the object for the specified vertex.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      vertex V
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/    
    public V getVertex(int index) 
    {
        return vertices.get(index);
    }

    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getIndex
     * Description  Return the index for the specified vertex object.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       vertex V
     * @return      index int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/    
    public int getIndex(V v) 
    {
        return vertices.indexOf(v);
    }

    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getNeighbors
     * Description  Return the neighbors of the specified vertex.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       index int
     * @return      neighbors List<Integer>
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/      
    public List<Integer> getNeighbors(int index) 
    {
        List<Integer> result = new ArrayList<>();
        for (Edge e: neighbors.get(index))
            result.add(e.v);

        return result;
    }

    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getDegree
     * Description  Return the degree for a specified vertex.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       v int
     * @return      degree int
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/    
    public int getDegree(int v) 
    {
        return neighbors.get(v).size();
    }

    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       printEdges
     * Description  Print the edges.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/     
    public void printEdges() 
    {
        for (int u = 0; u < neighbors.size(); u++) 
        {
            System.out.print(getVertex(u) + " (" + u + "): ");
            for (Edge e: neighbors.get(u)) 
            {
                System.out.print("(" + getVertex(e.u) + ", " +
                    getVertex(e.v) + ") ");
            }
            System.out.println();
        }
    }

     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       displayEdges
     * Description  Return the edges as a String.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      output String
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/     
    @Override
    public String displayEdges() 
    {
        StringBuilder output = new StringBuilder();
        for (int u = 0; u < neighbors.size(); u++) 
        {
            //output.append(getVertex(u) + " (" + u + "): ");
            output.append(getVertex(u) + " (" + u + "): ");
            for (Edge e: neighbors.get(u)) 
            {
                output.append("(" + getVertex(e.u) + ", " +
                    getVertex(e.v) + ") ");
            }
            output.append('\n');
        }
        return output.toString();
    }
    
    @Override
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       clear
     * Description  Clear graph.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/      
    public void clear() 
    {
        vertices.clear();
        neighbors.clear();
    }

    @Override 
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       addVertex
     * Description  Add a vertex to the graph.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       vertex V
     * @return      true/false boolean
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/      
    public boolean addVertex(V vertex) 
    {
        if (!vertices.contains(vertex)) 
        {
            vertices.add(vertex);
            neighbors.add(new ArrayList<Edge>());
            return true;
        }
        else 
        {
            return false;
        }
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       addEdge
     * Description  Add an edge to the graph.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       e Edge
     * @return      true/false boolean
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/     
    protected boolean addEdge(Edge e) 
    {
        if (e.u < 0 || e.u > getSize() - 1)
            throw new IllegalArgumentException("No such index: " + e.u);

        if (e.v < 0 || e.v > getSize() - 1)
            throw new IllegalArgumentException("No such index: " + e.v);

        if (!neighbors.get(e.u).contains(e)) 
        {
            neighbors.get(e.u).add(e);
            return true;
        }
        else 
        {
            return false;
        }
    }

    @Override
    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       addEdge
     * Description  Add an edge to the graph.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       u int
     * @param       v int
     * @return      true/false boolean
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/  
    public boolean addEdge(int u, int v) 
    {
        return addEdge(new Edge(u, v));
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Class        Edge
     * Description  Edge nested class inside the AbstractGraph class.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
    *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/    
    public static class Edge 
    {
        public int u; // Starting vertex of the edge
        public int v; // Ending vertex of the edge

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Constructor   Edge()-default constructor
        * Description   Construct an edge for (u, v).
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @param        u int
        * @param        v int
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
        public Edge(int u, int v) 
        {
            this.u = u;
            this.v = v;
        }
        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         *<pre>
         * Method       equals()
         * Description  Overridden method to check equality between edges.
         * @return      true or flase boolean
         * @param       obj Object
         * @author      <i>Niko Culevski</i>
         * Date         5/10/2020
         * History Log  7/18/2018  
         *</pre>
        *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        public boolean equals(Object obj) 
        {
            return u == ((Edge)obj).u && v == ((Edge)obj).v; 
        }
    }

    @Override
    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       dfs
     * Description  Depth-First Search method. Obtain a DFS tree starting from 
     *              vertex vCreates a parents int array and boolean array for 
     *              the visited certices. Calls recursive overloaded dfs method. 
     *              Returns a Tree object.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       v int
     * @return      tree Tree
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    public Tree dfs(int v) 
    {
        List<Integer> searchOrder = new ArrayList<>();
        int[] parent = new int[vertices.size()];
        for (int i = 0; i < parent.length; i++)
            parent[i] = -1; // Initialize parent[i] to -1

        // Mark visited vertices
        boolean[] isVisited = new boolean[vertices.size()];

        // Recursively search
        dfs(v, parent, searchOrder, isVisited);

        // Return a search tree
        return new Tree(v, parent, searchOrder);
    }

    /** Recursive method for DFS search */
    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       dfs
     * Description  Overloaded Depth-First Search method.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       u int
     * @param       parent int[]
     * @param       searchOrder List<Integer>
     * @param       isVisited boolean[]
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    private void dfs(int u, int[] parent, List<Integer> searchOrder,
        boolean[] isVisited) 
    {
        // Store the visited vertex
        searchOrder.add(u);
        isVisited[u] = true; // Vertex v visited

        for (Edge e : neighbors.get(u)) 
        {
            if (!isVisited[e.v]) 
            {
                parent[e.v] = u; // The parent of vertex e.v is u
                dfs(e.v, parent, searchOrder, isVisited); // Recursive search
            }
        }
    }

    @Override /** Starting bfs search from vertex v */
    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       dfs
     * Description  Breadth-First Search method. Obtain a DFS tree starting from 
     *              vertex v. Creates a parents int array and boolean array for 
     *              the visited certices. Returns a Tree object.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       v int
     * @return      tree Tree
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/     
    public Tree bfs(int v) 
    {
        List<Integer> searchOrder = new ArrayList<>();
        int[] parent = new int[vertices.size()];
        for (int i = 0; i < parent.length; i++)
          parent[i] = -1; // Initialize parent[i] to -1

        java.util.LinkedList<Integer> queue =
            new java.util.LinkedList<>(); // list used as a queue
        boolean[] isVisited = new boolean[vertices.size()];
        queue.offer(v); // Enqueue v
        isVisited[v] = true; // Mark it visited

        while (!queue.isEmpty()) 
        {
            int u = queue.poll(); // Dequeue to u
            searchOrder.add(u); // u searched
            for (Edge e: neighbors.get(u)) 
            {
                if (!isVisited[e.v]) 
                {
                    queue.offer(e.v); // Enqueue w
                    parent[e.v] = u; // The parent of w is u
                    isVisited[e.v] = true; // Mark it visited
                }
            }
        }

        return new Tree(v, parent, searchOrder);
    }

    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Class        Tree
     * Description  Tree inner class inside the AbstractGraph class.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    public class Tree 
    {
        private int root; // The root of the tree
        private int[] parent; // Store the parent of each vertex
        private List<Integer> searchOrder; // Store the search order

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Constructor   Tree()
        * Description   Construct a tree with root, parent, and searchOrder.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @param        root int
        * @param        parent int[]
        * @param        searchOrder List
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
        public Tree(int root, int[] parent, List<Integer> searchOrder) 
        {
            this.root = root;
            this.parent = parent;
            this.searchOrder = searchOrder;
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method       getRoot
        * Description  Return the root of the tree.
        * Date         4/5/2021
        * History Log  7/18/2018, 5/7/2020
        * @author      <i>Niko Culevski</i>
        * @return      root int
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/  
        public int getRoot() 
        {
            return root;
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method        getParent
        * Description   Return the parent of vertex v.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @param        v int
        * @return       root int
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/          
        public int getParent(int v) 
        {
            return parent[v];
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method        getSearchOrder
        * Description   Return a list representing search order.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @return       list List
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/        
        public List<Integer> getSearchOrder() 
        {
            return searchOrder;
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method        getNumberOfVerticesFound
        * Description   Return number of vertices found.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @return       vertices int
       *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/          
        public int getNumberOfVerticesFound() 
        {
            return searchOrder.size();
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method        getPath
        * Description   Return the path of vertices from a vertex to the root.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @param        index int
        * @return       vertices List
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/          
        public List<V> getPath(int index) 
        {
            ArrayList<V> path = new ArrayList<>();

            do 
            {
                path.add(vertices.get(index));
                index = parent[index];
            }
            while (index != -1);

            return path;
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method        printPath
        * Description   Print a path from the root to vertex v.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        * @param        index int
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/         
        public void printPath(int index) 
        {
            List<V> path = getPath(index);
            System.out.print("A path from " + vertices.get(root) + " to " +
                vertices.get(index) + ": ");
            for (int i = path.size() - 1; i >= 0; i--)
                System.out.print(path.get(i) + " ");
        }

        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        *<pre>
        * Method        printTree
        * Description   Print the whole tree.
        * Date          4/5/2021
        * History Log   7/18/2018, 5/7/2020
        * @author       <i>Niko Culevski</i>
        *</pre>
       *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/        
        public void printTree() 
        {
            System.out.println("Root is: " + vertices.get(root));
            System.out.print("Edges: ");
            for (int i = 0; i < parent.length; i++) 
            {
                if (parent[i] != -1) 
                {
                    // Display an edge
                    System.out.print("(" + vertices.get(parent[i]) + ", " +
                        vertices.get(i) + ") ");
                }
            }
            System.out.println();
        }
    }
    
    //Addition to the OG Abstrat Graph From lab 6
    
            /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getHamiltornianPath
     * Description  REturns a Hamiltonian path from the specificed vertex label.
     *              Return a null if the graph does not contain a hamiltonian path.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       vertex V
     * @return      list List
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    public List<Integer> getHamiltonianPath(V vertex) {
        return getHamiltonianPath(getIndex(vertex));
    }
    
                /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getHamiltonianPath
     * Description  REturns a Hamiltonian path from the specificed vertex label.
     *              Return a null if the graph does not contain a hamiltonian path.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       vertex V
     * @return      list List
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    public List<Integer> getHamiltonianPath(int v) {
        int[] next = new int[getSize()];
        for (int i = 0; i > next.length; i++)
            next[i] = -1;
        
        boolean[] isVisited = new boolean[getSize()];
        
        List<Integer> result = null;
        
        for (int i = 0; i < getSize(); i++)
            reorderNeigborsBasedOnDegree(getNeighbors(i));
        
        if (getHamiltonianPath(v, next, isVisited)) {
            result = new ArrayList<Integer>();
            int vertex = v;
            while (vertex != -1) {
                result.add(vertex);
                vertex = next[vertex];
            }
        }
        return result;
    }
    
     /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       reorderNeighborsOnDegree
     * Description  Records the adjacency list in increasing order of degrees
     *              at 37:50
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      list List
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    private void reorderNeigborsBasedOnDegree(List<Integer>list) {
        for (int i = list.size() - 1; i >= 1; i--)
        {
            int currentMaxDegree = getDegree(list.get(0));
            int currentMaxIndex = 0;
            
            for (int j = 1; j <= i; j++)
            {
                if (currentMaxDegree < getDegree(list.get(j)))
                {
                    currentMaxDegree = getDegree(list.get(j));
                    currentMaxIndex = j;
                }
        }
            if (currentMaxIndex != 1)
            {
                int temp = list.get(currentMaxIndex);
                list.set(currentMaxIndex, list.get(i));
                list.set(i, temp);
            }
    }
    }
         /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       allVisited
     * Description  Returns true if all elements in array isVisited are true
     *              at 37:50
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      true/false boolean[]
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    private boolean allVisited(boolean[] isVisited)
    {
        boolean result = true;
        
        for (int i = 0; i < getSize(); i++)
            result = result && isVisited[i];
        
        return result;
    }
    
    
        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getHamiltonianPath
     * Description  search for a Hamiltonian path from v
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       isVisited boolean[][
     * @return      true/false boolean
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    private boolean getHamiltonianPath(int v, int[] next,
           boolean[] isVisited) {
        isVisited[v] = true;
        
        if (allVisited(isVisited))
            return true;
        
        for (int i = 0; i < getNeighbors(v).size(); i++)
        {
            int u = getNeighbors(v).get(i);
            if (!isVisited[u] && getHamiltonianPath(u, next, isVisited))
            {
                next[v] = u;
            return true;
        }
    }
            isVisited[v] = false;
            System.out.println("Backtrack at " + v);
            return false;
}
            /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getHamiltonianPath
     * Description  search for a Hamiltonian path from v
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @return      list List
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    public List<Integer> getHamiltonianCycle()
    {
        return getHamiltonianCycle(0);
    }
    
    /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getHamiltonianCycle
     * Description  search for a Hamiltonian path from v
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       v int
     * @return      list List
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    public List<Integer> getHamiltonianCycle(int v)
    {
        int[] next = new int[getSize()];
        for (int i = 0; i < next.length; i++)
            next[i] = 1;
        
        boolean[] isVisited = new boolean[getSize()];
        
        List<Integer> result = null;
        
        for (int i = 0; i < getSize(); i++)
        
            reorderNeigborsBasedOnDegree(getNeighbors(i));
        
        if (getHamiltonianCycle(v, next, isVisited))
        {
            result = new ArrayList<Integer>();
            int vertex = v;
            while (vertex != -1)
            {
                result.add(vertex);
                vertex = next[vertex];
            }
        }
            
            
        return result;
    }
    
        /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       getHamiltonianCycle
     * Description  return a hamiltonian cyle. return null if the graph does
     *              not contain a hamiltonian cycle.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       v int
     * @param       next int[]
     * @return      true/false boolean
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    public boolean getHamiltonianCycle(int v, int[] next, boolean[] isVisited)
    {
        isVisited[v] = true;
        
        if (allVisited(isVisited) && isCycle(v))
            return true;
        
        for (int i = 0; i < getNeighbors(v).size(); i++)
        {
            int u = getNeighbors(v).get(i);
            if (!isVisited[u] && getHamiltonianCycle(u, next, isVisited))
            {
                next[v] = u;
                return true;
            }
        }
        isVisited[v] = false;
        return false;
    }
    
            /**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *<pre>
     * Method       isCycle
     * Description  is is a cycle if the neighbors of v contan o.
     * Date         4/5/2021
     * History Log  7/18/2018, 5/7/2020
     * @author      <i>Niko Culevski</i>
     * @param       v int
  
     * @return      true/false boolean
     *</pre>
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    private boolean isCycle(int v)
    {
        return getNeighbors(v).contains(0);
    }
    
    
}
